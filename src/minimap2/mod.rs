//! The primary `mini_consensus::Resolver` struct is the consensus resolver 
//! that uses rammap/minimap2 to locate identical spans between multiple input 
//! sequences and anchored POA to resolve conflicts in the interval regions.

// modules
mod add_aln;
mod add_seq;

// imports
use std::iter::repeat_n;
use rammap::{Aligner, Preset};
use crate::poa::{PoaConfig, AlignmentMode, Poa};

// data type aliases
type BaseByte  = u8;
type SeqPos0 = usize;

/* -----------------------------------------------------------------------------
mini_consensis::Resolver - scaffolded multiple sequence alignment and consensus
----------------------------------------------------------------------------- */
/// A reusable minimap2-assisted POA consensus resolver.
pub struct Resolver {
    aligner_preset: Preset,
    base_capacity: usize,
    anchor_len: usize,
    poa: Poa,
    aligner: Option<Aligner>,
    seqs: Vec<Vec<BaseByte>>,
    scaffold: Vec<BaseByte>,
    is_identical: Vec<bool>, // map of scaffold whether all seqs are identical at this position
    seq_maps: Vec<Vec<SeqPos0>>, // outer=seq, inner=scaffold SeqPos0, value=seq SeqPos0
    scaffold_pos0: SeqPos0,
    seq_pos0: SeqPos0,
}
impl Resolver {

    /// Create a new reusable `Resolver` consensus engine with the indicated 
    /// sequence and sequence length capacity. The returned engine can be used 
    /// for all  consensus building that uses the same configuration parameters.  
    /// Capacity will grow and persist stably as needed if the number of 
    /// sequences or their lengths exceed the values provided here.
    /// 
    /// `poa_bandwidth` determines how much of the POA matrix is filled, with a 
    /// default of . Use 
    /// a higher number
    /// 
    /// `anchor_len` is the required number of bases identical among all
    /// sequences on each side of each POA-resolved span (default: 5 bases, thus
    /// identity spans less than 10 bases will be resolved by POA along with 
    /// the flanking variant spans).
    pub fn with_capacity(
        aligner_preset: Preset,
        n_seqs:  usize,
        n_bases: usize,
        anchor_len:    Option<usize>,
    ) -> Self {
        const DEFAULT_ANCHOR_LEN: usize = 5;
        let poa_config: PoaConfig = PoaConfig {
            band_width: 0, // instead, we update band_width based on seq len differences belows
            adaptive_band: false,
            alignment_mode: AlignmentMode::Global, // since all POA is anchored
            ..PoaConfig::default()
        };
        let seq_map = Vec::with_capacity(n_bases);
        Self {
            aligner_preset,
            base_capacity: n_bases,
            anchor_len: anchor_len.unwrap_or(DEFAULT_ANCHOR_LEN),
            poa: Poa::with_capacity(poa_config, n_seqs, 500),
            aligner: None,
            seqs:         Vec::with_capacity(n_seqs),
            scaffold:     Vec::with_capacity(n_bases),
            is_identical: Vec::with_capacity(n_bases),
            seq_maps: vec![seq_map; n_seqs],
            scaffold_pos0: 0,
            seq_pos0: 0,
        }
    }

    /// Reset a `Resolver` to begin a new consensus event by adding a new
    /// scaffold sequence. When initialized using `resolver.set_scaffold()`,
    /// all sequences must be pre-aligned to scaffold and added using 
    /// `resolver.add_aln()`. Use `set_scaffold_with_aligner()` if you will add
    /// sequences using `resolver.add_seq()`.
    pub fn set_scaffold(
        &mut self,
        scaffold: &[BaseByte],
    ) {
        self.aligner = None;
        self.seqs.clear();
        self.scaffold.clear();
        self.scaffold.extend_from_slice(scaffold);
        self.is_identical.clear();
        self.is_identical.extend(repeat_n(true, scaffold.len()));
    }

    /// Reset a `Resolver` to begin a new consensus event by adding a new
    /// scaffold sequence, where at least some sequences will be added using 
    /// `resolver.add_seq()` and thus require a minimap2 aligner.
    pub fn set_scaffold_with_aligner(
        &mut self,
        scaffold: &[BaseByte],
    ) {
        self.set_scaffold(scaffold);
        // initialize the scaffold pairwise aligner
        // ensure that the CIGAR string is reported with =/X ops instead of M
        let mut aligner = Aligner::from_seqs(
            vec![("scaffold".to_string(), scaffold.to_vec())], 
            self.aligner_preset
        );
        aligner.output_config_mut().eqx = true;
        aligner.output_config_mut().do_cs = true;
        self.aligner = Some(aligner);
    }

    /// Use the minimap2-assisted scaffold maps to isolate anchored variant
    /// spans and resolve them by multi-sequence POA to return the consensus 
    /// sequence. Scaffold get a vote, and is the deciding vote when all 
    /// sequences are different at a position.
    pub fn get_consensus(
        &mut self,
    ) -> Option<Vec<BaseByte>> {
        let scaffold_len = self.scaffold.len();
        let n_seqs = self.seqs.len();

        let mut chunk_pos0:  usize = 0; // leftmost pos0 of the next encountered chunk in scaffold coordinates
        let mut left_start0: usize = 0; // leftmost pos0 of the uncommitted match span left of POA span
        let mut left_end1:   usize = 0; // righmost pos1 of the uncommitted match span left of POA span
        let mut poa_pending: bool  = false; // if true, left is set and we have a span pending POA
        let mut consensus: Vec<BaseByte> = Vec::new();

        for chunk in self.is_identical.chunk_by(|a, b| a == b){
            let is_identical = chunk[0];
            let n_chunk_pos = chunk.len();

            println!("{} {}", is_identical, n_chunk_pos);

            // in a span where all (or all but one) seqs matched scaffold
            if is_identical {

                // initialize the first (and possibly only) matching span
                if !poa_pending {
                    left_start0 = chunk_pos0;
                    left_end1 = chunk_pos0 + n_chunk_pos;

                // process a variant span by POA if sufficient anchors
                // too-short matching anchors are included in the POA span 
                } else if n_chunk_pos >= self.anchor_len * 2 {
                    let scaffold_start0 = left_end1.saturating_sub(self.anchor_len);
                    let scaffold_end1 = (chunk_pos0 + self.anchor_len).min(scaffold_len);

                    // fill out any left overhang of the seed sequence with its bases
                    if scaffold_start0 > left_start0 {
                        consensus.extend(&self.scaffold[left_start0..scaffold_start0]);
                    }

                    // collect the starts and ends of each sequence in its own coordinate space
                    // find the length difference between the longest and shortest sequence
                    let mut min_len = scaffold_end1 - scaffold_start0;
                    let mut max_len = min_len;
                    let seq_ranges: Vec<_> = (0..n_seqs).map(|seq0| {
                        let start0 =  self.seq_maps[seq0][scaffold_start0];
                        let end1 = self.seq_maps[seq0][scaffold_end1 - 1] + 1;
                        let len = end1 - start0;
                        if len < min_len { min_len = len }
                        if len > max_len { max_len = len }
                        (start0, end1)
                    }).collect();
                    let bandwidth = max_len - min_len + 1;

                println!("{}\t{}\t{}\t{}\t{}\t{}", 
                    left_end1 - left_start0,
                    chunk_pos0 - left_end1, 
                    min_len,
                    scaffold_end1 - scaffold_start0,
                    max_len,
                    bandwidth
                );

                    // seed the POA graph
                    self.poa.seed_new_graph(
                        &self.scaffold[scaffold_start0..scaffold_end1], 
                        Some(bandwidth)
                    );

                    // add/align all other sequences to the graph
                    for seq0 in 0..n_seqs {
                        let range = seq_ranges[seq0];
                        self.poa.add_seq(&self.seqs[seq0][range.0..range.1]);                        
                    }

                    // resolve and append the local consensus by heaviest bundle
                    consensus.extend(self.poa.get_heaviest_path());

                    // jump the left flank to the current right flank
                    // not including the flanking bases committed with POA above
                    left_start0 = scaffold_end1; 
                    left_end1 = chunk_pos0 + n_chunk_pos;
                    poa_pending = false;
                } 

            // in a span where at least one sequence differed from the scaffold
            } else {

                // handle flanking variant gaps, commit as scaffold bases
                if chunk_pos0 + n_chunk_pos == scaffold_len {
                    // do nothing, handled below after loop terminates
                } else if chunk_pos0 == 0 {
                    consensus.extend(&self.scaffold[0..n_chunk_pos]);
                
                // in middle, flag that we have a variant span pending a right anchor for POA
                } else {
                    poa_pending = true;
                }
            }
            chunk_pos0 += n_chunk_pos;
        }

        // fill out any right overhang of the seed sequence with its bases
        if left_start0 < scaffold_len {
            consensus.extend(&self.scaffold[left_start0..scaffold_len]);
        }

        // return our result
        Some(consensus)
    }
}
