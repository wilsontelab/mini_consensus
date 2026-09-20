//! The primary `mini_consensus::Resolver` struct is the consensus resolver 
//! that uses rammap/minimap2 to locate identical spans between multiple input 
//! sequences and anchored POA to resolve conflicts in the interval regions.

// modules
mod add_seq;

// imports
use rayon::prelude::*;
use rammap::{Aligner, Preset};
use crate::poa::{PoaConfig, AlignmentMode, Poa, PoaPool};

// data type aliases
type BaseByte  = u8;
type SeqPos0 = usize;
type SeqPos1 = usize;

/* -----------------------------------------------------------------------------
ResolverConfig
----------------------------------------------------------------------------- */
/// Configuration parameters for a consensus resolver. These parameters are used
/// as provided for all iterative calls to an instantiated resolver.
#[derive(Debug, Clone)]
pub struct ResolverConfig {
    /// A `rammap` (`minimap2`) alignment preset. 
    /// 
    /// Default is `MapHifi`.
    pub aligner_preset: Preset,

    /// The required number of identical bases on each side of a POA block to
    /// anchor local consensus resolution. 
    /// 
    /// Default is 5, thus, identity spans less than 10 bases will be resolved 
    /// by POA along with the flanking variant spans.
    pub anchor_len: usize,  

    /// Whether input sequences are all end-to-end on the scaffold, or, 
    /// alternatively, cover random spans of the scaffold.
    /// 
    /// End clips on sequences always contribute to the final consensus when
    /// `is_end_to_end` is `true`, and never do when it is `false`.
    /// 
    /// Default is `false`, i.e., consensuses are resolved in random coverage
    /// mode without using end clips.
    pub is_end_to_end: bool,
}
impl Default for ResolverConfig {
    /// Return a new `ResolverConfig` with default parameters that map only the 
    /// aligned poritions of HiFi reads to return an untrimmed consensus.
    fn default() -> Self {
        ResolverConfig {
            aligner_preset: Preset::MapHifi,
            anchor_len:     5,
            is_end_to_end:  false,
        }
    }
}

/* -----------------------------------------------------------------------------
Resolver - support structs and enums
----------------------------------------------------------------------------- */
/// A span, i.e., a region, on a chromosome.
pub struct ReferenceSpan {
    /// The BAM reference, i.e., target sequence id, as returned by noodles 
    /// `reference_sequence_id()`.
    pub reference_sequence_id: usize,

    /// 1-based inclusive span start, as returned by noodles `alignment_start()`.
    pub start1: SeqPos1,

    /// 1-based inclusive span end, as returned by noodles `alignment_end()`.
    pub end1:   SeqPos1,
}
impl Default for ReferenceSpan {
    /// Return an empty reference span with all values set to 0, where 
    /// `end1` == 0 is a reliable test for whether a span is empty.
    fn default() -> Self {
        ReferenceSpan {
            reference_sequence_id: 0,
            start1: 0,
            end1:   0,
        }
    }
}
impl ReferenceSpan {
    /// Determine if a `ReferenceSpan` has not been initialized, i.e., if 
    /// `end1` == 0.
    pub fn is_empty(&self) -> bool {
        self.end1 == 0
    }
}

/// Coverage and identity counts for sequences on scaffold positions.
#[derive(Clone, Copy)]
pub struct Coverage {
    n_seqs: u16,
    n_identical: u16,
}

/* -----------------------------------------------------------------------------
Resolver - scaffolded multiple sequence alignment and consensus
----------------------------------------------------------------------------- */
/// A reusable minimap2-assisted POA consensus resolver.
pub struct Resolver {
    // Fields set per resolver by `with_capacity()`.
    cfg: ResolverConfig,
    base_capacity: usize,

    // Fields set per scaffold by `set_scaffold[_with_aligner]()`.
    scaffold: Vec<BaseByte>,
    scaffold_span: ReferenceSpan,
    scaffold_len: usize,
    coverage: Vec<Coverage>, // map of which seqs are identical at this position
    aligner: Option<Aligner>,

    // Sequence Vecs filled per scaffold by `add_seq()`. `seqs` are cleared and 
    // re-instantiated per scaffold, `seq_maps` are re-used in place.
    seqs: Vec<Vec<BaseByte>>,
    seq_maps: Vec<Vec<Option<SeqPos0>>>, // outer=seq, inner=scaffold SeqPos0, value=seq SeqPos0

    // Fields used during seq_map construction.
    seq_pos0: SeqPos0,
    scaffold_pos0: SeqPos0,
}
impl Resolver {

    /* -------------------------------------------------------------------------
    Resolver - public API
    ------------------------------------------------------------------------- */
    /// Create a new reusable `Resolver` consensus engine with the indicated 
    /// sequence and sequence length capacity. The returned engine can be used 
    /// for all consensus resolutions that uses the same configuration 
    /// parameters. Capacity will grow and persist stably as needed if the 
    /// number of sequences or their lengths exceed the values provided here.
    pub fn with_capacity(
        cfg: ResolverConfig,
        n_seqs:  usize,
        n_bases: usize,
    ) -> (Self, PoaPool) {
        let poa_config: PoaConfig = PoaConfig {
            band_width: 0, // instead, we update band_width based on seq len differences below
            adaptive_band: false,
            alignment_mode: AlignmentMode::Global, // since all POA is anchored
            ..PoaConfig::default()
        };
        let seq_map = Vec::with_capacity(n_bases);
        let n_poa = n_bases / 10;
        (
            Self {
                cfg,
                base_capacity: n_bases,
                scaffold:     Vec::with_capacity(n_bases),
                scaffold_span: ReferenceSpan::default(),
                scaffold_len: 0,
                coverage: Vec::with_capacity(n_bases),
                aligner: None,
                seqs:         Vec::with_capacity(n_seqs), // `seqs` is cleared and re-instantiated for each scaffold
                seq_maps: vec![seq_map; n_seqs], // individual `seq_maps` are reused for all scaffolds
                seq_pos0: 0,
                scaffold_pos0: 0,
            },
            PoaPool::with_capacity(poa_config, n_seqs, n_bases, n_poa),
        )
    }

    /// Reset a `Resolver` to begin a new consensus event by adding a new
    /// scaffold sequence. When initialized using `resolver.set_scaffold()`,
    /// all sequences must be pre-aligned to scaffold and added using 
    /// `resolver.add_aln()`. Use `set_scaffold_with_aligner()` if you will add
    /// uanligned sequences using `resolver.add_seq()`. Scaffolds are coerced to 
    /// uppercase ACGTN upstream of further analysis. A `scaffold_span` is
    /// required to match alignment cs tags to the scaffold coordinate space.
    pub fn set_scaffold(
        &mut self,
        scaffold: &[BaseByte],
        scaffold_span: ReferenceSpan,
    ) {
        self.scaffold.clear();
        self.scaffold.extend(scaffold.iter().map(|base| match base {
            b'A' | b'a' => b'A',
            b'T' | b't' => b'T',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            _ => b'N',
        }));
        self.scaffold_span = scaffold_span;
        self.scaffold_len = scaffold.len();
        self.coverage.clear();
        self.coverage.resize(scaffold.len(), Coverage{n_seqs: 0, n_identical: 0});
        self.aligner = None;
        self.seqs.clear(); // seqs is cleared to reset, seq_maps is not
    }

    /// Reset a `Resolver` to begin a new consensus event by adding a new
    /// scaffold sequence, where at least some sequences will be added using 
    /// `resolver.add_seq()` and thus require a minimap2 aligner.
    pub fn set_scaffold_with_aligner(
        &mut self,
        scaffold: &[BaseByte],
    ) {
        self.set_scaffold(scaffold, ReferenceSpan::default());
        // initialize the scaffold pairwise aligner
        let mut aligner = Aligner::from_seqs(
            vec![("scaffold".to_string(), self.scaffold.clone())], 
            self.cfg.aligner_preset
        );
        // ensure that the CIGAR string is reported with =/X ops instead of M
        aligner.output_config_mut().eqx = true;
        self.aligner = Some(aligner);
    }

    /// Use minimap2-assisted scaffold maps to isolate anchored variant spans 
    /// and resolve them by multi-sequence POA to return the consensus sequence. 
    /// The scaffold gets a vote, and is the deciding vote when all sequences 
    /// are different at a position.
    pub fn get_consensus(
        &mut self,
        poa_pool: &mut PoaPool,
    ) -> Vec<BaseByte> {

        // abort and return scaffold if too few valid alignments exist to override it
        // this is not considered an error condition
        if self.seqs.len() < 2 { return self.scaffold.clone() }

        // determine which scaffold positions must be resolved by POA
        let is_poa: Vec<_> = self.coverage.iter().map(|coverage|{
            // as above, it takes at least two other votes to override scaffold
            coverage.n_seqs >= 2 && 
            // most cases resolve to scaffold without POA if simple majority is identical
            coverage.n_identical <= coverage.n_seqs / 2 
        }).collect();
        
        // establish POA and non-POA chunks, where non-POA chunks must have sufficient anchor length
        let chunks: Vec<_> = is_poa
            .chunk_by(|a, b| a == b)
            .map(|is_poa| {
                let chunk_len = is_poa.len();
                (
                    is_poa[0] || chunk_len < self.cfg.anchor_len * 2, 
                    chunk_len
                )
            })
            .collect();

        // collapse POA chunks when they included too-short non-POA chunks
        // establish the chunk map into scaffold
        let mut offset = 0;
        let chunks: Vec<_> = chunks
            .chunk_by(|a, b| a.0 == b.0)
            .map(|chunks| {
                let chunk_len: usize = chunks.iter().map(|c| c.1).sum();
                let chunk = (chunks[0].0, offset, chunk_len);
                offset += chunk_len;
                chunk
            })
            .collect();

        // ensure sufficient POA capacity for parallel processing
        let n_chunks = chunks.len();
        poa_pool.fill_to(n_chunks / 2 + 1);

        // solve chunks in parallel and flatten to the output consensus
        let max_chunk_i0 = n_chunks - 1;
        // eprintln!("  solve parallel {}", n_chunks);
        chunks.par_iter()
            .zip(&mut poa_pool.poas)
            .enumerate()
            .map(|(chunk0, ((is_poa, offset, chunk_len), poa))|{
                if *is_poa {
                    let scaffold_start0 = offset.saturating_sub(self.cfg.anchor_len);
                    let scaffold_end1 = (offset + chunk_len + self.cfg.anchor_len).min(self.scaffold_len);
                    self.execute_poa(poa, scaffold_start0, scaffold_end1)
                } else {
                    let scaffold_start0 = if chunk0 == 0 {
                        0
                    } else {
                        offset + self.cfg.anchor_len
                    };
                    let scaffold_end1 = if chunk0 == max_chunk_i0 {
                        self.scaffold_len
                    } else {
                        offset + chunk_len - self.cfg.anchor_len
                    };
                    self.scaffold[scaffold_start0..scaffold_end1].to_vec()
                }
            }).flatten().collect()
    }
    /* -------------------------------------------------------------------------
    Resolver - internal function called by get_consensus
    ------------------------------------------------------------------------- */
    fn execute_poa(
        &self,
        poa: &mut Poa,
        scaffold_start0: usize,
        scaffold_end1:   usize,
    ) -> Vec<BaseByte> {

        // collect the starts and ends of each sequence in its own coordinate space
        // find the length difference between the longest and shortest sequence
        let mut min_len = scaffold_end1 - scaffold_start0;
        let mut max_len = min_len;
        let seq_ranges: Vec<_> = self.seqs.iter().enumerate()
            .map(|(seq0, _seq)| { // don't iterate over seq_maps, which is not cleared per scaffold
                let seq_map = &self.seq_maps[seq0];
                let Some(start0) = seq_map[scaffold_start0]   else { return None; };
                let Some(end0)   = seq_map[scaffold_end1 - 1] else { return None; };
                let len = end0 - start0 + 1;
                if min_len > len { min_len = len }
                if max_len < len { max_len = len }
                Some((start0, end0 + 1))
            }).collect();

        // set the required bandwidth
        let bandwidth = max_len - min_len + 1;

        // seed the POA graph
        // eprintln!("scf {}", std::str::from_utf8(&self.scaffold[scaffold_start0..scaffold_end1]).unwrap());
        poa.seed_new_graph( 
            &self.scaffold[scaffold_start0..scaffold_end1], 
            Some(bandwidth)
        );

        // add/align all other sequences to the graph
        for (seq0, seq) in self.seqs.iter().enumerate(){
            // eprintln!("{seq0}");
            if let Some(range) = seq_ranges[seq0] {
                // eprintln!("seq {}", std::str::from_utf8(&seq[range.0..range.1]).unwrap());
                poa.add_seq(&seq[range.0..range.1]);
            }     
        }
        poa.get_heaviest_path()
    }
    /* -------------------------------------------------------------------------
    Resolver - internal functions called by add_seq, add_aln
    ------------------------------------------------------------------------- */
    /// Add an owned copy of a new sequence from a u8 iterator to a Resolver 
    /// as ACGTN.
    fn add_sequence<'a, I: Iterator<Item = u8>>(&mut self, bytes: I) -> usize {
        let seq0 = self.seqs.len();
        self.seqs.push(bytes.map(|base| match base {
            b'A' | b'a' => b'A',
            b'T' | b't' => b'T',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            _ => b'N',
        }).collect());
        seq0
    }

    /// Prepare to add a new sequence alignment to a Resolver scaffold map by 
    /// resetting the indexed entry in `self.seq_maps`.
    fn reset_scaffold_map(&mut self, seq0: usize) {
        if let Some(seq_map) = self.seq_maps.get_mut(seq0){
            seq_map.clear();
            seq_map.resize(self.scaffold_len, None);
        } else {
            let base_capacity = self.scaffold_len.max(self.base_capacity);
            let mut seq_map: Vec<Option<SeqPos0>> = Vec::with_capacity(base_capacity);
            seq_map.resize(self.scaffold_len, None);
            self.seq_maps.push(seq_map);
        }
    }

    /// Include any left-clipped bases of a sequence in the consensus in
    /// end-to-end mode. 
    fn fill_left_clip(&mut self, seq0: usize) {
        // o--O------           overhang beyond scaffold, aligned continuously to reference (`add_aln` only)
        // xxxo------           overhang beyond scaffold, clipped, e.g., an SV junction
        //  xxXxxo--------      overhang beyond scaffold, clipped at internal position
        //    ...Xxo------      end clip on internal alignment, no extension beyond scaffold
        //    ...o-------       continuous internal alignment to scaffold
        //    o------------     continuous complete alignment to scaffold
        // +++=============+++  scaffold without reference flanks 

        // if clip is too far from the left end of scaffold, never use the clip
        //  - internal clips are ignored
        //  - expect other reads with random ends to cover the potential variant span
        // if clip is close enough to the left end, force variant span to scaffold start

        // ==P===X====A====   P=X====A====        P=A====
        //       XxxxxA----     XxxxxA----     XxxxxA----
        // =P===X== ==A====   PX== ==A====        P=A====
        //       XxxxxA----     XxxxxA----     XxxxxA----
        // ==P===X====A====   P=X====A====        P=A====
        //      Xxx xxA----    Xxx xxA----    Xxx xxA----
        if self.cfg.is_end_to_end && self.seq_pos0 > 0 {
            if self.scaffold_pos0 > 0 {
                for coverage in &mut self.coverage[0..self.scaffold_pos0]{
                    coverage.n_seqs += 1;
                }
            }
            let n_scaffold_bases = self.scaffold_pos0.max(1);
            self.seq_maps[seq0][0..n_scaffold_bases].fill(Some(0));
        }
    }

    /// Include any right-clipped bases of a sequence in the consensus in
    /// end-to-end mode. 
    fn fill_right_clip(&mut self, seq0: usize){ 
        //          ------O--o  overhang beyond scaffold, aligned continuously to reference (`add_aln` only)
        //          ------oxxx  overhang beyond scaffold, clipped, e.g., an SV junction
        //     --------oxxXxx   overhang beyond scaffold, clipped at internal position
        //     ------oxX...     end clip on internal alignment, no extension beyond scaffold
        //       ------o...     continuous internal alignment to scaffold
        //    ------------o     continuous complete alignment to scaffold
        // +++=============+++  scaffold with reference flanks 
        if !self.cfg.is_end_to_end { return }
        let seq_len = self.seqs[seq0].len(); 
        if self.seq_pos0 < seq_len {
            if self.scaffold_pos0 < self.scaffold_len {
                for coverage in &mut self.coverage[self.scaffold_pos0..self.scaffold_len]{
                    coverage.n_seqs += 1;
                }
            }
            let scaffold_start0 = self.scaffold_pos0.min(self.scaffold_len - 1);
            let seq_clip_end0 = Some(seq_len - 1);
            self.seq_maps[seq0][scaffold_start0..self.scaffold_len].fill(seq_clip_end0);
        }
    }

}
