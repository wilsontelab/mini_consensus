//! The primary `mini_consensus::Resolver` struct is the consensus resolver 
//! that uses rammap/minimap2 to locate identical spans between multiple input 
//! sequences and anchored POA to resolve conflicts in the interval regions.

// modules
mod add_aln;
mod add_seq;

// imports
use rammap::{Aligner, Preset};
use crate::poa::{PoaConfig, AlignmentMode, Poa};

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

/* -----------------------------------------------------------------------------
Resolver - scaffolded multiple sequence alignment and consensus
----------------------------------------------------------------------------- */
/// A reusable minimap2-assisted POA consensus resolver.
pub struct Resolver {
    // Fields set per resolver by `with_capacity()`.
    cfg: ResolverConfig,
    base_capacity: usize,
    poa: Poa,

    // Fields set per scaffold by `set_scaffold[_with_aligner]()`.
    scaffold: Vec<BaseByte>,
    scaffold_span: ReferenceSpan,
    scaffold_len: usize,
    is_identical: Vec<bool>, // map of scaffold whether all seqs are identical at this position
    aligner: Option<Aligner>,

    // Sequence Vecs filled per scaffold by `add_seq()` or `add_aln()`. `seqs`
    // are cleared and re-instantiated per scaffold, `seq_maps` are re-used in 
    // place.
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
    ) -> Self {
        let poa_config: PoaConfig = PoaConfig {
            band_width: 0, // instead, we update band_width based on seq len differences below
            adaptive_band: false,
            alignment_mode: AlignmentMode::Global, // since all POA is anchored
            ..PoaConfig::default()
        };
        let seq_map = Vec::with_capacity(n_bases);
        Self {
            cfg,
            base_capacity: n_bases,
            poa: Poa::with_capacity(poa_config, n_seqs, 500),
            scaffold:     Vec::with_capacity(n_bases),
            scaffold_span: ReferenceSpan::default(),
            scaffold_len: 0,
            is_identical: Vec::with_capacity(n_bases),
            aligner: None,
            seqs:         Vec::with_capacity(n_seqs), // `seqs` is cleared and re-instantiated for each scaffold
            seq_maps: vec![seq_map; n_seqs], // individual `seq_maps` are reused for all scaffolds
            seq_pos0: 0,
            scaffold_pos0: 0,
        }
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
        self.scaffold.extend(scaffold.iter().map(|base|{
            match base {
                b'A' | b'a' => b'A',
                b'T' | b't' => b'T',
                b'C' | b'c' => b'C',
                b'G' | b'g' => b'G',
                _ => b'N',
            }
        }));
        self.scaffold_span = scaffold_span;
        self.scaffold_len = scaffold.len();
        self.is_identical.clear();
        self.is_identical.resize(scaffold.len(), true);
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
    ) -> Option<Vec<BaseByte>> {
        let mut chunk_pos0:  usize = 0; // leftmost pos0 of the next encountered chunk in scaffold coordinates
        let mut left_start0: usize = 0; // leftmost pos0 of the uncommitted match span left of POA span
        let mut left_end1:   usize = 0; // righmost pos1 of the uncommitted match span left of POA span
        let mut poa_pending: bool  = false; // if true, left is set and we have a span pending POA
        let mut consensus: Vec<BaseByte> = Vec::new();

        let chunks: Vec<_> = self.is_identical
            .chunk_by(|a, b| a == b)
            .map(|chunk| (chunk[0], chunk.len()))
            .collect();
        // eprintln!("{:?}", chunks);

        for (is_identical, n_chunk_pos) in chunks {

            // in a span where all (or all but one) seqs matched scaffold
            if is_identical {

                // initialize a first (and possibly only) matching span
                if !poa_pending {
                    left_start0 = chunk_pos0;
                    left_end1 = chunk_pos0 + n_chunk_pos;

                // process a prior variant span by POA if sufficient anchors or in last chunk
                // too-short matching anchors are included in the POA span 
                } else if n_chunk_pos >= self.cfg.anchor_len * 2 || 
                          chunk_pos0 + n_chunk_pos == self.scaffold_len {
                    let scaffold_end1 = (chunk_pos0 + self.cfg.anchor_len).min(self.scaffold_len);
                    self.execute_poa(
                        &mut consensus,
                        left_start0,
                        left_end1,
                        scaffold_end1,
                    );

                    // jump the left flank to the current right flank
                    // not including the flanking bases committed with POA above
                    left_start0 = scaffold_end1; 
                    left_end1 = chunk_pos0 + n_chunk_pos;
                    poa_pending = false;
                } 

            // in a span where at least one sequence differed from the scaffold
            } else {
                poa_pending = true;
            }
            chunk_pos0 += n_chunk_pos;
        }

        // fill out any identical right-side base spans not used as POA anchors using scaffold bases
        if !poa_pending {
            if left_start0 < self.scaffold_len {
                consensus.extend(&self.scaffold[left_start0..self.scaffold_len]);
            }

        // handle situation where the last chunk is variant 
        } else {
            self.execute_poa(
                &mut consensus,
                left_start0,
                left_end1,
                self.scaffold_len,
            );
        }

        // return our result
        Some(consensus)
    }
    /* -------------------------------------------------------------------------
    Resolver - internal function called by get_consensus
    ------------------------------------------------------------------------- */
    /// Perform anchored partial order alignment on a span where at least one 
    /// sequence differed from the others in a local region. 
    fn execute_poa(
        &mut self,
        consensus: &mut Vec<BaseByte>,
        left_start0:   usize,
        left_end1:     usize,
        scaffold_end1: usize,
    ) -> usize {

        // fill out any identical left-side base spans not used as POA anchors using scaffold bases
        let scaffold_start0 = left_end1.saturating_sub(self.cfg.anchor_len);
        if scaffold_start0 > left_start0 {
            consensus.extend(&self.scaffold[left_start0..scaffold_start0]);
        }

        // collect the starts and ends of each sequence in its own coordinate space
        // find the length difference between the longest and shortest sequence
        let mut min_len = scaffold_end1 - scaffold_start0;
        let mut max_len = min_len;
        let seq_ranges: Vec<_> = self.seqs.iter().enumerate()
            .map(|(seq0, _seq)| { // don't iterate over seq_maps, which is not cleared per scaffold
                let seq_map = &self.seq_maps[seq0];
                // eprintln!("{scaffold_start0} {:?} {scaffold_end1} {:?}", seq_map[scaffold_start0], seq_map[scaffold_end1 - 1]);
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
        self.poa.seed_new_graph( 
            &self.scaffold[scaffold_start0..scaffold_end1], 
            Some(bandwidth)
        );

        // add/align all other sequences to the graph
        for (seq0, seq) in self.seqs.iter().enumerate(){
            // eprintln!("{seq0}");
            if let Some(range) = seq_ranges[seq0] {
                // eprintln!("seq {}", std::str::from_utf8(&seq[range.0..range.1]).unwrap());
                self.poa.add_seq(&seq[range.0..range.1]);
            }     
        }

        // resolve and append the local consensus by heaviest bundle
        consensus.extend(self.poa.get_heaviest_path());
        scaffold_end1
    }

    /* -------------------------------------------------------------------------
    Resolver - internal functions called by add_seq, add_aln
    ------------------------------------------------------------------------- */
    /// Add an owned copy of a new sequence from a u8 iterator to a Resolver 
    /// as ACGTN.
    fn add_sequence<'a, I: Iterator<Item = u8>>(&mut self, bytes: I) -> usize {
        let seq0 = self.seqs.len();
        self.seqs.push(bytes.map(|base|{
            match base {
                b'A' | b'a' => b'A',
                b'T' | b't' => b'T',
                b'C' | b'c' => b'C',
                b'G' | b'g' => b'G',
                _ => b'N',
            }
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
        // if self.scaffold_pos0 > 0 && self.seq_pos0 > 0 {
        //     let scaffold_clip_start0 = self.scaffold_pos0 as isize - self.seq_pos0 as isize;

        //     eprintln!("{} {} {}", self.scaffold_pos0, self.seq_pos0, scaffold_clip_start0);

        //     // the < 0 check isn't working because GATTTCAT C AAGAGC was changed to GATTTCAT AAGAGC
        //     //  X
        //     // GATTTCATcAAGAGC
        //     //          ||||||
        //     // -GATTTCATAAGAGC
        //     // GATTTCAT-AAGAGC


        //     let (scaffold_clip_start0, seq_clip_start0) = if scaffold_clip_start0 < 0 {
        //         match end_clip_mode {
        //             EndClipMode::UseSome(n_additional_bases) => (
        //                 0,
        //                 self.seq_pos0.saturating_sub(self.scaffold_pos0 + n_additional_bases)
        //             ),
        //             _ => (0, 0), // only UseAll matches here
        //         }
        //     } else {
        //         (scaffold_clip_start0 as usize, 0)
        //     };

        //     eprintln!("{} {}", scaffold_clip_start0, seq_clip_start0);


            // let scaffold_clip_start0 = self.scaffold_pos0.saturating_sub(self.seq_pos0);
            // let seq_clip_start0 = if scaffold_clip_start0 == 0 {
            //     match end_clip_mode {
            //         EndClipMode::UseSome(n_additional_bases) => {
            //             self.seq_pos0.saturating_sub(self.scaffold_pos0 + n_additional_bases)
            //         },
            //         _ => 0, // only UseAll matches here
            //     }
            // } else {
            //     self.seq_pos0.saturating_sub(self.scaffold_pos0)
            // };

            // // let scaffold_clip_start0 = self.scaffold_pos0.saturating_sub(self.seq_pos0);
            // let (scaffold_clip_start0, seq_clip_start0) = match end_clip_mode {
            //     EndClipMode::UseSome(n_additional_bases) => (
            //         self.scaffold_pos0.saturating_sub(self.seq_pos0 + n_additional_bases),
            //         self.seq_pos0.saturating_sub(self.scaffold_pos0 + n_additional_bases)
            //     ),
            //     _ => (0, 0), // only UseAll matches here
            // };

            // scaffold_clip_start0 poa_start0 seq_clip_start0
        if self.cfg.is_end_to_end && self.seq_pos0 > 0 {
            let n_scaffold_bases = self.scaffold_pos0.max(1);
            self.is_identical[0..n_scaffold_bases].fill(false);
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
            let scaffold_start0 = self.scaffold_pos0.min(self.scaffold_len - 1);
            let seq_clip_end0 = Some(seq_len - 1);
            self.is_identical[scaffold_start0..self.scaffold_len].fill(false);
            self.seq_maps[seq0][scaffold_start0..self.scaffold_len].fill(seq_clip_end0);
        }


        // if self.scaffold_pos0 < self.scaffold_len && self.seq_pos0 < seq_len {
        //     let clip_len = seq_len - self.seq_pos0; // all bases, including those overhanging scaffold
        //     let scaffold_clip_end1 = (self.scaffold_pos0 + clip_len).min(self.scaffold_len);
        //     let clip_len = scaffold_clip_end1 - self.scaffold_pos0; // dropping any scaffold overhang
        //     let seq_clip_end0 = Some(match end_clip_mode {
        //         EndClipMode::UseSome(n_additional_bases) => {
        //             (self.seq_pos0 + clip_len + n_additional_bases).min(seq_len)
        //         },
        //         _ => seq_len, // only UseAll matches here
        //     } - 1);
        //     self.is_identical[self.scaffold_pos0..scaffold_clip_end1].fill(false);
        //     self.seq_maps[seq0][self.scaffold_pos0..scaffold_clip_end1].fill(seq_clip_end0);
        // }
    }

}
