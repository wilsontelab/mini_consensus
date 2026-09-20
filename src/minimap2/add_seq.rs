//! Resolver support for adding raw sequences as consensus inputs, where 
//! sequences are aligned to the scaffold by this crate using minimap2.

// imports
use thiserror::Error;
use rammap::{Strand, CigarOp, Mapping};
use super::*;

// constants
// 0 => 'M', 1 => 'I', 2 => 'D', 3 => 'N', 4 => 'S', 5 => 'H', 7 => '=', 8 => 'X', _ => '?'
const MATCH:     u8 = 7;
const MISMATCH:  u8 = 8;
const INSERTION: u8 = 1;
const DELETION:  u8 = 2;

/// Errors encountered while aligning a requested sequence to the scaffold.
#[derive(Error, Debug)]
pub enum AddSeqError {
    #[error("sequence did not align to the scaffold")]
    NoAlignment,
    #[error("sequence had more than one primary alignment to the scaffold")]
    MultiPrimaryAlignments,
    #[error("sequence mapping failed caller-provided validation")]
    FailedValidation,
    #[error("sequence mapping unexpectedly failed to yield CigarOps")]
    MissingCigarOps,
}

impl Resolver {

    /// Use `rammap` (`minimap2`) to establish a map of each sequence on 
    /// scaffold. Like scaffolds, sequences are coerced to uppercase ACGTN 
    /// upstream of further analysis. Sequences are reverse-complemented 
    /// internally to match the scaffold as needed.
    /// 
    /// The `validate` closure allows callers to apply custom methods to accept
    /// or reject the sequence alignment to scaffold before it is added.
    pub fn add_seq<V>(
        &mut self,
        seq: &[BaseByte],
        validate: V,
    ) -> Result<(), AddSeqError> 
    where V: Fn(&Mapping) -> bool
    {
        // create an owned copy of seq
        let seq0 = self.add_sequence(seq.iter().copied());
        let seq_mut = &mut self.seqs[seq0];

        // align seq to scaffold
        let aligner = self.aligner.as_ref().expect(
            "Must call `set_scaffold_with_aligner()` before calling `add_seq()`."
        );
        let map_result = &aligner.map_seq("seq", seq_mut);

        // check mapping validity
        // expect exactly one primary alignment span per input sequence
        if map_result.mappings.len() == 0 { return self.abort_add_seq(AddSeqError::NoAlignment) }
        if map_result.mappings.len() > 1 &&
           map_result.mappings[1].is_primary { return self.abort_add_seq(AddSeqError::MultiPrimaryAlignments) }
        let mapping = &map_result.mappings[0];
        let Some(cigar_ops) = &mapping.cigar_ops else { return self.abort_add_seq(AddSeqError::MissingCigarOps) };
        if !validate(mapping) { return self.abort_add_seq(AddSeqError::FailedValidation) }

        // determine the starting alignment position on scaffold and sequencee
        // as needed, reverse complement seq to scaffold orientation
        for coverage in &mut self.coverage[mapping.target_start..mapping.target_end]{
            coverage.n_seqs += 1;
        }
        self.scaffold_pos0 = mapping.target_start;
        self.seq_pos0 = if mapping.strand == Strand::Reverse {
            for base in seq_mut.iter_mut() {
                *base = match base {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    b'G' => b'C',
                    _    => b'N',
                };               
            }
            seq_mut.reverse();
            seq_mut.len() - mapping.query_end
        } else {
            mapping.query_start
        };

        // reset the scaffold map
        self.reset_scaffold_map(seq0);

        // if end-to-end, fill any left-side clips in the alignment
        self.fill_left_clip(seq0);

        // map the aligned portion of seq
        // eprintln!("{}S, {:?}", mapping.query_start, mapping.cigar);
        for op in cigar_ops {
            self.process_cigar_op_eqx(seq0, op);                
        }    

        // if end-to-end, fill any right-side clips in the alignment
        self.fill_right_clip(seq0);
        
        // return success
        Ok(())
    }

    /// Check whether a sequence is a productive alignment to the scaffold. 
    /// Expect exactly one primary alignment span per input sequence. This 
    /// function does not add the sequence to the resolver.
    pub fn check_seq<V>(
        &self,
        seq: &[BaseByte],
        validate: V,
    ) -> bool
    where V: Fn(&Mapping) -> bool,
    {
        let aligner = self.aligner.as_ref().expect(
            "Must call `set_scaffold_with_aligner()` before calling `check_seq()`."
        );
        let map_result = &aligner.map_seq("seq", seq);
        if map_result.mappings.len() == 0 { return false }
        if map_result.mappings.len() > 1 &&
            map_result.mappings[1].is_primary { return false }
        let mapping = &map_result.mappings[0];
        mapping.cigar_ops.is_some() && validate(mapping)
    }

    /// Return the appropriate error state for a failed sequence addition after
    /// removing the sequence from the buffer. Note that `seq_map` has not yet
    /// been extended, and the scaffold arrays remain as is while awaiting the
    /// next sequence.
    fn abort_add_seq(&mut self, error: AddSeqError) -> Result<(), AddSeqError> {
        self.seqs.pop();
        Err(error)
    }

    /// Process a single eqx CIGAR operation to build a scaffold map. Only 
    /// =, X, +, and - operations are expected and processed. In particular,
    /// outer clipped bases are not present in `rammap::mapping.cigar_ops`.
    #[inline(always)]
    fn process_cigar_op_eqx(&mut self, seq0: usize, op: &CigarOp){
        let op_len = op.len as usize;
        match op.op {
            MATCH => { 
                let mut seq_pos = self.seq_pos0..self.seq_pos0 + op_len;
                self.seq_maps[seq0][self.scaffold_pos0..self.scaffold_pos0 + op_len]
                    .fill_with(|| seq_pos.next());
                for coverage in &mut self.coverage[self.scaffold_pos0..self.scaffold_pos0 + op_len] { 
                    coverage.n_identical += 1; 
                }
                self.scaffold_pos0 += op_len;
                self.seq_pos0      += op_len;
            },
            MISMATCH => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                let mut seq_pos = self.seq_pos0..self.seq_pos0 + op_len;
                self.seq_maps[seq0][self.scaffold_pos0..self.scaffold_pos0 + op_len]
                    .fill_with(|| seq_pos.next());
                self.scaffold_pos0 += op_len;
                self.seq_pos0      += op_len;
            },
            INSERTION => {
                //    *III 
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                for coverage in &mut self.coverage[self.scaffold_pos0 - 1..=self.scaffold_pos0] { 
                    coverage.n_identical = coverage.n_identical.saturating_sub(1); 
                }
                self.seq_pos0 += op_len;
            },
            DELETION => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                self.seq_maps[seq0][self.scaffold_pos0..self.scaffold_pos0 + op_len]
                    .fill(Some(self.seq_pos0 - 1));
                self.scaffold_pos0 += op_len;
            },
            _ => panic!("Unexpected CIGAR operation: {:?}", op),
        }
    }
}
