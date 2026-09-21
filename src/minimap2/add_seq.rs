//! Resolver support for adding raw sequences as consensus inputs, where 
//! sequences are aligned to the scaffold by this crate using minimap2.

// imports
use thiserror::Error;
use rammap::{Strand, Mapping};
use super::*;

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
        let seq_len = seq_mut.len();

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
        let mut scaffold_pos0 = mapping.target_start;
        let mut seq_pos0 = if mapping.strand == Strand::Reverse {
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
            seq_len - mapping.query_end
        } else {
            mapping.query_start
        };

        // reset the scaffold map
        self.reset_scaffold_map(seq0);

        // if end-to-end, fill any left-side clips in the alignment
        if self.cfg.is_end_to_end {
            fill_left_clip(
                scaffold_pos0,
                seq_pos0,
                &mut self.seq_maps[seq0],
                &mut self.coverage,
            );
        } 

        // map the aligned portion of seq
        // eprintln!("{}S, {:?}", mapping.query_start, mapping.cigar);
        for op in cigar_ops {
            process_cigar_op_eqx(
                &mut scaffold_pos0, 
                &mut seq_pos0, 
                &mut self.seq_maps[seq0],
                &mut self.coverage,
                op,
            );                
        }    

        // if end-to-end, fill any right-side clips in the alignment
        if self.cfg.is_end_to_end {
            fill_right_clip(
                scaffold_pos0,
                seq_pos0,
                &mut self.seq_maps[seq0],
                &mut self.coverage,
                seq_len,
                self.scaffold_len,
            );
        } 
        
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
}
