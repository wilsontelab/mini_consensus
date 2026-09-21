//! Parallized version of `add_seq()`.

// imports
use rammap::{Strand, Mapping};
use super::*;

struct ParResult {
    aligned: bool,
    seq: Vec<BaseByte>,
    seq_map: Vec<Option<SeqPos0>>,
    coverage: Vec<Coverage>,
}

impl Resolver {

    /// Parallized version of `add_seq()`.
    pub fn par_set_seqs<V>(
        &mut self,
        seqs: Vec<Vec<BaseByte>>,
        validate: V,
    ) -> Vec<bool> 
    where V: Fn(&Mapping, usize) -> bool + Send + Sync
    {
        let aligner = self.aligner.as_ref().expect(
            "Must call `set_scaffold_with_aligner()` before calling `par_set_seqs()`."
        );
        self.seqs.clear();
        self.seq_maps.clear();
        let mut was_aligned: Vec<bool> = Vec::with_capacity(seqs.len());

        // align and map the sequences in parallel
        let scaffold_len = self.scaffold_len;
        let is_end_to_end = self.cfg.is_end_to_end;
        let par_results: Vec<ParResult> = seqs.into_par_iter().map(|seq|{
            let seq_len = seq.len();
            let mut par_result = ParResult{
                aligned: false,
                seq: seq.into_iter().map(|base| match base {
                    b'A' | b'a' => b'A',
                    b'T' | b't' => b'T',
                    b'C' | b'c' => b'C',
                    b'G' | b'g' => b'G',
                    _ => b'N',
                }).collect(),
                seq_map: vec![None; scaffold_len],
                coverage: vec![Coverage{ n_seqs: 0, n_identical: 0 }; scaffold_len],
            };

            // align seq to scaffold
            let map_result = &aligner.map_seq("seq", &par_result.seq);

            // check mapping validity
            // expect exactly one primary alignment span per input sequence
            if map_result.mappings.len() == 0 { return par_result }
            if map_result.mappings.len() > 1 &&
               map_result.mappings[1].is_primary { return par_result }
            let mapping = &map_result.mappings[0];
            let Some(cigar_ops) = &mapping.cigar_ops else { return par_result };
            if !validate(mapping, seq_len) { return par_result }
            par_result.aligned = true;

            // determine the starting alignment position on scaffold and sequencee
            // as needed, reverse complement seq to scaffold orientation
            for coverage in &mut par_result.coverage[mapping.target_start..mapping.target_end]{
                coverage.n_seqs += 1;
            }
            let mut scaffold_pos0 = mapping.target_start;
            let mut seq_pos0 = if mapping.strand == Strand::Reverse {
                for base in par_result.seq.iter_mut() {
                    *base = match base {
                        b'A' => b'T',
                        b'T' => b'A',
                        b'C' => b'G',
                        b'G' => b'C',
                        _    => b'N',
                    };               
                }
                par_result.seq.reverse();
                seq_len - mapping.query_end
            } else {
                mapping.query_start
            };

            // if end-to-end, fill any left-side clips in the alignment
            if is_end_to_end {
                fill_left_clip(
                    scaffold_pos0,
                    seq_pos0,
                    &mut par_result.seq_map,
                    &mut par_result.coverage,
                );
            } 

            // map the aligned portion of seq
            // eprintln!("{}S, {:?}", mapping.query_start, mapping.cigar);
            for op in cigar_ops {
                process_cigar_op_eqx(
                    &mut scaffold_pos0, 
                    &mut seq_pos0, 
                    &mut par_result.seq_map,
                    &mut par_result.coverage,
                    op,
                );                
            }    

            // if end-to-end, fill any right-side clips in the alignment
            if is_end_to_end {
                fill_right_clip(
                    scaffold_pos0,
                    seq_pos0,
                    &mut par_result.seq_map,
                    &mut par_result.coverage,
                    seq_len,
                    scaffold_len,
                );
            } 
            par_result
        }).collect();

        // collate the alignment results
        for par_result in par_results {
            if par_result.aligned {
                self.seqs.push(par_result.seq);
                self.seq_maps.push(par_result.seq_map);
                assert_eq!(self.coverage.len(), par_result.coverage.len());
                for (rc, prc) in self.coverage.iter_mut()
                    .zip(par_result.coverage) {
                        rc.add_assign(prc);
                }
            }
            was_aligned.push(par_result.aligned);
        }
        was_aligned
    }
}