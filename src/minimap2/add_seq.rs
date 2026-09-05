//! Resolver support for adding raw sequences as consensus inputs, where 
//! sequences are aligned to the scaffold by this crate using minimap2.

// imports
use rammap::{Strand, CigarOp};
use super::*;

// constants
// 0 => 'M', 1 => 'I', 2 => 'D', 3 => 'N', 4 => 'S', 5 => 'H', 7 => '=', 8 => 'X', _ => '?'
const MATCH:     u8 = 7;
const MISMATCH:  u8 = 8;
const INSERTION: u8 = 1;
const DELETION:  u8 = 2;

impl Resolver {

    /// Use `rammap` (`minimap2`) to establish a map of each sequence on 
    /// scaffold. Like scaffolds, sequences are coerced to uppercase ACGTN 
    /// upstream of further analysis. Sequences are reverse-complemented 
    /// internally to match the scaffold as needed.
    pub fn add_seq(
        &mut self,
        seq: &[BaseByte],
    ) -> Option<()> {

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
        if map_result.mappings.len() == 0 { return None; }
        if map_result.mappings.len() > 1 &&
           map_result.mappings[1].is_primary { return None; }
        let mapping = &map_result.mappings[0];
        let Some(cigar_ops) = &mapping.cigar_ops else { return None };

        // determine the starting alignment position on scaffold and sequencee
        // as needed, reverse complement seq to scaffold orientation
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
        Some(())
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
                self.scaffold_pos0 += op_len;
                self.seq_pos0      += op_len;
            },
            MISMATCH => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                self.is_identical[self.scaffold_pos0..self.scaffold_pos0 + op_len].fill(false);
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
                self.is_identical[self.scaffold_pos0 - 1..=self.scaffold_pos0].fill(false); // force both flanking bases to POA
                self.seq_pos0 += op_len;
            },
            DELETION => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                self.is_identical[self.scaffold_pos0..self.scaffold_pos0 + op_len].fill(false);
                self.seq_maps[seq0][self.scaffold_pos0..self.scaffold_pos0 + op_len]
                    .fill(Some(self.seq_pos0 - 1));
                self.scaffold_pos0 += op_len;
            },
            _ => panic!("Unexpected CIGAR operation: {:?}", op),
        }
    }
}
