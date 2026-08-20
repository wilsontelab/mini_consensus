//! Resolver support for adding raw sequences as consensus inputs, where 
//! sequences are aligned to the scaffold by this crate.

// imports
use rammap::{Strand, CigarOp};
use super::*;

impl Resolver {

    /// Use rammap/minimap2 to establish a map of each sequence on scaffold. 
    /// Internally reverse-complement sequences to match scaffold as needed.
    pub fn add_seq(
        &mut self,
        seq: &[BaseByte],
    ) -> Option<()> {

        // create an owned copy of seq
        let seq0 = self.seqs.len();
        self.seqs.push(seq.to_vec());
        let seq = &mut self.seqs[seq0];

        // align seq to scaffold
        let aligner = self.aligner.as_ref().expect(
            "Must call `set_scaffold_with_aligner()` before calling `add_seq()`."
        );
        let map_result = &aligner.map_seq("seq", seq);

        // check mapping validity
        if map_result.mappings.len() == 0 { return None; }
        let n_primary = map_result.mappings.iter()
            .filter(|m| m.is_primary).count();
        if n_primary > 1 { return None; }
        let mapping = &map_result.mappings[0];
        let Some(cigar_ops) = &mapping.cigar_ops else { return None; };

        // ends do get trimmed with errors near termini
        println!(
            "{}-{} {}-{} {:?}", 
            mapping.query_start + 1,
            mapping.query_end,
            mapping.target_start + 1,
            mapping.target_end,
            &mapping.cigar
        );

        // as needed, reverse complement seq to scaffold orientation
        self.seq_pos0 = if mapping.strand == Strand::Reverse {
            seq.iter_mut().for_each(|base| *base = match base {
                b'A' | b'a' => b'T', // all bases are returned as uppercase
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'N' | b'n' => b'N', // S, W, N are self-complementary
                b'S' | b's' => b'S',
                b'W' | b'w' => b'W',
                b'R' | b'r' => b'Y', // complemented IUPAC codes (R↔Y, K↔M, B↔V, D↔H)
                b'Y' | b'y' => b'R',
                b'K' | b'k' => b'M',
                b'M' | b'm' => b'K',
                b'B' | b'b' => b'V',
                b'V' | b'v' => b'B',
                b'D' | b'd' => b'H',
                b'H' | b'h' => b'D',
                _ => *base,
            });
            seq.reverse();
            seq.len() - mapping.query_end
        } else {
            mapping.query_start
        };

        // let base_capacity = scaffold_len.max(self.base_capacity);
        let scaffold_len = self.is_identical.len();
        let base_capacity = scaffold_len.max(self.base_capacity);
        if let Some(seq_map) = self.seq_maps.get_mut(seq0){
            seq_map.clear();
            seq_map.reserve(base_capacity);
        } else {
            self.seq_maps.push(Vec::with_capacity(base_capacity));
        } 

        // fill any left-side gaps in the alignment
        self.scaffold_pos0 = mapping.target_start;
        if self.scaffold_pos0 > 0 {
            (0..self.scaffold_pos0).for_each(|scaffold_pos0| {
                self.seq_maps[seq0].push(0);
                // self.is_identical[scaffold_pos0] = false;
            });
        }

        // map the aligned portion of seq
        for op in cigar_ops{
            self.process_cigar_op_eqx(seq0, op);                
        }    

        // fill any right-side gaps in the alignment
        while self.scaffold_pos0 < self.is_identical.len() {
            self.seq_maps[seq0].push(0);
            // self.is_identical[self.scaffold_pos0] = false;
            self.scaffold_pos0 += 1;
        }
        Some(())
    }

    /// Process a single minimap2 eqx CIGAR operation to build the scaffold maps.
    fn process_cigar_op_eqx(&mut self, seq0: usize, op: &CigarOp){
        // 0 => 'M', 1 => 'I', 2 => 'D', 3 => 'N', 4 => 'S', 5 => 'H', 7 => '=', 8 => 'X', _ => '?'
        const MATCH:     u8 = 7;
        const MISMATCH:  u8 = 8;
        const INSERTION: u8 = 1;
        const DELETION:  u8 = 2;
        match op.op {
            MATCH => { 
                (0..op.len).for_each(|_| {
                    self.seq_maps[seq0].push(self.seq_pos0);
                    self.seq_pos0 += 1;
                    self.scaffold_pos0 += 1;
                });
            },
            MISMATCH => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                (0..op.len).for_each(|_| {
                    self.seq_maps[seq0].push(self.seq_pos0);
                    self.is_identical[self.scaffold_pos0] = false;
                    self.seq_pos0 += 1;
                    self.scaffold_pos0 += 1;
                });
            },
            INSERTION => {
                //    *III 
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                self.is_identical[self.scaffold_pos0 - 1] = false; // force both flanking bases to POA
                self.is_identical[self.scaffold_pos0]     = false;
                self.seq_pos0 += op.len as usize;
            },
            DELETION => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                (0..op.len).for_each(|_| {
                    self.seq_maps[seq0].push(self.seq_pos0 - 1);
                    self.is_identical[self.scaffold_pos0] = false;
                    self.scaffold_pos0 += 1;
                });
            },
            _ => panic!("Unexpected CIGAR operation: {:?}", op),
        }
    }
}
