//! Resolver support for adding BamRecords with pre-existing cs tags as 
//! consensus inputs.

// imports
use std::str::from_utf8_unchecked;
use noodles::bam::{Record as BamRecord};
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::{
    data::field::Value as TagValue,
    cigar::op::Kind as CigarOp,
};
use super::*;

// constants
const MATCH:     u8 = b':';
const MISMATCH:  u8 = b'*';
const INSERTION: u8 = b'+';
const DELETION:  u8 = b'-';

impl Resolver {

    /// Use a prior `minimap2` `cs` tag to establish a map of each sequence on 
    /// scaffold. Like scaffolds, sequences are coerced to uppercase ACGTN 
    /// upstream of further analysis. The sequence is taken from the Bam record,
    /// so it is already reverse-complemented to scaffold orientation.
    pub fn add_aln(
        &mut self,
        aln: BamRecord,
    ) -> Option<()> {

        // create an owned copy of seq
        let seq0 = self.add_sequence(aln.sequence().iter());

        // check mapping validity
        // must be a productive alignment that overlaps the scaffold reference span
        let Some(Ok(TagValue::String(cs)))  = aln.data().get(b"cs")  else { return None };
        let Some(Ok(first_cigar_op))           = aln.cigar().iter().next()   else { return None };
        let Some(Ok(reference_sequence_id)) = aln.reference_sequence_id() else { return None };
        let Some(Ok(reference_start1))   = aln.alignment_start()       else { return None };
        let Some(Ok(reference_end1))     = aln.alignment_end()         else { return None };
        if reference_sequence_id != self.scaffold_span.reference_sequence_id ||
           reference_start1.get() > self.scaffold_span.end1 ||
           reference_end1.get()   < self.scaffold_span.start1 {
            return None
        }


        // o--O------           overhang beyond scaffold, aligned continuously to reference
        //          ------O--o  overhang beyond scaffold, aligned continuously to reference
        // +++=============+++  scaffold without reference flanks 


        // determine the starting alignment position on scaffold and sequence
        self.scaffold_pos0 = reference_start1.get().saturating_sub(self.scaffold_span.start1);
        self.seq_pos0 = match first_cigar_op.kind() {
            CigarOp::SoftClip => first_cigar_op.len(),
            _ => 0,
        };

        // reset the scaffold map
        self.reset_scaffold_map(seq0);

        // if end-to-end, fill any left-side clips in the alignment
        self.fill_left_clip(seq0);

        // handle the case where sequence alignment begins before scaffold start
        let mut bytes = cs.iter();
        let mut op = bytes.next().unwrap();
        let mut i0 = 1;
        let mut val_start0 = 1;
        if reference_start1.get() < self.scaffold_span.start1 {
            while let Some(byte) = bytes.next() {
                if !byte.is_ascii_alphanumeric() {

                    self.process_cs_op(seq0, cs, op, val_start0, i0);
                    op = byte;
                    val_start0 = i0 + 1;
                }
                i0 += 1;
            }
        }

        // map the aligned portion of seq
        while let Some(byte) = bytes.next() {
            // :9*ag:10
            if !byte.is_ascii_alphanumeric() {
                self.process_cs_op(seq0, cs, op, val_start0, i0);
                op = byte;
                val_start0 = i0 + 1;
            }
            i0 += 1;
        }
        self.process_cs_op(seq0, cs, op, val_start0, i0);

        // if end-to-end, fill any right-side clips in the alignment
        self.fill_right_clip(seq0);

        // return success
        Some(())
    }


    // /// Process a single minimap2 cs:Z: tag operation to build a scaffold map.
    // /// Only :, *, +, and - operations are expected and processed.
    // #[inline(always)] 
    // fn process_cs_op(
    //     &mut self, 
    //     seq0: usize, 
    //     cs:   &[u8], 
    //     op:   &u8, 
    //     val_start0: usize,
    //     val_end1:   usize,
    // ){
    //     match *op {
    //         // :[0-9]+   Identical sequence length
    //         MATCH => {
    //             let op_val = unsafe { from_utf8_unchecked(&cs[val_start0..val_end1]) };
    //             let op_len = op_val.parse::<usize>().unwrap();
    //             self.seq_maps[seq0].extend((self.seq_pos0..self.seq_pos0 + op_len).map(Some));
    //             self.scaffold_pos0 += op_len;
    //             self.seq_pos0 += op_len;
    //         },
    //         _   => {},
    //     }
    // }

    /// Process a single minimap2 cs:Z: tag operation to build a scaffold map.
    /// Only :, *, +, and - operations are expected and processed.
    #[inline(always)] 
    fn process_cs_op(
        &mut self, 
        seq0: usize, 
        cs:   &[u8], 
        op:   &u8, 
        val_start0: usize,
        val_end1:   usize,
    ){
        match *op {

            // :[0-9]+   Identical sequence length
            MATCH => {
                let op_val = unsafe { from_utf8_unchecked(&cs[val_start0..val_end1]) };
                let op_len = op_val.parse::<usize>().unwrap();
                let mut seq_pos = self.seq_pos0..self.seq_pos0 + op_len;
                self.seq_maps[seq0][self.scaffold_pos0..self.scaffold_pos0 + op_len]
                    .fill_with(|| seq_pos.next());
                self.scaffold_pos0 += op_len;
                self.seq_pos0      += op_len;
            },

            // *[acgtn][acgtn]   Substitution: target to query
            MISMATCH => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                self.is_identical[self.scaffold_pos0] = false;
                self.seq_maps[seq0][self.scaffold_pos0] = Some(self.seq_pos0);
                self.scaffold_pos0 += 1;
                self.seq_pos0      += 1;
            },

            // +[acgtn]+   Insertion to the target
            INSERTION => {
                //    *III 
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                let op_len = val_end1 - val_start0;
                self.is_identical[self.scaffold_pos0 - 1..=self.scaffold_pos0].fill(false); // force both flanking bases to POA
                self.seq_pos0 += op_len;
            },

            // -[acgtn]+   Deletion from the target
            DELETION => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                let op_len = val_end1 - val_start0;
                self.is_identical[self.scaffold_pos0..self.scaffold_pos0 + op_len].fill(false);
                self.seq_maps[seq0][self.scaffold_pos0..self.scaffold_pos0 + op_len]
                    .fill(Some(self.seq_pos0 - 1));
                self.scaffold_pos0 += op_len;
            },
            _   => panic!("Unexpected operation in cs tag: {}", op),
        }
    }

}
