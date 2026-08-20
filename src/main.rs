
// use mini_consensus::poa::Poa;
// use rammap::{Aligner, Preset, Strand};
use rammap::Preset;
use rand::prelude::*;
use mini_consensus::*;

const BASES: &[u8] = b"ACGT";

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let n_bases = 250;
    let n_seqs = 10;
    let n_variants = 5;
    let (reference, variants) = generate_dna_with_variants(
        n_bases,
        n_seqs,
        n_variants
    );
    println!("{}", reference);
    for variant in &variants {
        println!("{}", variant);
    }

    // println!("initialize resolver");
    let mut resolver = Resolver::with_capacity(
        Preset::MapHifi, 
        n_seqs, 
        n_bases,
        None
    );
    // println!("setting scaffold");
    resolver.set_scaffold_with_aligner(reference.as_bytes());
    // println!("adding seqs");
    for variant in variants {
    // println!("adding seq");
        resolver.add_seq(variant.as_bytes());
    }
    // println!("getting consensus");
    let consensus = resolver.get_consensus();
    if let Some(consensus) = consensus {
        let consensus = std::str::from_utf8(&consensus).unwrap();
        println!("{}", consensus);
        println!("{}", consensus.len());
        println!("{}", consensus == reference);
        println!("{}", reference);
        println!("{}", consensus);
    } else {
        println!("FAILED!");
    }
// CCGGAACTGTCCTGGGCTTGAGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGGGCTTGCGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATATGCCGGACCATTGCTAAGCTCTGTCGGTAAGTCATGGGCCGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGGGCCTGAGAATGAACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACCTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGCTCGCATAATATGGTGACGTATAGAATAAAGGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGGGCTTGAGAATGTACATCGATCCAGATAACTCTGGTAACCTCTAATAGAGTGCAATTCCGGACCACGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATTAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGTGAACTGTCCTGGGCTTGAGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACTGTCGGTAAGTCATCGGGCCGGGGGTCTTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTCGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGGGCTTGAGAATGTCATCGTCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACTGTCGGTAAGTCATGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCTCCCGCCTTCAAACGTTTT
// CCGGAACTCGTCCTGGGCTTGAGAATGTACATCGATCCGATAACTCTAGTAACCTCTAATTAGAGTGCAATTCCGGACATTGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGAATACCCCGGAGTAGTATAGGTCGCTTAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGGGCTTGAGAATGGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACGTCGGTAAGTCGATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGGCGCATAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGGGCTTGAGAATGTACGTCGATCCGATAAACTCTGGTAACCTCTAATTAGAGTGCAATTCTGGACCATTGCTAAGCACTGTCGGTAAGTCATCGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGGGCTTGAGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCGCGGACCATTGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGATTAGTATAGGTCGCATAATATGGTTACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCTTCAATCGTTTT
// CCGGAACGTGTCCTGGCTTGAGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGTCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAACGAACGTCCCCCCCTTCAAACGTTTT
// 11      38      40      52      15
// 10      41      43      45      5
// 11      25      26      28      4
// 10      17      17      18      2
// 14      35      35      36      2
// 35      11      12      12      2
// CCGGAACTGTCCTGTGGGCTTGAGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// 252
// false
// CCGGAACTGTCCTGGGCTTGAGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT
// CCGGAACTGTCCTGTGGGCTTGAGAATGTACATCGATCCGATAACTCTGGTAACCTCTAATTAGAGTGCAATTCCGGACCATTGCTAAGCACTGTCGGTAAGTCATGGGCCGGGGGTCTTATTTCCCGCAGCTCTGCATACCCCGGAGTAGTATAGGTCGCATAATATGGTGACGTATAGAATAAACGGAGGTCACACAACTAATGATAAGTGGCCAGATTAGAAACGAACGTCCCCCCCTTCAAACGTTTT


    

    // // Reference sequence (20 bp) ACAGAAATTACAGAAGATACCAGATACACGAATAGACCGAGGATAGGAATCCTAGACGCATCG
    // let ref_seq     =  b"ACTAGCAGATTACGCTATTAGGACCCGTATAGGACCGAGACAATGACTTGACTGGACTGACAGAAATTACAGAGATACCAGATACACGAATAGACCGAGGATAGGAATCCTAGACGCATCG";
    // let read_seq    = b"CACTAGCAGATTACGATATTAGGACCCGTATAGGACCGAGACAATGACTTGACTGGACTGACAGAAATTACAGAAGATACCAGATACACGATAGACCGAGGATAGGAATCCTAGACGCATCGCCC";
    // let read_seq_rc = b"GGGCGATGCGTCTAGGATTCCTATCCTCGGTCTATCGTGTATCTGGTATCTTCTGTAATTTCTGTCAGTCCAGTCAAGTCATTGTCTCGGTCCTATACGGGTCCTAATATCGTAATCTGCTAGTG";

// 1-122 Forward
// CigarOp { len: 14, op: 7 }
// CigarOp { len: 1, op: 8 }
// CigarOp { len: 57, op: 7 }
// CigarOp { len: 1, op: 1 }
// CigarOp { len: 16, op: 7 }
// CigarOp { len: 1, op: 2 }
// CigarOp { len: 32, op: 7 }
// 3-124 Reverse
// CigarOp { len: 14, op: 7 }
// CigarOp { len: 1, op: 8 }
// CigarOp { len: 57, op: 7 }
// CigarOp { len: 1, op: 1 }
// CigarOp { len: 16, op: 7 }
// CigarOp { len: 1, op: 2 }
// CigarOp { len: 32, op: 7 }

    //             7 => '=', // Sequence match
    //             8 => 'X', // Sequence mismatch
    //             1 => 'I', // Insertion to target
    //             2 => 'D', // Deletion from target
    
    //             0 => 'M', // Match or Mismatch
    //             3 => 'N', // Ref skip / intron
    //             4 => 'S', // Soft clip
    //             5 => 'H', // Hard clip
    //             6 => 'P', // Padding

    // let mut aligner = Aligner::from_seqs(
    //     vec![("ref_seq".to_string(), ref_seq.to_vec())], 
    //     Preset::Sr
    // );
    // aligner.output_config_mut().eqx = true;

    // let result = &aligner.map_seq(
    //     "read_seq", 
    //     read_seq,
    // ); 
    // let m = &result.mappings[0];
    // println!("{}-{} {:?}", m.query_start, m.query_end, m.strand);
    // let ops = m.cigar_ops.as_ref().unwrap();
    // for op in ops {
    //     println!("{:?}", op);
    // }
    
    // let result = &aligner.map_seq(
    //     "read_seq_rc", 
    //     read_seq_rc,
    // ); 
    // let m = &result.mappings[0];
    // println!("{}-{} {:?}", m.query_start, m.query_end, m.strand);
    // let ops = m.cigar_ops.as_ref().unwrap();
    // for op in ops {
    //     println!("{:?}", op);
    // }


    Ok(())
}


pub fn generate_dna_with_variants(n_bases: usize, n_seqs: usize, n_variants: usize) -> (String, Vec<String>) {
    let mut rng = rand::thread_rng();

    // 1. Generate a random reference sequence
    let ref_bytes: Vec<u8> = (0..n_bases)
        .map(|_| *BASES.choose(&mut rng).unwrap())
        .collect();

    let reference = String::from_utf8(ref_bytes.clone()).unwrap();

    // 2. Generate derivative sequences
    let mut variants = Vec::with_capacity(n_seqs);

    for _ in 0..n_seqs {
        let mut variant_bytes = ref_bytes.clone();

        for _ in 0..n_variants {
            if variant_bytes.is_empty() {
                break;
            }

            // Pick a random operation: 0 = Substitution, 1 = Insertion, 2 = Deletion
            let mutation_type = rng.gen_range(0..3);
            let pos = rng.gen_range(0..variant_bytes.len());

            match mutation_type {
                // Substitution
                0 => {
                    let current_base = variant_bytes[pos];
                    // Pick a base that is different from the current one
                    let new_base = *BASES
                        .iter()
                        .filter(|&&b| b != current_base)
                        .collect::<Vec<&u8>>()
                        .choose(&mut rng)
                        .unwrap();
                    variant_bytes[pos] = *new_base;
                }
                // Insertion
                1 => {
                    let new_base = *BASES.choose(&mut rng).unwrap();
                    variant_bytes.insert(pos, new_base);
                }
                // Deletion
                2 => {
                    variant_bytes.remove(pos);
                }
                _ => unreachable!(),
            }
        }

        variants.push(String::from_utf8(variant_bytes).unwrap());
    }

    (reference, variants)
}
