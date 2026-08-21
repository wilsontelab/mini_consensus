use rand::prelude::*;
use mini_consensus::*;
// use rammap::{Aligner, Strand};

const BASES: &[u8] = b"ACGT";

fn main() -> Result<(), Box<dyn std::error::Error>> {

    /* -------------------------------------------------------------------------
    direct test of POA over a specific challenge area
    ------------------------------------------------------------------------- */
    let cfg = PoaConfig { 
            band_width: 2, // instead, we update band_width based on seq len differences belows
            adaptive_band: false,
            alignment_mode: AlignmentMode::Global, // since all POA is anchored
        ..PoaConfig::default()
    };
    let mut poa = Poa::with_capacity(cfg, 3, 100);

    let reference = b"TCAGGAACCGCGTGATTAACAGAAATGGCT";
    let seq1      = b"TCAGGATCCGCGTATTAACAGAAATGGCT";
    let seq2      = b"TCAGGATCCGCGTCGATTAACAGAAATGGCT";
    let expected  = b"TCAGGATCCGCGTGATTAACAGAAATGGCT";

    poa.seed_new_graph(reference,None);
    poa.add_seq(seq1);
    poa.add_seq(seq2);
    let consensus = poa.get_heaviest_path();
    let consensus = std::str::from_utf8(&consensus).unwrap();
    let reference = std::str::from_utf8(reference).unwrap();
    let seq1 = std::str::from_utf8(seq1).unwrap();
    let seq2 = std::str::from_utf8(seq2).unwrap();
    let expected = std::str::from_utf8(expected).unwrap();
    eprintln!("{}", reference);
    eprintln!("{}", seq1);
    eprintln!("{}", seq2);
    eprintln!("{}", consensus);
    eprintln!("{}", expected);
    eprintln!("{}", consensus == expected);

    /* -------------------------------------------------------------------------
    run two test groups 
    - first has all variants relative to reference, so reference is the expected consensus
    - second has clonal variants introduced into variant0, so it is the expected consensus

    TODO: the second group often fails to give the expected consensus
    the issue seems to always be near one or the other terminus, where ref gets used
    most likely this is due to minimap2 trimming rather than POA
    may want to explore trim clamping? currently has a reference bias at termini,
    and thus a relative inability to call consensus variants near the ends
    ------------------------------------------------------------------------- */
    let n_iter = 50;
    // let n_bases = 4096 * 1;
    let n_bases: usize = 250;
    let n_seqs = 2;
    let n_variants = n_bases / 100;
    let mut resolver = Resolver::with_capacity(
        Preset::MapHifi, 
        n_seqs, 
        n_bases,
        None
    );
    for use_variant0_as_expected in vec![false, true] {
        for iter in 0..n_iter {
            eprintln!("{iter}");
            let (
                reference, 
                variants,
                expected
            ) = generate_dna_with_variants(
                n_bases,
                n_seqs,
                n_variants,
                use_variant0_as_expected
            );
            resolver.set_scaffold_with_aligner(reference.as_bytes());
            for variant in &variants {
                resolver.add_seq(variant.as_bytes());
            }
            let consensus = resolver.get_consensus();
            if let Some(consensus) = consensus {
                let consensus = std::str::from_utf8(&consensus).unwrap();
                if consensus != expected {
                    println!("iter {}\t{}\t{}", iter, consensus.len(), consensus == reference);
                    println!("{}", reference);
                    for variant in &variants {
                        println!("{}", variant);
                    }
                    println!("{}", consensus);
                    println!("{}", expected);
                    println!("");
                }
            } else {
                println!("CONSENSUS FAILED!");
            }
        }        
    }

    /* -------------------------------------------------------------------------
    confirm the behavior of rammap CIGAR, etc.
    ------------------------------------------------------------------------- */
    // let ref_seq     =  b"ACTAGCAGATTACGCTATTAGGACCCGTATAGGACCGAGACAATGACTTGACTGGACTGACAGAAATTACAGAGATACCAGATACACGAATAGACCGAGGATAGGAATCCTAGACGCATCG";
    // let read_seq    = b"CACTAGCAGATTACGATATTAGGACCCGTATAGGACCGAGACAATGACTTGACTGGACTGACAGAAATTACAGAAGATACCAGATACACGATAGACCGAGGATAGGAATCCTAGACGCATCGCCC";
    // let read_seq_rc = b"GGGCGATGCGTCTAGGATTCCTATCCTCGGTCTATCGTGTATCTGGTATCTTCTGTAATTTCTGTCAGTCCAGTCAAGTCATTGTCTCGGTCCTATACGGGTCCTAATATCGTAATCTGCTAGTG";
    // let mut aligner = Aligner::from_seqs(
    //     vec![("ref_seq".to_string(), ref_seq.to_vec())], 
    //     Preset::Sr
    // );
    // aligner.output_config_mut().eqx = true;

    // let result = &aligner.map_seq("read_seq", read_seq); 
    // let m = &result.mappings[0];
    // println!("{}-{} {:?}", m.query_start, m.query_end, m.strand);
    // let ops = m.cigar_ops.as_ref().unwrap();
    // for op in ops { println!("{:?}", op); }
    
    // let result = &aligner.map_seq("read_seq_rc", read_seq_rc); 
    // let m = &result.mappings[0];
    // println!("{}-{} {:?}", m.query_start, m.query_end, m.strand);
    // let ops = m.cigar_ops.as_ref().unwrap();
    // for op in ops { println!("{:?}", op); }

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

    Ok(())
}

/// Generate random sequences with variants.
fn generate_dna_with_variants(
    n_bases: usize, 
    n_seqs: usize, 
    n_variants: usize,
    use_variant0_as_expected: bool,
) -> (String, Vec<String>, String) {
    let mut rng = rand::thread_rng();
    let ref_bytes: Vec<u8> = (0..n_bases)
        .map(|_| *BASES.choose(&mut rng).unwrap())
        .collect();
    let reference = String::from_utf8(ref_bytes.clone()).unwrap();
    let mut variants = Vec::with_capacity(n_seqs);
    let expected = if use_variant0_as_expected {
        let variant0= modify_variant(ref_bytes.clone(), 2, &mut rng);
        let ref_bytes = variant0.as_bytes().to_vec();
        variants.push(variant0);
        for _ in 1..n_seqs {
            let variant= modify_variant(ref_bytes.clone(), n_variants, &mut rng);
            variants.push(variant);
        }
        variants[0].clone()
    } else {
        for _ in 0..n_seqs {
            let variant= modify_variant(ref_bytes.clone(), n_variants, &mut rng);
            variants.push(variant);
        }
        reference.clone()
    };
    (reference, variants, expected)
}

/// Introduce random variants into a sequence.
fn modify_variant(
    mut variant_bytes: Vec<u8>, 
    n_variants: usize, 
    rng: &mut ThreadRng
) -> String {
    for _ in 0..n_variants {

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
                    .choose(rng)
                    .unwrap();
                variant_bytes[pos] = *new_base;
            }
            // Insertion
            1 => {
                let new_base = *BASES.choose(rng).unwrap();
                variant_bytes.insert(pos, new_base);
            }
            // Deletion
            2 => {
                variant_bytes.remove(pos);
            }
            _ => unreachable!(),
        }
    }
    String::from_utf8(variant_bytes).unwrap()
}
