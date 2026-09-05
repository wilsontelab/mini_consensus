use rand::prelude::*;
use mini_consensus::*;

const BASES: &[u8] = b"ACGT";
const N_ITER: usize = 100;
const CLONAL_VARIANT_RATE: f64 = 0.01;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    debug_poa();
    // debug_resolver();
    // run_benchmarking();
    Ok(())
}

fn debug_poa(){
    let cfg = PoaConfig { 
            band_width: 5, 
            adaptive_band: false,
            alignment_mode: AlignmentMode::Global, // since all POA is anchored
        ..PoaConfig::default()
    };
    let mut poa = Poa::with_capacity(cfg, 3, 100);
    poa.seed_new_graph(b"CAAGAATCCAATTTCCTAGCTTA", None);
    poa.add_seq(b"CAAGAATCCAATTCCTAGCTTA");
    poa.add_seq(b"CAAGAATCCAATTTCTAGCTTA");
    // poa.get_heaviest_path();
    poa.xxxx();
    // this example doesn't align the way I would do it, see debug.txt
    // increasing the mismatch score will push it to the two deletion graph
}

fn debug_resolver(){
    let ref0 = b"GAAATAAGAACCGGCAAATCCTACACTAATCCCCTCCACACCCAACATTGAAGACTGATGTACCCCGTCTGCGGAGGGGAGACGCGGAGGCATATGTGGTAAAGAAACAGTAATGAATGGTCCACCTCTAAAAACAGACTGGCTCTGCACGGAGAGTCCACGTGGACTACATATACGAATAGATTACCCGACGGGGGTGCGCATCGAGGTTACCCAAATTTCAAGAACCACAGCTGGACTCAGCGATACCCTAGGAGGCAACTCCTTTTGACGGCCAGTAGTACTATTTCGTGTTCAAACCCCGAATAACGCCAGCAAATGTTCGTGACTTTGCAGCCCATCTACCGTAGCGGGAAAGTGGATACATTAGAACAACAACGCCGGACTGGGCCACGGGCTGACATTGGTAAACCCGTTTCCATATCGCATCAGCAAACGATTGGACCCCACCATATCTGCCCCCTCGATGCGTTGTCATCTCGGCAGGGGTTTACCCGCGG";

    let seq0 = b"CGCGATTTCATAAGAGCTCGGACACGGAATCGGGTACTTATCACTGTATACTGCGTAGTCCGAGAGTAATACACAAGCTATAATCGACAACCCGAGAGAGGCCAGTTAATTAGCGATGAAATAAGAACCGGCAAATCCTACACTAATCCCCTCCACACCCAACATTGAAGACTGATGTACCCCTGTCTGCGGAGGGGAGAACGCGGAGGCATATGTGGGTAAAGAAACAGTAATGAATGGTCCACCTCTAAAAACAGACTGGCTCTGCACGGAGAGTCCACGTGGACTACATATACGAATAGATTACCCGACGGGGGTGCGCATCGAGGTTACCCAAATTTCAAGAACCACAGCTGGACTCAGCGATACCCTAGGAGGCAACTCCTTTTGACGGCCAGTAGTACTATTTCGTGTTCAAACCCCGAATAACGCCAGCAAATGTTCGTGACTTTGCAGCCCATCTACCGTAGCGGGAAAGTGGATACATTAGAACAACAACGCCGGACTGGCCACGGGCTGACATTGGTAAACCCGTTTCCATATCGCATCAGCAAACGATTGGACCCCACCATATCTGCCCCCTCGATGCGTTGTCATCTCGGCAGGGGTTTACCCGCGCGAGATAAACTGAGTTCTATTTGGCCGATCAAGCAGCGCCTTGTGCAAGCAACCTGATAATAAGTACCGAGATACGGCTTGTACAGTACGTTGAAGTTTTACGACCAACTCCAACGGATATACAACCCATGTTGATTT";
    let seq1 = b"CGCGATTTCATAAGAGCTCGGACACGGAATCGGGTACTTATCACTGTATACTGCGTAGTCCGAGAGTAATACACAAGCTATAATCGACAACCCGAGAGAGGCCAGTTAATTAGCGATCGAGATGAAATAAGAACCGGCAAATCCTACACTAATCCCCTCCACACCCAACATTGAAGACTGATGTACCCCTGTCTGCGGAGGGGAGAACGCGGAGGCATATGTGGGTAAAGAAACAGTAATGAATGGTCCACCTCTAAAAACAGACTGGCTCTGCACGGAGAGTCCACGTGGACTACATATACGAATAGATTACCCGACGGGGGTGCGCATCGAGGTTACCCAAATTTCAAGAACCACAGCTGGACTCAGCGATACCCTAGGAGGCAACTCCTTTTGACGGCCAGTAGTACTATTTCGTGTTCAAACCCCGAATAACGCCAGCAAATGTTCGTGACTTTGCAGCCCATCTACCGTAGCGGGAAAGTGGATACATTAGAACAACAACGCCGGACTGGCCACGGGCTGACATTGGTAAACCCGTTTCCATATCGCATCAGCAAACGATTGGACCCCACCATATCTGCCCCCTCGATGCGTTGTCATCTCGGCAGGGGTTTACCCGCGAAACTGAGTTCTATTTGGCCGATCAAGCAGCGCCTTGTGCAAGCAACCTGATAATAAGTACCGAGATAACGGCTTGTACAGTACGTTGAAGTTTTACGACCAACTCCAGACGGATATACAACCCATGTTGATTT";
    let seq2 = b"GCGATTTCATAAGAGCTCGGACACGGAATCGGGTACTTATCACTGTATACTGCGTAGTCCGAGAGTAATACACAAGCTATAATCGACAACCCGAGAGAGGCCAGTTAATTAGCGATCGAGATAAACTGAGTTCTATTTGGCCGATCAAGCAGCGCCTTGTGCAAGCAACCTGATAATAAGTACCGAGATACTGGCTTGTACAGTACGTTGAAGTTTTCGACCAACTCCAACGGATATACAACCCATGTTGATTT";
    let seq3 = b"GCGATTTCATAAGAGCTCGGACACGGAATCGGGTACTTATCACTGTATACTGCGTAGTCCGAGAGTAATACACAAGCTATAATCGACAACCCGAGAGAGGCCAGTTAAATTAGCGATCCGAGATAAACTGAGTTCTATTTGGCCGATCAAGCAGCGCCTTGTGCAAGCAACCTGATAATAAGTACCGAGATACGGCTTGTACAGTACGTTGAAGTTTTACGACCAACTCCAACGGATATACAACCCATGTTGATTT";
    let seq4 = b"CGCGATTTCATAAGAGCTCGGACACGGAATCGGGTACTTATCACTGTATACTGCGTAGTCCGAAGTCAATACACAAGCTATAATCGACAACCCGAGAGAGGCCAGTTAATTAGCGATCGAGATAAACTGAGTTCTATTTGGCCGATCAAGCAGCGCCTTGTGCAAGCAACCTGATAATAAGTACCGAGATACGGCTTGTACAGTACGTTGAAGTTTTACGACCAACTCCAACGGATATACAACCCATGTTGATTT";

    let mut resolver = Resolver::with_capacity(
        ResolverConfig {
            is_end_to_end: true,
            ..ResolverConfig::default()
        },
        5, 
        500
    );
    resolver.set_scaffold_with_aligner(ref0);
    resolver.add_seq(seq0).expect("failed to add seq0");
    resolver.add_seq(seq1).expect("failed to add seq1");
    // resolver.add_seq(seq2).expect("failed to add seq2");
    // resolver.add_seq(seq3).expect("failed to add seq3");
    // resolver.add_seq(seq4).expect("failed to add seq4");
    let consensus = resolver.get_consensus().unwrap();
    eprintln!("{}", consensus == ref0);
    eprintln!("{}", consensus == seq0);
    let consensus = std::str::from_utf8(&consensus).unwrap();
    eprintln!("consensus {}", consensus);
}

fn run_benchmarking() {
    println!(
        "{},{},{},{},{},{}",
        "n_bases",
        "n_seqs",
        "error_rate",
        "has_clonal_variants",
        "ms_per_consensus",
        "frac_expected",
    );
    let mut rng = rand::thread_rng();
    for n_bases in vec![500, 5000, 25000] {
        let n_clonal_variants = (n_bases as f64 * CLONAL_VARIANT_RATE) as usize;
        for n_seqs in vec![2, 3, 15] { // 30
            for error_rate in vec![0.001, 0.01] { //, 0.05
                let n_errors = (n_bases as f64 * error_rate) as usize;
                run_iterations(
                    n_bases,
                    n_clonal_variants,
                    n_seqs,
                    error_rate,
                    n_errors,
                    false,
                    &mut rng,
                );
                run_iterations(
                    n_bases,
                    n_clonal_variants,
                    n_seqs,
                    error_rate,
                    n_errors,
                    true,
                    &mut rng,
                );
            }
        }
    }
}

fn run_iterations(
    n_bases: usize,
    n_clonal_variants: usize,
    n_seqs: usize,
    error_rate: f64,
    n_errors: usize,
    has_clonal_variants: bool,
    rng: &mut ThreadRng,
){
    let mut resolver = Resolver::with_capacity(
        ResolverConfig {
            is_end_to_end: true,
            ..ResolverConfig::default()
        },
        n_seqs, 
        n_bases
    );
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(n_seqs);
    let mut n_expected = 0_usize;
    let timer = std::time::Instant::now();
    for _ in 0..N_ITER {
        let scaffold: Vec<u8> = (0..n_bases)
            .map(|_| *BASES.choose(rng).unwrap())
            .collect();
        seqs.clear();
        seqs.push(modify_seq(
            &scaffold,
            n_clonal_variants,
            rng,
        ));
        let expected = if has_clonal_variants { seqs[0].clone() }  else { scaffold.clone() };
        for _ in 1..n_seqs {
            seqs.push(modify_seq(
                &expected,
                n_errors,
                rng,
            ));
        }
        resolver.set_scaffold_with_aligner(&scaffold);
        for seq in &seqs { resolver.add_seq(seq); } //.expect("SEQ ERROR")
        let consensus = resolver.get_consensus().unwrap_or_else(|| Vec::new());
        if consensus == expected {
            n_expected += 1;
        } else {
            println!("==========");
            println!("FAIL, has_clonal_variants: {}", has_clonal_variants);
            println!("{}", std::str::from_utf8(&scaffold).unwrap());
            println!("----------");
            for seq in &seqs {
                println!("{}", std::str::from_utf8(&seq).unwrap());
            }
            println!("----------");
            println!("{}", std::str::from_utf8(&consensus).unwrap());
            println!("");
            panic!("");
        }
    }
    let elapsed_time = timer.elapsed().as_micros();
    println!(
        "{},{},{},{},{},{}",
        n_bases,
        n_seqs,
        error_rate,
        has_clonal_variants,
        (elapsed_time as f64 / N_ITER as f64) / 1000.0,
        n_expected as f64 / N_ITER as f64,
    );
}

fn modify_seq(
    source_bytes: &[u8], 
    n_variants: usize, 
    rng: &mut ThreadRng,
) -> Vec<u8> {
    let mut seq_bytes = source_bytes.to_vec();
    for _ in 0..n_variants {

        // Pick a random operation: 0 = Substitution, 1 = Insertion, 2 = Deletion
        let mutation_type = rng.gen_range(0..3);

        // !!! clamping the scaffold termini !!!
        // let pos = rng.gen_range(0..seq_bytes.len());
        let pos = rng.gen_range(10..seq_bytes.len() - 10);

        match mutation_type {
            // Substitution
            0 => {
                let current_base = seq_bytes[pos];
                // Pick a base that is different from the current one
                let new_base = *BASES
                    .iter()
                    .filter(|&&b| b != current_base)
                    .collect::<Vec<&u8>>()
                    .choose(rng)
                    .unwrap();
                seq_bytes[pos] = *new_base;
            }
            // Insertion
            1 => {
                let new_base = *BASES.choose(rng).unwrap();
                seq_bytes.insert(pos, new_base);
            }
            // Deletion
            2 => {
                seq_bytes.remove(pos);
            }
            _ => unreachable!(),
        }
    }
    seq_bytes
}
