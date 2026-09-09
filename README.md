# mini-consensus

A pure-Rust library for resolving consensus sequences from a set of DNA inputs 
where (i) sequence alignments to a scaffold sequence identify identical spans as 
anchors, and (ii) banded partial order alignment (POA) resolves spans between 
the anchors where the sequences differ.

This crate imports 
[rammap](https://github.com/jwanglab/rammap) 
as its pure-Rust
[minimap2](https://github.com/lh3/minimap2)
implementation and carries POA code re-implemented from the 
[poa-consensus](https://github.com/Psy-Fer/poa-consensus)
crate.

Throughout, we say that we 'resolve' a consensus sequence from its inputs to 
avoid confusion with words like 'build', 'find', 'map', or 'align' that have 
other meanings in Rust and bioinformatics.

## Quick start

```toml
# cargo.toml
[dependencies]
mini_consensus = { git = "https://github.com/wilsontelab/mini_consensus", branch = "main" }
```

```rust
use mini_consensus::*;

// your sequences will be longer
let scaffold = b"GAAATAAGAACCGGCAAATCCTACACTAATCCCTCCACACCCAACATTGAAGACTGATGTA"; 
let seqs: Vec<&[u8]> = vec![
    b"GAAATAAGAACCGGCAAATCCTACACTAATCCCTCCACACCCAACATTGAAGACTGATGTA",
    b"GAAATAAGAACCGGCAAATCCTACACTAATCCCCTCCACACCCAACATTGAAGACTGATGTA",
];
let mut resolver = Resolver::with_capacity(
    ResolverConfig::default(),
    seqs.len(), 
    scaffold.len() * 2
);
resolver.set_scaffold_with_aligner(scaffold);
for seq in &seqs { 
    resolver.add_seq(seq);
}
if let Some(consensus) = resolver.get_consensus() {
    // use the consensus, which is Vec<u8>
}
// use the resolver iteratively with new scaffold and seqs
```

## Use cases

We initially developed `mini_consensus` around the first two use cases below and
are extending it to the others.

- three-strand long-read duplex error correction (one reference, two strands)
- local long-read contig assembly around complex transgene insertions
- extended reference-assisted whole-genome assembly
- de novo assembly

## Scaffold-assisted iterative local POA

POA aligns DNA sequences into a directed acyclic graph (DAG) using affine-gap 
dynamic programming and extracts a consensus by following the heaviest 
(most-supported) path. Insertions and deletions create branches that are 
resolved by all sequences together, which reduces consensus errors in regions 
with clustered differences that occur when building a consensus by sequential
pairwise alignment. However, POA slows considerably as sequence length grows. 
`mini_consensus` uses two strategies to keep it small and fast. 

First, code adapted from `poa_consensus` supports **adaptively banded POA** that 
only computes part of the POA matrix for input sequences that are mostly 
co-linear with only small insertions and deletions. 

More importantly, `mini_consensus` uses `minimap2` to 
**align each input sequence to a scaffold sequence to identify anchor spans**
where all sequences are identical. Anchor spans are committed as is, which 
**restricts POA to small local spans** between the anchors to yield a fast and 
accurate consensus. 

```
seq1       -------GA----------------TT-------
seq2       -------GA----------------  -------
seq3       -------GA----------------T -------
scaffold   =======AG================T =======
method     ==PPPPPPPPPPPP======PPPPPPPPPPPP==  '=' committed as is, 'P' subjected to POA
consensus  =======GA================T =======
```

## Relationship between scaffold and other sequences

`mini_consensus` algorithms are generalized and broadly applicable, but some
forethought is needed to obtain meaningful results.

### Scaffold selection

The choice of scaffold can influence the final result because 
**the scaffold gets a consensus vote**. Unlike `poa_consensus`, `mini_consensus` 
does not select a seed or scaffold for you as the proper choice will depend on
your application. 

Although `mini_consensus` does not depend on what your input sequences are, they 
are typically long HiFi sequencing reads. The scaffold might be one of those  
reads chosen to be representative, e.g., `poa_consensus` recommends the median  
length read. Alternatively, you may choose to use an external reference sequence  
as scaffold, noting that by design 
**scaffold base values are reported when all input sequences have a different value**, 
e.g., during three-strand error correction.

### Sequences can be unaligned or pre-aligned to scaffold

If you already aligned your sequences to the scaffold (e.g., reads to a  
reference genome), the resolver can work from your existing `minimap2` BAM 
records with `cs:Z:` tags using `resolver.add_aln()`, otherwise it can align 
reads to the scaffold for you using `resolver.add_seq()`.

### Sequence orientation

**Consensuses are reported in the strand orientation of the scaffold**. You do 
not need to pre-orient unaligned sequences because `rammap` (minimap2) 
inherently orients reads during alignment to the scaffold. If you use 
pre-existing BAM records, be sure they were aligned to a reference in the same 
orientation as your scaffold sequence so that the `cs:Z:` tags are already in 
scaffold orientation.

## Usage mode #1 - end-to-end sequences over scaffold

The `mini_consensus` resolver operates in two distinct modes based on parameter
`is_end_to_end`. In **end-to-end mode**, i.e, when `is_end_to_end` is `true`, 
the resolver expects that all sequences span the entire scaffold. All bases of 
all sequences contribute to the consensus, including any outer clips.

Because clips are used to resolve consensuses in end-to-end mode, it is best 
practice to **trim adapters and other irrelevant bases** from your input 
sequence ends to prevent them from contaminating the consensus result.

```
seq1      cgt-----A--------  'cgt' bases are used even though clipped by minimap2 
seq2      cgt-----A--------
seq3      cgt-----A------xx  'x' indicates non-consensus clipped bases
seq4      cgt-----A--------
scaffold  CGA=====G========
consensus CGT=====A========
```

Additionally, outer clips might cause the consensus to extend beyond the end of 
the scaffold in end-to-end mode depending on your inputs. Because missing bases 
are treated as a possible outcome, the consensus will extend to the point of 
median coverage 

```
seq1      --------------A--------
seq2      --------------A--------
seq3         -----------A--------
seq4            --------A--------
scaffold        ========G========
consensus    ---========A========
```

A typical use case for end-to-end consensus resolution is restriction enzyme 
fragments with known, fixed endpoints.

## Usage mode #2 - sequences with random coverage over scaffold

In **random coverage mode**, i.e, when `is_end_to_end` is `false`, the 
`mini_consensus` resolver makes no assumptions about the coverage span of any
given sequence on scaffold. Sequences may only align to a portion of the 
scaffold and only contribute to the portion of the consensus they overlap. Outer 
**sequence clips are ignored**, including any extensions past the end of the 
scaffold.

```
seq1           xx----A----        'x' indicates unused clipped bases
seq2                xa------      'xa' is clipped due to variant-end proximity
seq3             ----A--------xxx
seq4      xxx--------A--------xxx
scaffold     ========G========
coverage     22224444455554433     coverage including the voting scaffold
consensus    ========A========
```

Random coverage mode is inherently 
**resistant to adapter contamination and chimeric supplemental alignments** on 
the ends of read sequences as these bases are not used during consensus 
resolution.

Some end clips will occur at true variants when the alignment penalities yield
a higher score without a terminal segment, e.g., seq2 above. Sufficient coverage 
is therefore important to ensure that the remaining reads have sufficient flanks 
to productively align and contribute to proper variant resolution.

Consensuses always span the entire scaffold. If no sequences cover any scaffold
positions those positions are reported as scaffold base values to maintain
contiguity, including at the ends of the scaffold. If low or zero coverage 
consensus spans are not desirable, it is the responsibility of the caller to 
split scaffolds into chunks with sufficient coverage.

A typical use case for random coverage consensus resolution is when sequences
are reads and the scaffold is a portion of the reference genome selected to be
free of SV junctions, i.e., to resolve a haplotype consensus over a contig.

## Performance benchmarking

How fast and accurate is it? As one test, we iteratively generated a random 
DNA sequence as a scaffold as well as a set of sequences with randomly 
introduced base substitutions and indels relative to that scaffold. The resolved 
end-to-end consensus is expected to match the scaffold in most instances with
sufficent coverage, where the added changes are a model for sequencing errors.

We also created one sequence with random changes relative to the scaffold and 
used that derived sequence to generate additional sequences with further changes. 
Now the resolved consensus in most instances is expected to match the first
added sequence, not the scaffold. The changes added to the first sequence are 
a model for clonal variants that should be found in the consensus.

We performed these tests for 1K iterations at a range of scaffold sizes, 
sequence counts, and random variant densities and measured the elapsed time and 
frequency of consensuses that matched expectations above. Results are tabulated 
below (times include random sequence generation but this is fast relative to 
consensus resolution).

PENDING (most cases are sub-second resolution of long-read consensuses)

## Other differences between `poa_consensus` and `mini_consensus`

The iterative local POA calls made during anchor-assisted consensus resolution 
prompted minor implementation revisions around the core logic code taken from 
`poa-consensus`. 

First, `mini_consensus::Poa` uses a **single pre-built POA graph engine** with 
pre-allocated buffers that are iteratively reset with new scaffolds to promote
computational efficiency for genome-scale applications.

Second, `mini_consensus::Poa` resolves all graph branches using 
**heaviest bundle logic**, and **gives preference to the scaffold sequence**
when breaking ties.

Finally, `mini_consensus::Poa` internally adds **non-DNA clamp sequences** on 
the flanks of the scaffold and all sequences, which forces the outer ends of 
sequences to align to each other and encourages proper consensus resolution, 
especially when applying POA to variant spans at the end of the scaffold. The 
clamp bases are removed prior to returning the consensus.

```
seq1     EFFJFJJEFJ--------A--------LPPQPQQLPQ  'FFJFJJEFJ' and 'LPPQPQQLPQ' = non-IUPAC clamps
seq2     EFFJFJJEFJ  ------A--------LPPQPQQLPQ
seq3     EFFJFJJEFJ--  ----A--------LPPQPQQLPQ
scaffold EFFJFJJEFJ========G========LPPQPQQLPQ
```

If these changes are useful to you for performing POA without minimap2 support, 
you can access the `mini_consensus::Poa` module directly, starting with 
`let poa = mini_consensus::Poa::with_capacity(...)`.

## Licenses

MIT - James Ferguson 2026 (borrowed POA base code)  
MIT - Thomas Wilson 2026 (POA modifications and the remaining implementation)
