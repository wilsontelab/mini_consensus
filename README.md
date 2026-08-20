# mini-consensus

A pure-Rust library for resolving consensus sequences from a set of DNA inputs
where pairwise minimap2 alignments to a scaffold sequence identify identical 
spans as anchors and banded partial order alignment (POA) resolves spans between 
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
[dependencies]
mini_consensus = "0.1"
```

```rust
use mini_consensus::{Resolver, PoaConfig};

let reads: Vec<&[u8]> = vec![
    b"CATCATCAT",
    b"CATCATCAT",
    b"CATCGTCAT",
    b"CATCATCAT",
];
```

## Use cases

We developed `mini_consensus` around the first two use cases below, but imagine 
it could be applied to the others as well. 
- three-strand long-read duplex error correction (one reference, two strands)
- local long-read contig assembly around complex transgene insertions
- extended reference-assisted whole-genome assembly
- de novo assembly

## Analysis strategy

POA aligns DNA sequences into a directed acyclic graph (DAG) using affine-gap 
dynamic programming and extracts a consensus by following the heaviest 
(most-supported) path. Insertions and deletions create branches that are 
resolved by all sequences together, which reduces consensus errors in regions 
with clustered differences that occur when building a consensus by pairwise 
alignment.

POA slows considerably as sequence length grows. `mini_consensus` uses two 
strategies to keep it small and fast. First, code from `poa_consensus` supports 
adaptively banded POA which only computes part of the POA matrix for input 
sequences that are mostly co-linear with only small insertions and deletions. 
Second, `minimap2` is used to align each input sequence to a scaffold sequence 
to identify spans where all sequences are identical. These anchor spans are 
committed as is, which restricts POA to much smaller spans between the anchors
to yield a fast and accurate consensus. 

If you already aligned your sequences to a scaffold (e.g., reads to a reference 
genome), the resolver can work from your existing minimap2 cs tags using 
`resolver.add_aln()`, otherwise it can align reads to the scaffold for you using
`resolver.add_seq()`.

## Scaffold selection and requirements

Initial sequence scanning uses a single scaffold sequence to which all other
sequences are individually aligned. The choice of scaffold can influence the 
consensus result at the outer boundaries and in high-error regions. Unlike 
`poa_consensus`, `mini_consensus` does not select a POA seed or scaffold for you 
due to the use cases is was built around. 

Although `mini_consensus` does not depend on what your input sequences are, they 
are typically long sequencing reads. The scaffold might be one of those reads 
chosen to be representiative, e.g., `poa_consensus` recommends the median length 
read. Alternatively, you may choose to use an external reference sequence as 
scaffold, noting that by design the scaffold base values will be reported when 
all input sequences have a different value, e.g., during three-strand error 
correction.

Importantly, the scaffold sequence must be "end-to-end" over the entire expected
output consensus. Any sequence portions that overhang the scaffold will be 
trimmed. Also, if the very ends of the consensus are non-identical, the scaffold
base values will be reported out to the first and last scaffold bases. Generally,
this means that your calling code should execute consensus assembly across
regions that are anchored at their boundaries.

In contrast, input sequences do not all need to be end-to-end over the expected
consensus. Sequences only contribute to the portion of the consensus they overlap.

## Sequence orientation

Unlike `poa_consensus`, `mini_consensus` does not auto-orient reads because:
- `rammap` (minimap2) inherently orients reads during alignment to the scaffold
- many other POA use cases can exploit previously oriented sequences

Importantly, the minimap2-assisted consensus resolver will reverse-complement
input sequences as needed to make them match the forward strand of the scaffold,
with sequences passed in as mutable references.

## Other differences between `poa_consensus` and `mini_consensus`

The iterative POA calls made during anchor-assisted consensus resolution 
prompted minor implementation revisions around the core logic code taken from 
`poa-consensus`. First, `mini_consensus::Poa` uses a single pre-built POA graph 
engine with pre-allocated buffers that are iteratively reset with new sequences. 
Also, `mini_consensus` resolves all graph branches using heaviest bundle logic.

If these changes are useful to you for performing POA without minimap2 support, 
you can access the `mini_consensus::Poa` module directly. 

## Licenses

MIT - James Ferguson 2026 (borrowed POA base code)  
MIT - Thomas Wilson 2026 (POA modifications and the remaining implementation)
