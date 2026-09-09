// A pure-Rust library for resolving consensus sequences from a set of DNA inputs 
// where (i) sequence alignments to a scaffold sequence identify identical spans as 
// anchors, and (ii) banded partial order alignment (POA) resolves spans between 
// the anchors where the sequences differ.
//! 
//! ```toml
//! # cargo.toml
//! [dependencies]
//! mini-consensus = { git = "https://github.com/wilsontelab/mini_consensus", branch = "main" }
//! ```
//! 
//! ```rust
//! use mini_consensus::*;
//! 
//! // your sequences will be longer
//! let scaffold = b"GAAATAAGAACCGGCAAATCCTACACTAATCCCTCCACACCCAACATTGAAGACTGATGTA"; 
//! let seqs: Vec<&[u8]> = vec![
//!     b"GAAATAAGAACCGGCAAATCCTACACTAATCCCTCCACACCCAACATTGAAGACTGATGTA",
//!     b"GAAATAAGAACCGGCAAATCCTACACTAATCCCCTCCACACCCAACATTGAAGACTGATGTA",
//! ];
//! let mut resolver = Resolver::with_capacity(
//!     ResolverConfig::default(),
//!     seqs.len(), 
//!     scaffold.len() * 2
//! );
//! resolver.set_scaffold_with_aligner(scaffold);
//! for seq in &seqs { 
//!     resolver.add_seq(seq, |_mapping| true); // or validate mapping as needed
//! }
//! if let Some(consensus) = resolver.get_consensus() {
//!     // use the consensus, which is Vec<u8>
//! }
//! // use the resolver iteratively with new scaffold and seqs
//! ```

// modules
pub mod poa;
pub mod minimap2;

// re-exports
pub use poa::*;
pub use minimap2::*;
pub use rammap::Preset;
