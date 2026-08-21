//! A pure-Rust library for resolving consensus sequences from a set of DNA inputs 
//! where (i) pairwise minimap2 alignments to a scaffold sequence identify identical 
//! spans as anchors, and (ii) banded partial order alignment (POA) resolves spans 
//! between the anchors where the sequences differ.
//! 
//! ```toml
//! # cargo.toml
//! [dependencies]
//! mini_consensus = "0.1"
//! ```
//! 
//! ```rust
//! use mini_consensus::*;
//! 
//! let ref = b"CATCATCAT"; // your sequences will be longer
//! let seqs: Vec<&[u8]> = vec![
//!     b"CATCATTCAT",
//!     b"CATCATCAT",
//!     b"CATCGTCAT",
//!     b"CATCATCAT",
//! ];
//! 
//! let mut resolver = Resolver::with_capacity(Preset::MapHifi, seqs.len(), 100, None);
//! resolver.set_scaffold_with_aligner(ref);
//! for seq in &seqs { resolver.add_seq(seq); }
//! if let Some(consensus: Vec<u8>) = resolver.get_consensus() {
//!     // use the consensus
//! }
//! ```

// modules
pub mod poa;
pub mod minimap2;

// re-exports
pub use poa::*;
pub use minimap2::*;
pub use rammap::Preset;
