//! A pure-Rust library for resolving consensus sequences from a set of DNA inputs
//! where pairwise minimap2 alignments to a scaffold sequence identify identical 
//! spans as anchors and banded partial order alignment (POA) resolves spans between 
//! the anchors where the sequences differ.
//! 
//! ## Quick start
//! 
//! ```toml
//! [dependencies]
//! mini_consensus = "0.1"
//! ```
//! 
//! ```rust
//! use mini_consensus::{Resolver, PoaConfig};
//! 
//! let reads: Vec<&[u8]> = vec![
//!     b"CATCATCAT",
//!     b"CATCATCAT",
//!     b"CATCGTCAT",
//!     b"CATCATCAT",
//! ];
//! ```

pub mod poa;
pub mod minimap2;

pub use poa::*;
pub use minimap2::*;
