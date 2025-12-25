//! Protein data management and processing
//!
//! This module provides functionality for fetching, parsing, and managing
//! protein structures

pub mod fetch_rcsb;
pub mod structure;
pub mod bonds;

pub use structure::{ColorScheme, ProteinData, ProteinStore, Representation};
