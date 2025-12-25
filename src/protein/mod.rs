//! Protein data management and processing
//!
//! This module provides functionality for fetching, parsing, and managing
//! protein structures

pub mod fetch;
pub mod structure;

pub use structure::{ColorScheme, ProteinData, ProteinStore, Representation};
