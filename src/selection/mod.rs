//! Protein selection language implementation
//!
//! This module provides a PyMOL-style selection language for selecting atoms
//! based on properties like chain, residue name, residue number, element, etc.

use std::collections::HashSet;

pub mod parser;
pub mod evaluator;

/// A set of selected atom indices
#[derive(Debug, Clone, Default)]
pub struct SelectionSet {
    /// Indices of selected atoms in the PDB structure
    pub atom_indices: HashSet<usize>,
}

impl SelectionSet {
    pub fn new(atom_indices: HashSet<usize>) -> Self {
        Self { atom_indices }
    }

    pub fn count(&self) -> usize {
        self.atom_indices.len()
    }

    pub fn contains(&self, index: usize) -> bool {
        self.atom_indices.contains(&index)
    }

    pub fn union(&self, other: &SelectionSet) -> SelectionSet {
        let mut atom_indices = self.atom_indices.clone();
        atom_indices.extend(&other.atom_indices);
        SelectionSet { atom_indices }
    }

    pub fn intersection(&self, other: &SelectionSet) -> SelectionSet {
        let atom_indices = self.atom_indices.iter()
            .filter(|i| other.atom_indices.contains(i))
            .cloned()
            .collect();
        SelectionSet { atom_indices }
    }

    pub fn difference(&self, other: &SelectionSet) -> SelectionSet {
        let atom_indices = self.atom_indices.iter()
            .filter(|i| !other.atom_indices.contains(i))
            .cloned()
            .collect();
        SelectionSet { atom_indices }
    }
}

/// Abstract syntax tree for selection expressions
#[derive(Debug, Clone, PartialEq)]
pub enum SelectionExpression {
    /// Select all atoms in a specific chain
    Chain(String),
    /// Select all atoms in residues with a specific name (e.g., ALA)
    ResidueName(String),
    /// Select all atoms in a residue with a specific serial number
    ResidueNumber(isize),
    /// Select all atoms in a range of residue serial numbers
    ResidueRange(isize, isize),
    /// Select all atoms with a specific name (e.g., CA)
    AtomName(String),
    /// Select all atoms of a specific element (e.g., C)
    Element(String),
    /// Select all atoms in a specific model (1-based index)
    Model(usize),
    /// Select all backbone atoms (N, CA, C, O, OXT)
    Backbone,
    /// Select all non-backbone atoms
    Sidechain,
    /// Select all atoms in alpha helices
    Helix,
    /// Select all atoms in beta sheets
    Sheet,
    /// Select atoms within a certain distance of another selection
    Within(f32, Box<SelectionExpression>),
    /// Logical AND of two selections
    And(Box<SelectionExpression>, Box<SelectionExpression>),
    /// Logical OR of two selections
    Or(Box<SelectionExpression>, Box<SelectionExpression>),
    /// Logical NOT of a selection
    Not(Box<SelectionExpression>),
    /// Select all atoms
    All,
    /// Select no atoms
    None,
}
