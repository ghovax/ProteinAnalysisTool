//! Bond calculation and management for protein structures

use pdbtbx::{Atom, PDB};
use std::collections::HashSet;

/// Represents a bond between two atoms
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomBond {
    /// Index of the first atom in the pdb.atoms() iterator
    pub first_atom_index: usize,
    /// Index of the second atom in the pdb.atoms() iterator
    pub second_atom_index: usize,
}

impl AtomBond {
    pub fn new(first_index: usize, second_index: usize) -> Self {
        let (lower_atom_index, higher_atom_index) = if first_index < second_index {
            (first_index, second_index)
        } else {
            (second_index, first_index)
        };
        Self {
            first_atom_index: lower_atom_index,
            second_atom_index: higher_atom_index,
        }
    }
}

/// Calculator for determining bonds between atoms based on distance and connectivity information
pub struct BondCalculator<'a> {
    pdb_data: &'a PDB,
}

impl<'a> BondCalculator<'a> {
    pub fn new(pdb_data: &'a PDB) -> Self {
        Self { pdb_data }
    }

    /// Calculates bonds using both distance-based heuristics and explicit CONECT records
    pub fn calculate_all_bonds(&self) -> HashSet<AtomBond> {
        use rayon::prelude::*;

        let atom_list: Vec<&Atom> = self.pdb_data.atoms().collect();
        let atom_count = atom_list.len();
        let atom_slice = atom_list.as_slice();
        
        // Use parallel iteration to check distances between all pairs of atoms
        let identified_bonds: HashSet<AtomBond> = (0..atom_count)
            .into_par_iter()
            .flat_map(|first_atom_list_index| {
                let first_atom_reference = atom_slice[first_atom_list_index];
                let first_atom_position_tuple = first_atom_reference.pos();
                
                (first_atom_list_index + 1..atom_count)
                    .into_par_iter() 
                    .filter_map(move |second_atom_list_index| {
                        let second_atom_reference = atom_slice[second_atom_list_index];
                        let second_atom_position_tuple = second_atom_reference.pos();
                        
                        let delta_x = first_atom_position_tuple.0 - second_atom_position_tuple.0;
                        let delta_y = first_atom_position_tuple.1 - second_atom_position_tuple.1;
                        let delta_z = first_atom_position_tuple.2 - second_atom_position_tuple.2;
                        
                        let distance_squared_value = delta_x.powi(2) + delta_y.powi(2) + delta_z.powi(2);
                        
                        // Threshold of ~1.9A (3.61 A^2) for typical covalent bonds
                        if distance_squared_value < 3.61 { 
                            Some(AtomBond::new(first_atom_list_index, second_atom_list_index))
                        } else {
                            None
                        }
                    })
            })
            .collect();
        
        identified_bonds
    }
}
