//! Geometric detection of hydrogen bonds in protein structures

use glam::Vec3;
use crate::protein::structure::ProteinData;

/// Represents a detected hydrogen bond
#[derive(Debug, Clone, Copy)]
pub struct HydrogenBond {
    pub donor_atom_index: usize,
    pub acceptor_atom_index: usize,
    pub distance_in_angstroms: f32,
}

/// Identifies hydrogen bonds based on geometric criteria (distance and optional angle)
pub fn detect_hydrogen_bonds_in_protein(
    protein_data: &ProteinData,
) -> Vec<HydrogenBond> {
    use rayon::prelude::*;

    if let Some(pdb) = &protein_data.underlying_pdb_data {
        // 1. Identify potential donor and acceptor atoms
        // Acceptors: O, N (with lone pairs)
        // Donors: N, O (with hydrogens attached)
        let atom_info_list: Vec<(usize, Vec3, &str)> = pdb.atoms()
            .enumerate()
            .filter_map(|(atom_index, atom_reference)| {
                let element_symbol = atom_reference.element().map(|e| e.symbol()).unwrap_or("?");
                match element_symbol {
                    "O" | "N" => {
                        let pos = atom_reference.pos();
                        Some((atom_index, Vec3::new(pos.0 as f32, pos.1 as f32, pos.2 as f32), element_symbol))
                    },
                    _ => None,
                }
            })
            .collect();
            
        let atom_info_slice = atom_info_list.as_slice();

        // 2. Perform distance-based search (threshold < 3.5 A) in parallel
        let identified_hydrogen_bonds: Vec<HydrogenBond> = atom_info_slice.par_iter()
            .enumerate()
            .flat_map(|(list_index, &(donor_atom_index, donor_position, _))| {
                // Check against all other potential acceptors using another parallel iterator
                atom_info_slice.par_iter()
                    .skip(list_index + 1) // Avoid double counting and self-comparison
                    .filter_map(move |&(acceptor_atom_index, acceptor_position, _)| {
                        let distance_squared = donor_position.distance_squared(acceptor_position);
                        
                        if distance_squared < 12.25 { // 3.5 A squared
                            let distance_value = distance_squared.sqrt();
                            Some(HydrogenBond {
                                donor_atom_index,
                                acceptor_atom_index,
                                distance_in_angstroms: distance_value,
                            })
                        } else {
                            None
                        }
                    })
            })
            .collect();
        
        identified_hydrogen_bonds
    } else {
        Vec::new()
    }
}
