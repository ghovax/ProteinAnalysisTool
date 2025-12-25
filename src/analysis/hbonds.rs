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
    let mut identified_hydrogen_bonds = Vec::new();
    
    // 1. Identify potential donor and acceptor atoms
    // Acceptors: O, N (with lone pairs)
    // Donors: N, O (with hydrogens attached)
    
    let mut potential_acceptor_indices = Vec::new();
    let mut potential_donor_indices = Vec::new();
    
    for (atom_index, atom_reference) in protein_data.pdb.atoms().enumerate() {
        let element_symbol = atom_reference.element().map(|e| e.symbol()).unwrap_or("?");
        match element_symbol {
            "O" | "N" => {
                potential_acceptor_indices.push(atom_index);
                potential_donor_indices.push(atom_index);
            }
            _ => {}
        }
    }
    
    let atom_coordinates: Vec<Vec3> = protein_data.pdb.atoms()
        .map(|a| {
            let p = a.pos();
            Vec3::new(p.0 as f32, p.1 as f32, p.2 as f32)
        })
        .collect();
        
    // 2. Perform distance-based search (threshold < 3.5 A)
    // Note: A real implementation would use a spatial index (RTree)
    for &donor_index in &potential_donor_indices {
        let donor_position = atom_coordinates[donor_index];
        
        for &acceptor_index in &potential_acceptor_indices {
            if donor_index == acceptor_index { continue; }
            
            let acceptor_position = atom_coordinates[acceptor_index];
            let distance_squared = donor_position.distance_squared(acceptor_position);
            
            if distance_squared < 12.25 { // 3.5 A squared
                let distance = distance_squared.sqrt();
                identified_hydrogen_bonds.push(HydrogenBond {
                    donor_atom_index: donor_index,
                    acceptor_atom_index: acceptor_index,
                    distance_in_angstroms: distance,
                });
            }
        }
    }
    
    identified_hydrogen_bonds
}
