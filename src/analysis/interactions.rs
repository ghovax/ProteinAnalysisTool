//! Detection of physical and chemical interactions in proteins

use glam::Vec3;
use crate::protein::structure::ProteinData;
use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain};

/// Represents a detected salt bridge interaction
#[derive(Debug, Clone, Copy)]
pub struct SaltBridgeInteraction {
    pub donor_atom_index: usize,
    pub acceptor_atom_index: usize,
    pub inter_atomic_distance: f32,
}

/// Represents a detected pi-stacking interaction between aromatic rings
#[derive(Debug, Clone, Copy)]
pub struct PiStackingInteraction {
    pub first_ring_center: Vec3,
    pub second_ring_center: Vec3,
    pub inter_ring_distance: f32,
    pub relative_angle_degrees: f32,
}

/// Identifies salt bridges based on proximity of charged sidechain atoms
pub fn detect_salt_bridge_interactions(
    protein_data: &ProteinData,
    distance_threshold_angstroms: f32,
) -> Vec<SaltBridgeInteraction> {
    let mut identified_salt_bridges = Vec::new();

    // Positively charged atoms (Arg NH1/NH2, Lys NZ)
    let positively_charged_atoms: Vec<_> = protein_data.pdb.atoms_with_hierarchy()
        .filter(|h| {
            let res_name = h.residue().name().unwrap_or("");
            let atom_name = h.atom().name();
            (res_name == "ARG" && (atom_name == "NH1" || atom_name == "NH2")) ||
            (res_name == "LYS" && atom_name == "NZ")
        })
        .map(|h| (h.atom().pos(), h.atom().serial_number())) // Use serial number as index proxy if needed
        .collect();

    // Negatively charged atoms (Asp OD1/OD2, Glu OE1/OE2)
    let negatively_charged_atoms: Vec<_> = protein_data.pdb.atoms_with_hierarchy()
        .filter(|h| {
            let res_name = h.residue().name().unwrap_or("");
            let atom_name = h.atom().name();
            (res_name == "ASP" && (atom_name == "OD1" || atom_name == "OD2")) ||
            (res_name == "GLU" && (atom_name == "OE1" || atom_name == "OE2"))
        })
        .map(|h| (h.atom().pos(), h.atom().serial_number()))
        .collect();

    let threshold_squared = distance_threshold_angstroms * distance_threshold_angstroms;

    for (pos_coord, pos_idx) in &positively_charged_atoms {
        let pos_vec = Vec3::new(pos_coord.0 as f32, pos_coord.1 as f32, pos_coord.2 as f32);
        for (neg_coord, neg_idx) in &negatively_charged_atoms {
            let neg_vec = Vec3::new(neg_coord.0 as f32, neg_coord.1 as f32, neg_coord.2 as f32);
            let dist_sq = pos_vec.distance_squared(neg_vec);
            if dist_sq < threshold_squared {
                identified_salt_bridges.push(SaltBridgeInteraction {
                    donor_atom_index: *pos_idx as usize,
                    acceptor_atom_index: *neg_idx as usize,
                    inter_atomic_distance: dist_sq.sqrt(),
                });
            }
        }
    }

    identified_salt_bridges
}

/// Identifies pi-stacking interactions between aromatic residues (Phe, Tyr, Trp, His)
pub fn detect_pi_stacking_interactions(
    _protein_data: &ProteinData,
    _distance_threshold: f32,
) -> Vec<PiStackingInteraction> {
    // This would involve finding aromatic rings, calculating their centers and normals,
    // and checking for face-to-face or edge-to-face geometry.
    // For now, returning an empty set as a placeholder for this complex geometry task.
    Vec::new()
}
