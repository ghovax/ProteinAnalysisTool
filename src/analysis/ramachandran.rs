//! Calculation of Phi/Psi angles for Ramachandran plots

use crate::protein::structure::ProteinData;
use super::dihedrals::{calculate_dihedral_angle_between_points, BackboneDihedralAngles};

/// Represents a point on a Ramachandran plot
#[derive(Debug, Clone)]
pub struct RamachandranPoint {
    pub phi_angle: f32,
    pub psi_angle: f32,
    pub residue_name: String,
    pub residue_number: isize,
    pub chain_identifier: String,
}

/// Calculates Ramachandran (Phi/Psi) angles for all residues in a protein
pub fn calculate_ramachandran_angles_for_protein(
    protein_data: &ProteinData,
) -> Vec<RamachandranPoint> {
    let mut calculated_ramachandran_points = Vec::new();

    for current_chain_reference in protein_data.pdb.chains() {
        let chain_identifier = current_chain_reference.id().to_string();
        let residues_in_chain: Vec<_> = current_chain_reference.residues().collect();

        for residue_index in 0..residues_in_chain.len() {
            let current_residue = residues_in_chain[residue_index];
            
            // Need atoms: C(i-1), N(i), CA(i), C(i), N(i+1)
            let nitrogen_atom_position = find_atom_position_within_residue(current_residue, "N");
            let alpha_carbon_atom_position = find_atom_position_within_residue(current_residue, "CA");
            let carbonyl_carbon_atom_position = find_atom_position_within_residue(current_residue, "C");

            if let (Some(position_n), Some(position_ca), Some(position_c)) = (nitrogen_atom_position, alpha_carbon_atom_position, carbonyl_carbon_atom_position) {
                let mut backbone_dihedral_angles = BackboneDihedralAngles::new();

                // Phi: C(i-1) - N(i) - CA(i) - C(i)
                if residue_index > 0 {
                    let previous_residue = residues_in_chain[residue_index - 1];
                    if let Some(previous_carbonyl_carbon_position) = find_atom_position_within_residue(previous_residue, "C") {
                        backbone_dihedral_angles.phi_angle = Some(calculate_dihedral_angle_between_points(
                            previous_carbonyl_carbon_position, position_n, position_ca, position_c
                        ));
                    }
                }

                // Psi: N(i) - CA(i) - C(i) - N(i+1)
                if residue_index < residues_in_chain.len() - 1 {
                    let next_residue = residues_in_chain[residue_index + 1];
                    if let Some(next_nitrogen_atom_position) = find_atom_position_within_residue(next_residue, "N") {
                        backbone_dihedral_angles.psi_angle = Some(calculate_dihedral_angle_between_points(
                            position_n, position_ca, position_c, next_nitrogen_atom_position
                        ));
                    }
                }

                if let (Some(phi_angle_value), Some(psi_angle_value)) = (backbone_dihedral_angles.phi_angle, backbone_dihedral_angles.psi_angle) {
                    calculated_ramachandran_points.push(RamachandranPoint {
                        phi_angle: phi_angle_value,
                        psi_angle: psi_angle_value,
                        residue_name: current_residue.name().unwrap().to_string(),
                        residue_number: current_residue.serial_number(),
                        chain_identifier: chain_identifier.clone(),
                    });
                }
            }
        }
    }

    calculated_ramachandran_points
}

fn find_atom_position_within_residue(residue: &pdbtbx::Residue, target_atom_name: &str) -> Option<glam::Vec3> {
    for current_conformer in residue.conformers() {
        for current_atom in current_conformer.atoms() {
            if current_atom.name() == target_atom_name {
                let atom_position_tuple = current_atom.pos();
                return Some(glam::Vec3::new(atom_position_tuple.0 as f32, atom_position_tuple.1 as f32, atom_position_tuple.2 as f32));
            }
        }
    }
    None
}
