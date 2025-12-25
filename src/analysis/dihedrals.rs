//! Calculation of protein dihedral angles (phi, psi, omega)

use glam::Vec3;

/// Calculates the dihedral angle (torsion angle) between four points in degrees
pub fn calculate_dihedral_angle_between_points(
    first_point: Vec3,
    second_point: Vec3,
    third_point: Vec3,
    fourth_point: Vec3,
) -> f32 {
    let vector_from_first_to_second = second_point - first_point;
    let vector_from_second_to_third = third_point - second_point;
    let vector_from_third_to_fourth = fourth_point - third_point;

    let normal_vector_of_first_plane = vector_from_first_to_second.cross(vector_from_second_to_third);
    let normal_vector_of_second_plane = vector_from_second_to_third.cross(vector_from_third_to_fourth);

    let orthogonal_component_vector = normal_vector_of_first_plane.cross(normal_vector_of_second_plane);
    
    let sine_of_dihedral_angle = orthogonal_component_vector.dot(vector_from_second_to_third.normalize());
    let cosine_of_dihedral_angle = normal_vector_of_first_plane.dot(normal_vector_of_second_plane);

    sine_of_dihedral_angle.atan2(cosine_of_dihedral_angle).to_degrees()
}

/// Represents the backbone dihedral angles for a single residue
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BackboneDihedralAngles {
    pub phi_angle: Option<f32>,
    pub psi_angle: Option<f32>,
    pub omega_angle: Option<f32>,
}

impl BackboneDihedralAngles {
    pub fn new() -> Self {
        Self {
            phi_angle: None,
            psi_angle: None,
            omega_angle: None,
        }
    }
}

/// Calculates all backbone dihedral angles for a protein structure in parallel
pub fn calculate_all_backbone_dihedrals(
    protein_data_reference: &crate::protein::structure::ProteinData,
) -> Vec<(String, isize, BackboneDihedralAngles)> {
    use rayon::prelude::*;

    let mut calculated_dihedral_results_collection = Vec::new();

    for current_chain_reference in protein_data_reference.pdb.chains() {
        let chain_identifier_string = current_chain_reference.id().to_string();
        let residues_in_current_chain: Vec<_> = current_chain_reference.residues().collect();
        
        // Parallelize over residues within the chain
        let chain_results: Vec<(String, isize, BackboneDihedralAngles)> = (0..residues_in_current_chain.len())
            .into_par_iter()
            .map(|residue_list_index| {
                let current_residue = residues_in_current_chain[residue_list_index];
                let mut angles_for_current_residue = BackboneDihedralAngles::new();
                
                // Helper to find atom positions
                let find_atom_world_position = |target_residue: &pdbtbx::Residue, target_atom_name: &str| -> Option<Vec3> {
                    target_residue.conformers().next()?.atoms()
                        .find(|current_atom_reference| current_atom_reference.name() == target_atom_name)
                        .map(|current_atom_reference| {
                            let atom_coordinates_tuple = current_atom_reference.pos();
                            Vec3::new(
                                atom_coordinates_tuple.0 as f32,
                                atom_coordinates_tuple.1 as f32,
                                atom_coordinates_tuple.2 as f32
                            )
                        })
                };

                let nitrogen_atom_position = find_atom_world_position(current_residue, "N");
                let alpha_carbon_atom_position = find_atom_world_position(current_residue, "CA");
                let carbonyl_carbon_atom_position = find_atom_world_position(current_residue, "C");

                if let (Some(n_pos), Some(ca_pos), Some(c_pos)) = (nitrogen_atom_position, alpha_carbon_atom_position, carbonyl_carbon_atom_position) {
                    // Phi: C(i-1) - N(i) - CA(i) - C(i)
                    if residue_list_index > 0 {
                        if let Some(previous_carbonyl_carbon_position) = find_atom_world_position(residues_in_current_chain[residue_list_index - 1], "C") {
                            angles_for_current_residue.phi_angle = Some(calculate_dihedral_angle_between_points(previous_carbonyl_carbon_position, n_pos, ca_pos, c_pos));
                        }
                    }

                    // Psi: N(i) - CA(i) - C(i) - N(i+1)
                    if residue_list_index < residues_in_current_chain.len() - 1 {
                        if let Some(next_nitrogen_atom_position) = find_atom_world_position(residues_in_current_chain[residue_list_index + 1], "N") {
                            angles_for_current_residue.psi_angle = Some(calculate_dihedral_angle_between_points(n_pos, ca_pos, c_pos, next_nitrogen_atom_position));
                        }
                    }

                    // Omega: CA(i) - C(i) - N(i+1) - CA(i+1)
                    if residue_list_index < residues_in_current_chain.len() - 1 {
                        let next_residue_reference = residues_in_current_chain[residue_list_index + 1];
                        if let (Some(next_nitrogen_atom_position), Some(next_alpha_carbon_atom_position)) = (find_atom_world_position(next_residue_reference, "N"), find_atom_world_position(next_residue_reference, "CA")) {
                            angles_for_current_residue.omega_angle = Some(calculate_dihedral_angle_between_points(ca_pos, c_pos, next_nitrogen_atom_position, next_alpha_carbon_atom_position));
                        }
                    }
                }

                (chain_identifier_string.clone(), current_residue.serial_number(), angles_for_current_residue)
            })
            .collect();
        
        calculated_dihedral_results_collection.extend(chain_results);
    }

    calculated_dihedral_results_collection
}
