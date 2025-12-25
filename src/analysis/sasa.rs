//! Solvent Accessible Surface Area (SASA) calculation using Shrake-Rupley algorithm

use glam::Vec3;
use crate::protein::structure::ProteinData;

/// Calculates the Solvent Accessible Surface Area (SASA) of the protein
pub fn calculate_protein_solvent_accessible_surface_area(
    protein_data: &ProteinData,
    probe_radius_angstroms: f32,
    number_of_test_points: usize,
) -> f32 {
    use rayon::prelude::*;

    if let Some(pdb) = &protein_data.underlying_pdb_data {
        // 1. Pre-calculate atom positions and their VdW radii extended by probe radius
        let atom_data_collection: Vec<(Vec3, f32)> = pdb.atoms()
            .map(|atom_reference| {
            let position_tuple = atom_reference.pos();
            let element_symbol = atom_reference.element().map(|e| e.symbol()).unwrap_or("?");
            let vdw_radius = match element_symbol {
                "H" => 1.20,
                "C" => 1.70,
                "N" => 1.55,
                "O" => 1.52,
                "S" => 1.80,
                "P" => 1.80,
                _ => 1.50,
            };
            (
                Vec3::new(position_tuple.0 as f32, position_tuple.1 as f32, position_tuple.2 as f32),
                vdw_radius + probe_radius_angstroms
            )
        })
        .collect();

    // 2. Generate a set of test points uniformly distributed on a unit sphere
    let test_points_on_unit_sphere = generate_uniformly_distributed_points_on_sphere(number_of_test_points);

    // 3. For each atom, check how many of its test points are not covered by any other atom
    let total_accessible_surface_area: f32 = atom_data_collection.par_iter()
        .enumerate()
        .map(|(current_atom_index, &(atom_center_position, extended_radius))| {
            let mut accessible_points_count = 0;
            
            // Optimization: Find neighbors within a reasonable distance (max extended radius * 2)
            let neighboring_atoms: Vec<_> = atom_data_collection.iter()
                .enumerate()
                .filter(|&(other_index, &(other_center, other_radius))| {
                    if other_index == current_atom_index { return false; }
                    let distance_squared = atom_center_position.distance_squared(other_center);
                    let max_distance = extended_radius + other_radius;
                    distance_squared < max_distance * max_distance
                })
                .map(|(_, data)| data)
                .collect();

            for test_point_vector in &test_points_on_unit_sphere {
                let actual_test_point_position = atom_center_position + (*test_point_vector * extended_radius);
                let mut is_point_buried = false;

                for (neighbor_center, neighbor_extended_radius) in &neighboring_atoms {
                    if actual_test_point_position.distance_squared(*neighbor_center) < neighbor_extended_radius * neighbor_extended_radius {
                        is_point_buried = true;
                        break;
                    }
                }

                if !is_point_buried {
                    accessible_points_count += 1;
                }
            }

            let atom_surface_area = 4.0 * std::f32::consts::PI * extended_radius * extended_radius;
            (accessible_points_count as f32 / number_of_test_points as f32) * atom_surface_area
        })
        .sum();

        total_accessible_surface_area
    } else {
        0.0
    }
}

/// Generates a set of points uniformly distributed on a unit sphere using the Fibonacci spiral method
fn generate_uniformly_distributed_points_on_sphere(number_of_points: usize) -> Vec<Vec3> {
    let mut points_collection = Vec::with_capacity(number_of_points);
    let golden_ratio_phi = std::f32::consts::PI * (3.0 - (5.0f32).sqrt()); // Golden angle in radians

    for point_index in 0..number_of_points {
        let vertical_coordinate_y = 1.0 - (point_index as f32 / (number_of_points - 1) as f32) * 2.0; // y goes from 1 to -1
        let radius_at_y = (1.0 - vertical_coordinate_y * vertical_coordinate_y).sqrt(); // radius at y

        let azimuthal_angle_theta = golden_ratio_phi * point_index as f32; // golden angle increment

        let horizontal_coordinate_x = azimuthal_angle_theta.cos() * radius_at_y;
        let horizontal_coordinate_z = azimuthal_angle_theta.sin() * radius_at_y;

        points_collection.push(Vec3::new(horizontal_coordinate_x, vertical_coordinate_y, horizontal_coordinate_z));
    }
    points_collection
}
