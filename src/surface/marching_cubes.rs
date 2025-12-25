//! Parallelized Marching Cubes algorithm for isosurface extraction

use glam::Vec3;
use rayon::prelude::*;
use super::SurfaceVertex;
use super::distance_field::DistanceFieldGrid;

// The user will provide these tables in the code
extern "C" {
    // Placeholder for where the tables will be provided/linked
    // Since we are in Rust, we just assume they are defined elsewhere or will be added
}

/// Linearly interpolates between two corner positions based on their scalar values
pub fn calculate_interpolated_vertex_position(
    position_vector_a: Vec3,
    position_vector_b: Vec3,
    scalar_value_a: f32,
    scalar_value_b: f32,
    target_threshold_value: f32,
) -> Vec3 {
    if (scalar_value_a - scalar_value_b).abs() < 1e-6 {
        return position_vector_a;
    }
    let interpolation_factor = (target_threshold_value - scalar_value_a) / (scalar_value_b - scalar_value_a);
    position_vector_a + interpolation_factor * (position_vector_b - position_vector_a)
}

/// Generates a triangular mesh from a distance field grid using parallel processing
pub fn extract_isosurface_from_distance_field(
    distance_field_grid: &DistanceFieldGrid,
    isosurface_threshold: f32,
    edge_intersection_lookup_table: &[u32; 256],
    triangulation_lookup_table: &[[i32; 16]; 256],
) -> (Vec<SurfaceVertex>, Vec<u32>) {
    let grid_dimension_x = distance_field_grid.grid_dimensions[0];
    let grid_dimension_y = distance_field_grid.grid_dimensions[1];
    let grid_dimension_z = distance_field_grid.grid_dimensions[2];

    // Process Z-slices in parallel
    let generated_mesh_fragments: Vec<(Vec<SurfaceVertex>, Vec<u32>)> = (0..(grid_dimension_z - 1))
        .into_par_iter()
        .map(|coordinate_z| {
            let mut slice_vertices = Vec::new();
            let mut slice_indices = Vec::new();

            for coordinate_y in 0..(grid_dimension_y - 1) {
                for coordinate_x in 0..(grid_dimension_x - 1) {
                    let mut cube_configuration_index = 0usize;
                    
                    let corner_indices = [
                        (coordinate_x, coordinate_y, coordinate_z),
                        (coordinate_x + 1, coordinate_y, coordinate_z),
                        (coordinate_x + 1, coordinate_y, coordinate_z + 1),
                        (coordinate_x, coordinate_y, coordinate_z + 1),
                        (coordinate_x, coordinate_y + 1, coordinate_z),
                        (coordinate_x + 1, coordinate_y + 1, coordinate_z),
                        (coordinate_x + 1, coordinate_y + 1, coordinate_z + 1),
                        (coordinate_x, coordinate_y + 1, coordinate_z + 1),
                    ];

                    let corner_values: [f32; 8] = corner_indices.map(|(x, y, z)| {
                        distance_field_grid.get_distance_value_at_coordinates(x, y, z)
                    });

                    for i in 0..8 {
                        if corner_values[i] < isosurface_threshold {
                            cube_configuration_index |= 1 << i;
                        }
                    }

                    let edges_to_intersect = edge_intersection_lookup_table[cube_configuration_index];
                    if edges_to_intersect == 0 {
                        continue;
                    }

                    let corner_positions: [Vec3; 8] = corner_indices.map(|(x, y, z)| {
                        distance_field_grid.calculate_world_position_for_voxel(x, y, z)
                    });

                    let mut edge_vertex_positions = [Vec3::ZERO; 12];
                    let edge_to_vertices = [
                        (0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4), (0,4), (1,5), (2,6), (3,7)
                    ];

                    for edge_index in 0..12 {
                        if (edges_to_intersect & (1 << edge_index)) != 0 {
                            let (vertex_a_index, vertex_b_index) = edge_to_vertices[edge_index];
                            edge_vertex_positions[edge_index] = calculate_interpolated_vertex_position(
                                corner_positions[vertex_a_index],
                                corner_positions[vertex_b_index],
                                corner_values[vertex_a_index],
                                corner_values[vertex_b_index],
                                isosurface_threshold,
                            );
                        }
                    }

                    let mut table_index = 0;
                    while triangulation_lookup_table[cube_configuration_index][table_index] != -1 {
                        let edge_a = triangulation_lookup_table[cube_configuration_index][table_index] as usize;
                        let edge_b = triangulation_lookup_table[cube_configuration_index][table_index + 1] as usize;
                        let edge_c = triangulation_lookup_table[cube_configuration_index][table_index + 2] as usize;

                        let base_index = slice_vertices.len() as u32;
                        
                        for &edge_idx in &[edge_a, edge_b, edge_c] {
                            slice_vertices.push(SurfaceVertex {
                                world_space_position: edge_vertex_positions[edge_idx].to_array(),
                                surface_normal_vector: [0.0, 1.0, 0.0], // Placeholder
                                surface_color_rgb: [0.7, 0.7, 0.8],
                            });
                        }
                        
                        slice_indices.push(base_index);
                        slice_indices.push(base_index + 1);
                        slice_indices.push(base_index + 2);

                        table_index += 3;
                    }
                }
            }
            (slice_vertices, slice_indices)
        })
        .collect();

    // Combine all fragments into a single mesh
    let mut final_vertices = Vec::new();
    let mut final_indices = Vec::new();

    for (mut fragment_vertices, fragment_indices) in generated_mesh_fragments {
        let index_offset = final_vertices.len() as u32;
        final_vertices.append(&mut fragment_vertices);
        for &idx in &fragment_indices {
            final_indices.push(idx + index_offset);
        }
    }

    (final_vertices, final_indices)
}
