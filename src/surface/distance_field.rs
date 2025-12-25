//! Parallelized distance field calculation for molecular surfaces

use glam::Vec3;
use rayon::prelude::*;
use crate::protein::structure::ProteinData;

/// A 3D grid of distance values for isosurface extraction
pub struct DistanceFieldGrid {
    /// Flattened array of scalar values representing the distance field
    pub distance_values_collection: Vec<f32>,
    /// Number of voxels along each axis (X, Y, Z)
    pub grid_dimensions: [usize; 3],
    /// World space position of the grid's minimum corner
    pub grid_origin_position: Vec3,
    /// Edge length of a single cubic voxel in Angstroms
    pub voxel_size_in_angstroms: f32,
}

impl DistanceFieldGrid {
    pub fn new(
        origin_position: Vec3,
        dimensions: [usize; 3],
        voxel_size: f32,
    ) -> Self {
        let total_voxels_count = dimensions[0] * dimensions[1] * dimensions[2];
        Self {
            distance_values_collection: vec![0.0; total_voxels_count],
            grid_dimensions: dimensions,
            grid_origin_position: origin_position,
            voxel_size_in_angstroms: voxel_size,
        }
    }

    pub fn get_distance_value_at_coordinates(
        &self,
        coordinate_x: usize,
        coordinate_y: usize,
        coordinate_z: usize,
    ) -> f32 {
        let linear_index = coordinate_z * self.grid_dimensions[0] * self.grid_dimensions[1]
            + coordinate_y * self.grid_dimensions[0]
            + coordinate_x;
        self.distance_values_collection[linear_index]
    }

    pub fn calculate_world_position_for_voxel(
        &self,
        coordinate_x: usize,
        coordinate_y: usize,
        coordinate_z: usize,
    ) -> Vec3 {
        self.grid_origin_position
            + Vec3::new(
                coordinate_x as f32,
                coordinate_y as f32,
                coordinate_z as f32,
            ) * self.voxel_size_in_angstroms
    }
}

/// Calculates a Signed Distance Field for the protein using parallel processing
pub fn calculate_solvent_accessible_distance_field(
    protein_data_reference: &ProteinData,
    desired_voxel_size: f32,
    grid_padding_in_angstroms: f32,
) -> DistanceFieldGrid {
    let (minimum_bounding_bound, maximum_bounding_bound) = protein_data_reference.bounding_box();
    let padded_grid_minimum = minimum_bounding_bound - Vec3::splat(grid_padding_in_angstroms);
    let padded_grid_maximum = maximum_bounding_bound + Vec3::splat(grid_padding_in_angstroms);
    let total_grid_size_vector = padded_grid_maximum - padded_grid_minimum;

    let grid_dimension_x = (total_grid_size_vector.x / desired_voxel_size).ceil() as usize + 1;
    let grid_dimension_y = (total_grid_size_vector.y / desired_voxel_size).ceil() as usize + 1;
    let grid_dimension_z = (total_grid_size_vector.z / desired_voxel_size).ceil() as usize + 1;

    // Pre-collect atom positions and their associated Van der Waals radii
    let mut atom_world_positions_collection = Vec::new();
    let mut atom_vdw_radii_collection = Vec::new();

    for current_atom_reference in protein_data_reference.pdb.atoms() {
        let atom_position_tuple = current_atom_reference.pos();
        atom_world_positions_collection.push(Vec3::new(
            atom_position_tuple.0 as f32,
            atom_position_tuple.1 as f32,
            atom_position_tuple.2 as f32,
        ));

        let current_atom_element_symbol = current_atom_reference
            .element()
            .map(|element_reference| element_reference.symbol())
            .unwrap_or("?");

        let vdw_radius_value = match current_atom_element_symbol {
            "H" => 1.20,
            "C" => 1.70,
            "N" => 1.55,
            "O" => 1.52,
            "S" => 1.80,
            "P" => 1.80,
            _ => 1.50,
        };
        atom_vdw_radii_collection.push(vdw_radius_value);
    }

    // Use Rayon to calculate slices of the grid in parallel
    let total_elements_per_z_slice = grid_dimension_x * grid_dimension_y;
    let mut flattened_distance_values = vec![0.0f32; grid_dimension_z * total_elements_per_z_slice];

    flattened_distance_values
        .par_chunks_mut(total_elements_per_z_slice)
        .enumerate()
        .for_each(|(coordinate_z, current_z_slice)| {
            for coordinate_y in 0..grid_dimension_y {
                for coordinate_x in 0..grid_dimension_x {
                    let current_voxel_world_position = padded_grid_minimum
                        + Vec3::new(
                            coordinate_x as f32,
                            coordinate_y as f32,
                            coordinate_z as f32,
                        ) * desired_voxel_size;
                    
                    let mut minimum_distance_to_surface = f32::MAX;

                    for (atom_index, atom_position_vector) in
                        atom_world_positions_collection.iter().enumerate()
                    {
                        let distance_to_atom_center =
                            current_voxel_world_position.distance(*atom_position_vector);
                        let distance_to_atom_surface =
                            distance_to_atom_center - atom_vdw_radii_collection[atom_index];

                        if distance_to_atom_surface < minimum_distance_to_surface {
                            minimum_distance_to_surface = distance_to_atom_surface;
                        }
                    }
                    
                    current_z_slice[coordinate_y * grid_dimension_x + coordinate_x] = minimum_distance_to_surface;
                }
            }
        });

    DistanceFieldGrid {
        distance_values_collection: flattened_distance_values,
        grid_dimensions: [grid_dimension_x, grid_dimension_y, grid_dimension_z],
        grid_origin_position: padded_grid_minimum,
        voxel_size_in_angstroms: desired_voxel_size,
    }
}
