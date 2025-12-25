//! Molecular surface generation and management
//!
//! This module provides tools for calculating and rendering molecular surfaces,
//! such as Solvent-Excluded Surfaces (SES) using Marching Cubes.

pub mod marching_cubes;
pub mod distance_field;

/// Represents a vertex on the molecular surface mesh
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct SurfaceVertex {
    /// World space position of the vertex
    pub world_space_position: [f32; 3],
    /// Normal vector for lighting calculations
    pub surface_normal_vector: [f32; 3],
    /// Color of the surface at this vertex
    pub surface_color_rgb: [f32; 3],
}

/// A complete mesh representing the molecular surface
pub struct MolecularSurfaceMesh {
    /// List of vertices in the mesh
    pub surface_vertices_collection: Vec<SurfaceVertex>,
    /// List of indices for triangle drawing
    pub surface_indices_collection: Vec<u32>,
    /// Whether the surface is currently visible
    pub is_surface_visible: bool,
}

impl MolecularSurfaceMesh {
    pub fn new() -> Self {
        Self {
            surface_vertices_collection: Vec::new(),
            surface_indices_collection: Vec::new(),
            is_surface_visible: false,
        }
    }

    pub fn is_mesh_empty(&self) -> bool {
        self.surface_vertices_collection.is_empty()
    }
}