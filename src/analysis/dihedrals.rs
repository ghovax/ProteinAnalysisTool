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
