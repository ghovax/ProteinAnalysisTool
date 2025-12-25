//! Interactive 3D camera for protein visualization
//!
//! This module provides a spherical coordinate camera system that supports
//! orbiting, zooming, and focusing on specific points in 3D space

use glam::{Mat4, Vec2, Vec3};
use serde::{Deserialize, Serialize};

/// A 3D camera using spherical coordinates
#[derive(Serialize, Deserialize, Clone, Copy)]
pub struct Camera {
    /// The point the camera is looking at
    pub target: Vec3,
    /// Distance from the target
    pub distance: f32,
    /// Horizontal rotation in radians
    pub yaw: f32,
    /// Vertical rotation in radians
    pub pitch: f32,
    /// Field of view in radians
    pub fov: f32,
    /// Aspect ratio (width / height)
    pub aspect: f32,
    /// Near clipping plane
    pub near: f32,
    /// Far clipping plane
    pub far: f32,
}

impl Camera {
    /// Creates a new camera with default orientation and the specified aspect ratio
    pub fn new(aspect: f32) -> Self {
        Self {
            target: Vec3::ZERO,
            distance: 50.0,
            yaw: 0.0,
            pitch: 0.3,
            fov: 45.0_f32.to_radians(),
            aspect,
            near: 0.1,
            far: 1000.0,
        }
    }

    /// Calculates the world-space position of the camera
    pub fn position(&self) -> Vec3 {
        let camera_position_x = self.distance * self.pitch.cos() * self.yaw.sin();
        let camera_position_y = self.distance * self.pitch.sin();
        let camera_position_z = self.distance * self.pitch.cos() * self.yaw.cos();
        self.target + Vec3::new(camera_position_x, camera_position_y, camera_position_z)
    }

    /// Returns the view matrix
    pub fn view_matrix(&self) -> Mat4 {
        Mat4::look_at_rh(self.position(), self.target, Vec3::Y)
    }

    /// Returns the perspective projection matrix
    pub fn projection_matrix(&self) -> Mat4 {
        Mat4::perspective_rh(self.fov, self.aspect, self.near, self.far)
    }

    /// Returns the combined view-projection matrix
    pub fn view_projection_matrix(&self) -> Mat4 {
        self.projection_matrix() * self.view_matrix()
    }

    /// Rotates the camera by the given yaw and pitch deltas
    pub fn rotate(&mut self, delta_yaw: f32, delta_pitch: f32) {
        self.yaw += delta_yaw;
        self.pitch = (self.pitch + delta_pitch).clamp(-1.5, 1.5);
    }

    /// Adjusts the camera distance (zoom) based on the given delta
    pub fn zoom(&mut self, delta: f32) {
        self.distance = (self.distance * (1.0 - delta * 0.1)).clamp(1.0, 500.0);
    }

    /// Sets the camera's aspect ratio
    pub fn set_aspect(&mut self, aspect: f32) {
        self.aspect = aspect;
    }

    /// Focuses the camera on a specific center point and adjusts distance based on radius
    pub fn focus_on(&mut self, center: Vec3, radius: f32) {
        self.target = center;
        self.distance = radius * 2.5;
    }

    /// Translates the camera target position along the camera's local axes
    pub fn translate_camera_target_position(&mut self, horizontal_delta: f32, vertical_delta: f32) {
        let rotation_matrix = Mat4::from_euler(glam::EulerRot::YXZ, self.yaw, self.pitch, 0.0);
        let camera_right_vector = rotation_matrix.transform_vector3(Vec3::X);
        let camera_up_vector = rotation_matrix.transform_vector3(Vec3::Y);

        // Scale translation speed by distance to target for consistent feel
        let translation_speed_scaling_factor = self.distance * 0.0005;
        self.target += camera_right_vector * -horizontal_delta * translation_speed_scaling_factor;
        self.target += camera_up_vector * vertical_delta * translation_speed_scaling_factor;
    }

    /// Generates a ray in world space from screen coordinates
    pub fn calculate_ray_from_screen_coordinates(
        &self,
        screen_coordinates: Vec2,
        viewport_dimensions: Vec2,
    ) -> (Vec3, Vec3) {
        let normalized_device_coordinates = Vec2::new(
            (2.0 * screen_coordinates.x / viewport_dimensions.x) - 1.0,
            1.0 - (2.0 * screen_coordinates.y / viewport_dimensions.y),
        );
        let inverse_view_projection_matrix = self.view_projection_matrix().inverse();
        
        // glam's Mat4::perspective_rh maps the near plane to -1.0 and far to 1.0 (GL convention)
        let near_plane_world_point = inverse_view_projection_matrix.project_point3(Vec3::new(
            normalized_device_coordinates.x,
            normalized_device_coordinates.y,
            -1.0,
        ));
        let far_plane_world_point = inverse_view_projection_matrix.project_point3(Vec3::new(
            normalized_device_coordinates.x,
            normalized_device_coordinates.y,
            1.0,
        ));
        
        let ray_direction_unit_vector = (far_plane_world_point - near_plane_world_point).normalize();
        (near_plane_world_point, ray_direction_unit_vector)
    }
}
