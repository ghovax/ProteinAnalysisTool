use glam::{Mat4, Vec3, Vec2};

pub struct Camera {
    pub target: Vec3,
    pub distance: f32,
    pub yaw: f32,   // Horizontal rotation (radians)
    pub pitch: f32, // Vertical rotation (radians)
    pub fov: f32,
    pub aspect: f32,
    pub near: f32,
    pub far: f32,
}

impl Camera {
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

    pub fn position(&self) -> Vec3 {
        let camera_position_x = self.distance * self.pitch.cos() * self.yaw.sin();
        let camera_position_y = self.distance * self.pitch.sin();
        let camera_position_z = self.distance * self.pitch.cos() * self.yaw.cos();
        self.target + Vec3::new(camera_position_x, camera_position_y, camera_position_z)
    }

    pub fn view_matrix(&self) -> Mat4 {
        Mat4::look_at_rh(self.position(), self.target, Vec3::Y)
    }

    pub fn projection_matrix(&self) -> Mat4 {
        Mat4::perspective_rh(self.fov, self.aspect, self.near, self.far)
    }

    pub fn view_projection_matrix(&self) -> Mat4 {
        self.projection_matrix() * self.view_matrix()
    }

    pub fn rotate(&mut self, delta_yaw: f32, delta_pitch: f32) {
        self.yaw += delta_yaw;
        self.pitch = (self.pitch + delta_pitch).clamp(-1.5, 1.5);
    }

    pub fn zoom(&mut self, delta: f32) {
        self.distance = (self.distance * (1.0 - delta * 0.1)).clamp(1.0, 500.0);
    }

    pub fn set_aspect(&mut self, aspect: f32) {
        self.aspect = aspect;
    }

    pub fn focus_on(&mut self, center: Vec3, radius: f32) {
        self.target = center;
        self.distance = radius * 2.5;
    }

    pub fn ray_from_screen(&self, screen_coordinates: Vec2, screen_dimensions: Vec2) -> (Vec3, Vec3) {
        let normalized_device_coordinates = Vec2::new(
            (2.0 * screen_coordinates.x / screen_dimensions.x) - 1.0,
            1.0 - (2.0 * screen_coordinates.y / screen_dimensions.y),
        );
        let inverse_view_projection_matrix = self.view_projection_matrix().inverse();
        let near_plane_world_point = inverse_view_projection_matrix.project_point3(Vec3::new(normalized_device_coordinates.x, normalized_device_coordinates.y, -1.0));
        let far_plane_world_point = inverse_view_projection_matrix.project_point3(Vec3::new(normalized_device_coordinates.x, normalized_device_coordinates.y, 1.0));
        let ray_direction_vector = (far_plane_world_point - near_plane_world_point).normalize();
        (near_plane_world_point, ray_direction_vector)
    }
}