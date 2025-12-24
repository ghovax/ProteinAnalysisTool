use glam::{Mat4, Vec3};

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
        let x = self.distance * self.pitch.cos() * self.yaw.sin();
        let y = self.distance * self.pitch.sin();
        let z = self.distance * self.pitch.cos() * self.yaw.cos();
        self.target + Vec3::new(x, y, z)
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
}
