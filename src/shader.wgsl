struct Uniforms {
    view_proj: mat4x4<f32>,
    camera_pos: vec3<f32>,
    _padding: f32,
}

@group(0) @binding(0)
var<uniform> uniforms: Uniforms;

// ============================================
// SPHERE RENDERING (Billboard quads)
// ============================================

struct SphereInput {
    @location(0) instance_pos: vec3<f32>,
    @location(1) radius: f32,
    @location(2) color: vec3<f32>,
    @builtin(vertex_index) vertex_index: u32,
}

struct SphereOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec3<f32>,
    @location(1) local_pos: vec2<f32>,
    @location(2) world_pos: vec3<f32>,
    @location(3) sphere_center: vec3<f32>,
    @location(4) radius: f32,
}

@vertex
fn vs_sphere(input: SphereInput) -> SphereOutput {
    var out: SphereOutput;

    // Billboard quad vertices (triangle strip: 0,1,2,3)
    var local: vec2<f32>;
    switch input.vertex_index {
        case 0u: { local = vec2<f32>(-1.0, -1.0); }
        case 1u: { local = vec2<f32>(1.0, -1.0); }
        case 2u: { local = vec2<f32>(-1.0, 1.0); }
        case 3u: { local = vec2<f32>(1.0, 1.0); }
        default: { local = vec2<f32>(0.0, 0.0); }
    }
    out.local_pos = local;

    // Calculate billboard axes (camera-facing)
    let to_camera = normalize(uniforms.camera_pos - input.instance_pos);
    let right = normalize(cross(vec3<f32>(0.0, 1.0, 0.0), to_camera));
    let up = cross(to_camera, right);

    // Offset from sphere center
    let offset = (right * local.x + up * local.y) * input.radius;
    let world_pos = input.instance_pos + offset;

    out.clip_position = uniforms.view_proj * vec4<f32>(world_pos, 1.0);
    out.color = input.color;
    out.world_pos = world_pos;
    out.sphere_center = input.instance_pos;
    out.radius = input.radius;

    return out;
}

@fragment
fn fs_sphere(input: SphereOutput) -> @location(0) vec4<f32> {
    // Discard pixels outside the sphere (circle in 2D)
    let dist_sq = input.local_pos.x * input.local_pos.x + input.local_pos.y * input.local_pos.y;
    if dist_sq > 1.0 {
        discard;
    }

    // Calculate sphere normal for lighting
    let z = sqrt(1.0 - dist_sq);
    let normal = vec3<f32>(input.local_pos.x, input.local_pos.y, z);

    // Simple lighting
    let light_dir = normalize(vec3<f32>(0.5, 1.0, 0.8));
    let ambient = 0.3;
    let diffuse = max(dot(normal, light_dir), 0.0) * 0.7;

    // Specular highlight
    let view_dir = normalize(uniforms.camera_pos - input.world_pos);
    let reflect_dir = reflect(-light_dir, normal);
    let specular = pow(max(dot(view_dir, reflect_dir), 0.0), 32.0) * 0.3;

    let lighting = ambient + diffuse + specular;
    let final_color = input.color * lighting;

    return vec4<f32>(final_color, 1.0);
}

// ============================================
// LINE RENDERING (Backbone trace)
// ============================================

struct LineInput {
    @location(0) position: vec3<f32>,
    @location(1) color: vec3<f32>,
}

struct LineOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec3<f32>,
}

@vertex
fn vs_line(input: LineInput) -> LineOutput {
    var out: LineOutput;
    out.clip_position = uniforms.view_proj * vec4<f32>(input.position, 1.0);
    out.color = input.color;
    return out;
}

@fragment
fn fs_line(input: LineOutput) -> @location(0) vec4<f32> {
    return vec4<f32>(input.color, 1.0);
}
