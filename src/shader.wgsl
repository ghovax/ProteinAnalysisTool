struct Uniforms {
    view_projection_matrix: mat4x4<f32>,
    camera_world_position: vec3<f32>,
    _padding: f32,
}

@group(0) @binding(0)
var<uniform> global_uniforms: Uniforms;

// Sphere rendering (billboard quads)

struct SphereInput {
    @location(0) instance_world_position: vec3<f32>,
    @location(1) radius: f32,
    @location(2) color: vec3<f32>,
    @location(3) selection_factor: f32,
    @builtin(vertex_index) builtin_vertex_index: u32,
}

struct SphereOutput {
    @builtin(position) clip_space_position: vec4<f32>,
    @location(0) color: vec3<f32>,
    @location(1) local_billboard_position: vec2<f32>,
    @location(2) world_space_position: vec3<f32>,
    @location(3) sphere_center_position: vec3<f32>,
    @location(4) radius: f32,
    @location(5) selection_factor: f32,
}

@vertex
fn vs_sphere(input: SphereInput) -> SphereOutput {
    var shader_output_data: SphereOutput;

    // Billboard quad vertices (triangle strip: 0,1,2,3)
    var local_vertex_coordinates: vec2<f32>;
    switch input.builtin_vertex_index {
        case 0u: { local_vertex_coordinates = vec2<f32>(-1.0, -1.0); }
        case 1u: { local_vertex_coordinates = vec2<f32>(1.0, -1.0); }
        case 2u: { local_vertex_coordinates = vec2<f32>(-1.0, 1.0); }
        case 3u: { local_vertex_coordinates = vec2<f32>(1.0, 1.0); }
        default: { local_vertex_coordinates = vec2<f32>(0.0, 0.0); }
    }
    shader_output_data.local_billboard_position = local_vertex_coordinates;

    // Calculate billboard axes (camera-facing)
    let direction_to_camera = normalize(global_uniforms.camera_world_position - input.instance_world_position);
    let billboard_right_vector = normalize(cross(vec3<f32>(0.0, 1.0, 0.0), direction_to_camera));
    let billboard_up_vector = cross(direction_to_camera, billboard_right_vector);

    // Offset from sphere center
    let vertex_world_offset = (billboard_right_vector * local_vertex_coordinates.x + billboard_up_vector * local_vertex_coordinates.y) * input.radius;
    let world_space_position = input.instance_world_position + vertex_world_offset;

    shader_output_data.clip_space_position = global_uniforms.view_projection_matrix * vec4<f32>(world_space_position, 1.0);
    shader_output_data.color = input.color;
    shader_output_data.world_space_position = world_space_position;
    shader_output_data.sphere_center_position = input.instance_world_position;
    shader_output_data.radius = input.radius;
    shader_output_data.selection_factor = input.selection_factor;

    return shader_output_data;
}

@fragment
fn fs_sphere(input: SphereOutput) -> @location(0) vec4<f32> {
    // Discard pixels outside the sphere (circle in 2D)
    let squared_distance_from_center = input.local_billboard_position.x * input.local_billboard_position.x + input.local_billboard_position.y * input.local_billboard_position.y;
    if squared_distance_from_center > 1.0 {
        discard;
    }

    // Calculate sphere normal for lighting
    let surface_z_coordinate = sqrt(1.0 - squared_distance_from_center);
    let fragment_surface_normal = vec3<f32>(input.local_billboard_position.x, input.local_billboard_position.y, surface_z_coordinate);

    // Simple lighting
    let light_source_direction = normalize(vec3<f32>(0.5, 1.0, 0.8));
    let ambient_lighting_factor = 0.3;
    let diffuse_lighting_factor = max(dot(fragment_surface_normal, light_source_direction), 0.0) * 0.7;

    // Specular highlight
    let view_direction_vector = normalize(global_uniforms.camera_world_position - input.world_space_position);
    let reflection_direction_vector = reflect(-light_source_direction, fragment_surface_normal);
    let specular_lighting_factor = pow(max(dot(view_direction_vector, reflection_direction_vector), 0.0), 32.0) * 0.3;

    let total_lighting_intensity = ambient_lighting_factor + diffuse_lighting_factor + specular_lighting_factor;
    var final_fragment_color = input.color * total_lighting_intensity;

    // Highlight selected atoms
    if input.selection_factor > 0.5 {
        // Add a white/yellow tint and an outline-like effect at the edges
        let selection_tint = vec3<f32>(1.0, 0.9, 0.4);
        final_fragment_color = mix(final_fragment_color, selection_tint, 0.4);
        
        // Edge highlight (simple outline)
        if squared_distance_from_center > 0.8 {
            final_fragment_color = selection_tint;
        }
    }

    return vec4<f32>(final_fragment_color, 1.0);
}

// Line rendering (backbone trace)

struct LineInput {
    @location(0) position: vec3<f32>,
    @location(1) color: vec3<f32>,
}

struct LineOutput {
    @builtin(position) clip_space_position: vec4<f32>,
    @location(0) color: vec3<f32>,
}

@vertex
fn vs_line(input: LineInput) -> LineOutput {
    var shader_output_data: LineOutput;
    shader_output_data.clip_space_position = global_uniforms.view_projection_matrix * vec4<f32>(input.position, 1.0);
    shader_output_data.color = input.color;
    return shader_output_data;
}

@fragment
fn fs_line(input: LineOutput) -> @location(0) vec4<f32> {
    return vec4<f32>(input.color, 1.0);
}
