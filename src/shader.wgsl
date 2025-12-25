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

// Cylinder rendering (for sticks and bonds)

struct CylinderInput {
    @location(0) start_position: vec3<f32>,
    @location(1) end_position: vec3<f32>,
    @location(2) radius: f32,
    @location(3) color: vec3<f32>,
    @location(4) selection_factor: f32,
    @builtin(vertex_index) builtin_vertex_index: u32,
}

struct CylinderOutput {
    @builtin(position) clip_space_position: vec4<f32>,
    @location(0) color: vec3<f32>,
    @location(1) local_billboard_coordinates: vec2<f32>,
    @location(2) world_space_fragment_position: vec3<f32>,
    @location(3) cylinder_axis_vector: vec3<f32>,
    @location(4) cylinder_start_point: vec3<f32>,
    @location(5) cylinder_radius: f32,
    @location(6) selection_factor: f32,
}

@vertex
fn vs_cylinder(instance_input: CylinderInput) -> CylinderOutput {
    var shader_output_data: CylinderOutput;

    let axis_vector = instance_input.end_position - instance_input.start_position;
    let cylinder_length = length(axis_vector);
    let normalized_axis_vector = axis_vector / cylinder_length;

    // Billboard quad indices: 0, 1, 2, 3
    var local_vertex_coordinates: vec2<f32>;
    switch instance_input.builtin_vertex_index {
        case 0u: { local_vertex_coordinates = vec2<f32>(-1.0, 0.0); }
        case 1u: { local_vertex_coordinates = vec2<f32>(1.0, 0.0); }
        case 2u: { local_vertex_coordinates = vec2<f32>(-1.0, 1.0); }
        case 3u: { local_vertex_coordinates = vec2<f32>(1.0, 1.0); }
        default: { local_vertex_coordinates = vec2<f32>(0.0, 0.0); }
    }
    
    // Calculate an orthogonal vector for the billboard width
    let camera_view_direction = normalize(global_uniforms.camera_world_position - instance_input.start_position);
    let billboard_width_vector = normalize(cross(normalized_axis_vector, camera_view_direction));

    // World position of the vertex
    let vertex_world_position = instance_input.start_position 
        + billboard_width_vector * local_vertex_coordinates.x * instance_input.radius 
        + normalized_axis_vector * local_vertex_coordinates.y * cylinder_length;

    shader_output_data.clip_space_position = global_uniforms.view_projection_matrix * vec4<f32>(vertex_world_position, 1.0);
    shader_output_data.color = instance_input.color;
    shader_output_data.local_billboard_coordinates = local_vertex_coordinates;
    shader_output_data.world_space_fragment_position = vertex_world_position;
    shader_output_data.cylinder_axis_vector = normalized_axis_vector;
    shader_output_data.cylinder_start_point = instance_input.start_position;
    shader_output_data.cylinder_radius = instance_input.radius;
    shader_output_data.selection_factor = instance_input.selection_factor;

    return shader_output_data;
}

@fragment
fn fs_cylinder(input_fragment: CylinderOutput) -> @location(0) vec4<f32> {
    // Simple Ray-Cylinder intersection in fragment shader for perfect rounded appearance
    let camera_to_fragment_vector = normalize(input_fragment.world_space_fragment_position - global_uniforms.camera_world_position);
    let ray_origin_relative_to_start = global_uniforms.camera_world_position - input_fragment.cylinder_start_point;
    
    // Project ray onto plane orthogonal to cylinder axis
    let ray_direction_orthogonal = camera_to_fragment_vector - input_fragment.cylinder_axis_vector * dot(camera_to_fragment_vector, input_fragment.cylinder_axis_vector);
    let origin_orthogonal = ray_origin_relative_to_start - input_fragment.cylinder_axis_vector * dot(ray_origin_relative_to_start, input_fragment.cylinder_axis_vector);
    
    // Solve quadratic: |origin_orthogonal + t * ray_direction_orthogonal|^2 = radius^2
    let quadratic_a = dot(ray_direction_orthogonal, ray_direction_orthogonal);
    let quadratic_b = 2.0 * dot(origin_orthogonal, ray_direction_orthogonal);
    let quadratic_c = dot(origin_orthogonal, origin_orthogonal) - input_fragment.cylinder_radius * input_fragment.cylinder_radius;
    
    let discriminant_value = quadratic_b * quadratic_b - 4.0 * quadratic_a * quadratic_c;
    if discriminant_value < 0.0 {
        discard;
    }
    
    let intersection_distance_t = (-quadratic_b - sqrt(discriminant_value)) / (2.0 * quadratic_a);
    let intersection_world_position = global_uniforms.camera_world_position + camera_to_fragment_vector * intersection_distance_t;
    
    // Check if the intersection is within the cylinder caps
    let projection_along_axis = dot(intersection_world_position - input_fragment.cylinder_start_point, input_fragment.cylinder_axis_vector);
    // Note: We don't discard based on length here because stick segments are capped by atoms usually
    
    // Calculate normal for lighting
    let point_on_axis = input_fragment.cylinder_start_point + input_fragment.cylinder_axis_vector * projection_along_axis;
    let surface_normal_vector = normalize(intersection_world_position - point_on_axis);
    
    // Lighting calculation
    let light_source_direction = normalize(vec3<f32>(0.5, 1.0, 0.8));
    let diffuse_lighting_intensity = max(dot(surface_normal_vector, light_source_direction), 0.0) * 0.7 + 0.3;
    
    var final_fragment_color = input_fragment.color * diffuse_lighting_intensity;
    
    if input_fragment.selection_factor > 0.5 {
        let selection_highlight_color = vec3<f32>(1.0, 0.9, 0.4);
        final_fragment_color = mix(final_fragment_color, selection_highlight_color, 0.4);
    }
    
    return vec4<f32>(final_fragment_color, 1.0);
}
