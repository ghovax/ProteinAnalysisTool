//! WGPU renderer implementation for protein visualization
//!
//! This module handles the creation of render pipelines, management of GPU buffers,
//! and the actual drawing of protein structures (spheres and lines)

use glam::Vec4Swizzles;
use pdbtbx::{ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain};
use std::collections::HashMap;
use std::sync::Arc;
use wgpu::util::DeviceExt;

use super::Camera;
use crate::protein::{ColorScheme, ProteinStore, Representation};

/// Instance data for rendering a sphere (atom)
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct SphereInstance {
    /// Position of the sphere in world space
    pub position: [f32; 3],
    /// Radius of the sphere
    pub radius: f32,
    /// Color of the sphere (RGB)
    pub color: [f32; 3],
    /// Selection factor (0.0 = normal, 1.0 = selected)
    pub selection_factor: f32,
}

/// Instance data for rendering a cylinder (bond or stick)
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct CylinderInstance {
    /// World position of the start of the cylinder
    pub start_position: [f32; 3],
    /// World position of the end of the cylinder
    pub end_position: [f32; 3],
    /// Radius of the cylinder
    pub radius: f32,
    /// Color of the cylinder (RGB)
    pub color: [f32; 3],
    /// Selection factor (0.0 = normal, 1.0 = selected)
    pub selection_factor: f32,
}

/// Vertex data for rendering a line (backbone)
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct LineVertex {
    /// Position of the vertex in world space
    pub position: [f32; 3],
    /// Color of the vertex (RGB)
    pub color: [f32; 3],
}

/// Global uniforms shared across all shaders
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Uniforms {
    /// Combined view and projection matrix
    view_projection_matrix: [[f32; 4]; 4],
    /// World position of the camera
    camera_world_position: [f32; 3],
    /// Explicit padding for alignment
    _padding: f32,
}

/// Default colors for protein chains
const CHAIN_COLORS: &[[f32; 3]] = &[
    [0.2, 0.6, 1.0], // Blue
    [1.0, 0.4, 0.4], // Red
    [0.4, 0.9, 0.4], // Green
    [1.0, 0.8, 0.2], // Yellow
    [0.9, 0.5, 0.9], // Magenta
    [0.5, 0.9, 0.9], // Cyan
    [1.0, 0.6, 0.3], // Orange
    [0.7, 0.7, 0.9], // Light purple
];

/// Colors assigned to chemical elements
const ELEMENT_COLORS: &[(&str, [f32; 3])] = &[
    ("C", [0.5, 0.5, 0.5]), // Carbon - gray
    ("N", [0.2, 0.2, 1.0]), // Nitrogen - blue
    ("O", [1.0, 0.2, 0.2]), // Oxygen - red
    ("S", [1.0, 1.0, 0.2]), // Sulfur - yellow
    ("H", [1.0, 1.0, 1.0]), // Hydrogen - white
    ("P", [1.0, 0.5, 0.0]), // Phosphorus - orange
];

/// Returns the color associated with a given element symbol
fn _get_element_color(element_symbol: &str) -> [f32; 3] {
    for (current_element_symbol, element_color) in ELEMENT_COLORS {
        if *current_element_symbol == element_symbol {
            return *element_color;
        }
    }
    [0.8, 0.8, 0.8] // Default gray
}

/// Maps a B-factor value to a color using a Blue-White-Red gradient
fn bfactor_to_color(
    temperature_factor: f32,
    minimum_bfactor: f32,
    maximum_bfactor: f32,
) -> [f32; 3] {
    let normalized_factor = if maximum_bfactor > minimum_bfactor {
        ((temperature_factor - minimum_bfactor) / (maximum_bfactor - minimum_bfactor))
            .clamp(0.0, 1.0)
    } else {
        0.5
    };
    // Blue (low) -> White (mid) -> Red (high)
    if normalized_factor < 0.5 {
        let interpolation_step = normalized_factor * 2.0;
        [interpolation_step, interpolation_step, 1.0]
    } else {
        let interpolation_step = (normalized_factor - 0.5) * 2.0;
        [1.0, 1.0 - interpolation_step, 1.0 - interpolation_step]
    }
}

/// The main protein renderer
pub struct Renderer {
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    sphere_render_pipeline: wgpu::RenderPipeline,
    cylinder_render_pipeline: wgpu::RenderPipeline,
    line_render_pipeline: wgpu::RenderPipeline,
    global_uniform_buffer: wgpu::Buffer,
    global_uniform_bind_group: wgpu::BindGroup,
    sphere_instance_buffer: wgpu::Buffer,
    sphere_count: u32,
    cylinder_instance_buffer: wgpu::Buffer,
    cylinder_count: u32,
    line_vertex_buffer: wgpu::Buffer,
    line_vertex_count: u32,
    depth_stencil_texture_view: wgpu::TextureView,

    // Text rendering
    font_system: glyphon::FontSystem,
    swash_cache: glyphon::SwashCache,
    text_atlas: glyphon::TextAtlas,
    text_renderer: glyphon::TextRenderer,
    text_buffer: glyphon::Buffer,
    viewport: glyphon::Viewport,
    render_format: wgpu::TextureFormat,
}

/// A text label at a 3D position
pub struct TextLabel {
    pub position: glam::Vec3,
    pub text: String,
    pub color: [f32; 4],
}

impl Renderer {
    /// Creates a new `Renderer` with initialized GPU pipelines
    pub fn new(
        device: Arc<wgpu::Device>,
        queue: Arc<wgpu::Queue>,
        surface_format: wgpu::TextureFormat,
        width: u32,
        height: u32,
    ) -> Self {
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("../shader.wgsl").into()),
        });

        let global_uniform_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Uniform Buffer"),
            size: std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let global_uniform_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Uniform Bind Group Layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });

        let global_uniform_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Uniform Bind Group"),
            layout: &global_uniform_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: global_uniform_buffer.as_entire_binding(),
            }],
        });

        let render_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("Render Pipeline Layout"),
                bind_group_layouts: &[&global_uniform_bind_group_layout],
                push_constant_ranges: &[],
            });

        // Sphere pipeline (billboard quads)
        let sphere_render_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("Sphere Pipeline"),
                layout: Some(&render_pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &shader_module,
                    entry_point: "vs_sphere",
                    buffers: &[wgpu::VertexBufferLayout {
                        array_stride: std::mem::size_of::<SphereInstance>() as u64,
                        step_mode: wgpu::VertexStepMode::Instance,
                        attributes: &[
                            wgpu::VertexAttribute {
                                offset: 0,
                                shader_location: 0,
                                format: wgpu::VertexFormat::Float32x3,
                            },
                            wgpu::VertexAttribute {
                                offset: 12,
                                shader_location: 1,
                                format: wgpu::VertexFormat::Float32,
                            },
                            wgpu::VertexAttribute {
                                offset: 16,
                                shader_location: 2,
                                format: wgpu::VertexFormat::Float32x3,
                            },
                            wgpu::VertexAttribute {
                                offset: 28,
                                shader_location: 3,
                                format: wgpu::VertexFormat::Float32,
                            },
                        ],
                    }],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &shader_module,
                    entry_point: "fs_sphere",
                    targets: &[Some(wgpu::ColorTargetState {
                        format: surface_format,
                        blend: Some(wgpu::BlendState::REPLACE),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleStrip,
                    ..Default::default()
                },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: wgpu::TextureFormat::Depth32Float,
                    depth_write_enabled: true,
                    depth_compare: wgpu::CompareFunction::Less,
                    stencil: wgpu::StencilState::default(),
                    bias: wgpu::DepthBiasState::default(),
                }),
                multisample: wgpu::MultisampleState::default(),
                multiview: None,
                cache: None,
            });

        // Cylinder pipeline
        let cylinder_render_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("Cylinder Pipeline"),
                layout: Some(&render_pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &shader_module,
                    entry_point: "vs_cylinder",
                    buffers: &[wgpu::VertexBufferLayout {
                        array_stride: std::mem::size_of::<CylinderInstance>() as u64,
                        step_mode: wgpu::VertexStepMode::Instance,
                        attributes: &[
                            wgpu::VertexAttribute {
                                offset: 0,
                                shader_location: 0,
                                format: wgpu::VertexFormat::Float32x3,
                            },
                            wgpu::VertexAttribute {
                                offset: 12,
                                shader_location: 1,
                                format: wgpu::VertexFormat::Float32x3,
                            },
                            wgpu::VertexAttribute {
                                offset: 24,
                                shader_location: 2,
                                format: wgpu::VertexFormat::Float32,
                            },
                            wgpu::VertexAttribute {
                                offset: 28,
                                shader_location: 3,
                                format: wgpu::VertexFormat::Float32x3,
                            },
                            wgpu::VertexAttribute {
                                offset: 40,
                                shader_location: 4,
                                format: wgpu::VertexFormat::Float32,
                            },
                        ],
                    }],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &shader_module,
                    entry_point: "fs_cylinder",
                    targets: &[Some(wgpu::ColorTargetState {
                        format: surface_format,
                        blend: Some(wgpu::BlendState::REPLACE),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleStrip,
                    ..Default::default()
                },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: wgpu::TextureFormat::Depth32Float,
                    depth_write_enabled: true,
                    depth_compare: wgpu::CompareFunction::Less,
                    stencil: wgpu::StencilState::default(),
                    bias: wgpu::DepthBiasState::default(),
                }),
                multisample: wgpu::MultisampleState::default(),
                multiview: None,
                cache: None,
            });

        // Line pipeline
        let line_render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Line Pipeline"),
            layout: Some(&render_pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader_module,
                entry_point: "vs_line",
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<LineVertex>() as u64,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &[
                        wgpu::VertexAttribute {
                            offset: 0,
                            shader_location: 0,
                            format: wgpu::VertexFormat::Float32x3,
                        },
                        wgpu::VertexAttribute {
                            offset: 12,
                            shader_location: 1,
                            format: wgpu::VertexFormat::Float32x3,
                        },
                    ],
                }],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader_module,
                entry_point: "fs_line",
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::LineList,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        let sphere_instance_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Sphere Buffer"),
            size: 1024 * 1024,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let cylinder_instance_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Cylinder Buffer"),
            size: 1024 * 1024,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let line_vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Line Buffer"),
            size: 1024 * 1024,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let depth_stencil_texture_view = Self::create_depth_texture(&device, width, height);

        // glyphon setup
        let mut font_system = glyphon::FontSystem::new();
        let swash_cache = glyphon::SwashCache::new();
        let cache = glyphon::Cache::new(&device);
        let mut text_atlas = glyphon::TextAtlas::new(&device, &queue, &cache, surface_format);
        let text_renderer = glyphon::TextRenderer::new(
            &mut text_atlas,
            &device,
            wgpu::MultisampleState::default(),
            None,
        );
        let text_buffer = glyphon::Buffer::new(&mut font_system, glyphon::Metrics::new(32.0, 40.0));
        let viewport = glyphon::Viewport::new(&device, &cache);

        Self {
            device,
            queue,
            sphere_render_pipeline,
            cylinder_render_pipeline,
            line_render_pipeline,
            global_uniform_buffer,
            global_uniform_bind_group,
            sphere_instance_buffer,
            sphere_count: 0,
            cylinder_instance_buffer,
            cylinder_count: 0,
            line_vertex_buffer,
            line_vertex_count: 0,
            depth_stencil_texture_view,
            font_system,
            swash_cache,
            text_atlas,
            text_renderer,
            text_buffer,
            viewport,
            render_format: surface_format,
        }
    }

    /// Internal helper to create the depth texture for Z-buffering
    fn create_depth_texture(device: &wgpu::Device, width: u32, height: u32) -> wgpu::TextureView {
        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Depth Texture"),
            size: wgpu::Extent3d {
                width: width.max(1),
                height: height.max(1),
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        texture.create_view(&wgpu::TextureViewDescriptor::default())
    }

    /// Resizes the depth texture when the window is resized
    pub fn resize(&mut self, width: u32, height: u32) {
        self.depth_stencil_texture_view = Self::create_depth_texture(&self.device, width, height);
        self.viewport
            .update(&self.queue, glyphon::Resolution { width, height });
    }

    /// Updates the GPU buffers with the latest protein instance data
    pub fn update_instances(
        &mut self,
        protein_data_store: &ProteinStore,
        currently_selected_atoms: &[(String, usize)],
        active_measurement_pairs: &[(usize, usize)],
    ) {
        let mut sphere_instances_collection = Vec::new();
        let mut cylinder_instances_collection = Vec::new();
        let mut line_vertices_collection = Vec::new();

        let mut atom_world_positions_lookup_table = HashMap::new();

        for protein_shared_reference in protein_data_store.iter() {
            let protein_locked_data = protein_shared_reference.read().unwrap();
            if !protein_locked_data.visible {
                continue;
            }

            let protein_identifier_name = &protein_locked_data.name;
            let available_chain_identifiers_collection: Vec<_> = protein_locked_data.chain_ids();
            let active_representation_mode = protein_locked_data.representation;
            let active_color_scheme_mode = protein_locked_data.color_scheme;

            // Get B-factor range for normalization during coloring
            let (minimum_bfactor_value, maximum_bfactor_value) =
                protein_locked_data.calculate_bfactor_range();

            // Helper to get color for a specific chain identifier
            let get_color_for_chain_identifier = |chain_id_string: &str| -> [f32; 3] {
                match active_color_scheme_mode {
                    ColorScheme::ByChain => {
                        let chain_index_position = available_chain_identifiers_collection
                            .iter()
                            .position(|c| c == chain_id_string)
                            .unwrap_or(0);
                        CHAIN_COLORS[chain_index_position % CHAIN_COLORS.len()]
                    }
                    ColorScheme::ByElement => [0.5, 0.5, 0.5], // Default for Alpha Carbon (Carbon)
                    ColorScheme::ByBFactor => [0.5, 0.5, 0.5], // Will be overridden per-atom
                    ColorScheme::BySecondary => [0.7, 0.7, 0.7], // Will be overridden per-residue
                    ColorScheme::Uniform(uniform_rgb_color) => uniform_rgb_color,
                }
            };

            let get_secondary_structure_color = |secondary_structure_type: crate::protein::structure::SecondaryStructureType| -> [f32; 3] {
                match secondary_structure_type {
                    crate::protein::structure::SecondaryStructureType::Helix => [1.0, 0.4, 1.0], // Magenta/Pink for helix
                    crate::protein::structure::SecondaryStructureType::Sheet => [1.0, 1.0, 0.2], // Yellow for sheet
                    crate::protein::structure::SecondaryStructureType::Other => [0.6, 0.6, 0.6], // Gray for loop/other
                }
            };

            // Pre-calculate atom data
            let mut protein_atom_data_list = Vec::new();
            for (atom_index, atom_hierarchy) in protein_locked_data.pdb.atoms_with_hierarchy().enumerate() {
                let atom_reference = atom_hierarchy.atom();
                let position_tuple = atom_reference.pos();
                let world_position = glam::Vec3::new(position_tuple.0 as f32, position_tuple.1 as f32, position_tuple.2 as f32);
                
                atom_world_positions_lookup_table.insert(
                    (protein_identifier_name.clone(), atom_index),
                    world_position,
                );

                let chain_id = atom_hierarchy.chain().id();
                let residue_number = atom_hierarchy.residue().serial_number();
                
                let secondary_structure_type = if protein_locked_data.pdb.is_residue_in_helix(chain_id, residue_number) {
                    crate::protein::structure::SecondaryStructureType::Helix
                } else if protein_locked_data.pdb.is_residue_in_sheet(chain_id, residue_number) {
                    crate::protein::structure::SecondaryStructureType::Sheet
                } else {
                    crate::protein::structure::SecondaryStructureType::Other
                };

                let color = match active_color_scheme_mode {
                    ColorScheme::ByChain => get_color_for_chain_identifier(chain_id),
                    ColorScheme::ByElement => _get_element_color(atom_reference.element().map(|e| e.symbol()).unwrap_or("?")),
                    ColorScheme::ByBFactor => bfactor_to_color(atom_reference.b_factor() as f32, minimum_bfactor_value, maximum_bfactor_value),
                    ColorScheme::BySecondary => get_secondary_structure_color(secondary_structure_type),
                    ColorScheme::Uniform(rgb) => rgb,
                };

                let is_selected = currently_selected_atoms.contains(&(protein_identifier_name.clone(), atom_index));

                protein_atom_data_list.push((world_position, color, is_selected, atom_reference.name().to_string(), atom_reference.element().map(|e| e.symbol().to_string())));
            }

            // Generate spheres
            match active_representation_mode {
                Representation::Spheres | Representation::BackboneAndSpheres => {
                    // Only Alpha Carbons (CA)
                    for (atom_index, (atom_world_position, atom_color, is_atom_selected, atom_name, _)) in protein_atom_data_list.iter().enumerate() {
                        if atom_name == "CA" {
                            sphere_instances_collection.push(SphereInstance {
                                position: atom_world_position.to_array(),
                                radius: 1.5,
                                color: *atom_color,
                                selection_factor: if *is_atom_selected { 1.0 } else { 0.0 },
                            });
                        }
                    }
                }
                Representation::BallAndStick => {
                    // All atoms as small spheres
                    for (atom_world_position, atom_color, is_atom_selected, _, _) in &protein_atom_data_list {
                        sphere_instances_collection.push(SphereInstance {
                            position: atom_world_position.to_array(),
                            radius: 0.4,
                            color: *atom_color,
                            selection_factor: if *is_atom_selected { 1.0 } else { 0.0 },
                        });
                    }
                }
                Representation::SpaceFilling => {
                    // All atoms as VdW spheres
                    for (atom_world_position, atom_color, is_atom_selected, _, _) in &protein_atom_data_list {
                        sphere_instances_collection.push(SphereInstance {
                            position: atom_world_position.to_array(),
                            radius: 1.7, // Average VdW radius
                            color: *atom_color,
                            selection_factor: if *is_atom_selected { 1.0 } else { 0.0 },
                        });
                    }
                }
                _ => {}
            }

            // Generate cylinders/lines for bonds
            match active_representation_mode {
                Representation::Sticks | Representation::BallAndStick => {
                    for current_atom_bond in &protein_locked_data.identified_atom_bonds {
                        let (start_atom_position, start_atom_color, is_start_atom_selected, _, _) = protein_atom_data_list[current_atom_bond.first_atom_index];
                        let (end_atom_position, end_atom_color, is_end_atom_selected, _, _) = protein_atom_data_list[current_atom_bond.second_atom_index];
                        
                        let cylinder_radius_value = if active_representation_mode == Representation::Sticks { 0.3 } else { 0.15 };
                        
                        // Create two cylinders meeting in the middle for per-atom coloring
                        let middle_position_vector = (start_atom_position + end_atom_position) * 0.5;
                        
                        cylinder_instances_collection.push(CylinderInstance {
                            start_position: start_atom_position.to_array(),
                            end_position: middle_position_vector.to_array(),
                            radius: cylinder_radius_value,
                            color: start_atom_color,
                            selection_factor: if is_start_atom_selected { 1.0 } else { 0.0 },
                        });
                        cylinder_instances_collection.push(CylinderInstance {
                            start_position: middle_position_vector.to_array(),
                            end_position: end_atom_position.to_array(),
                            radius: cylinder_radius_value,
                            color: end_atom_color,
                            selection_factor: if is_end_atom_selected { 1.0 } else { 0.0 },
                        });
                    }
                }
                Representation::Lines => {
                    for current_atom_bond in &protein_locked_data.identified_atom_bonds {
                        let (start_atom_position, start_atom_color, _, _, _) = protein_atom_data_list[current_atom_bond.first_atom_index];
                        let (end_atom_position, end_atom_color, _, _, _) = protein_atom_data_list[current_atom_bond.second_atom_index];
                        
                        line_vertices_collection.push(LineVertex { position: start_atom_position.to_array(), color: start_atom_color });
                        line_vertices_collection.push(LineVertex { position: end_atom_position.to_array(), color: end_atom_color });
                    }
                }
                _ => {}
            }

            // Generate backbone backbone trace lines if needed
            if active_representation_mode == Representation::Backbone
                || active_representation_mode == Representation::BackboneAndSpheres
            {
                if active_color_scheme_mode == ColorScheme::BySecondary {
                    let backbone_segment_collection =
                        protein_locked_data.get_backbone_segments_with_secondary_structure();
                    for (
                        segment_start_position,
                        segment_end_position,
                        _,
                        secondary_structure_type,
                    ) in backbone_segment_collection
                    {
                        let final_line_segment_color =
                            get_secondary_structure_color(secondary_structure_type);

                        line_vertices_collection.push(LineVertex {
                            position: [
                                segment_start_position.x,
                                segment_start_position.y,
                                segment_start_position.z,
                            ],
                            color: final_line_segment_color,
                        });
                        line_vertices_collection.push(LineVertex {
                            position: [
                                segment_end_position.x,
                                segment_end_position.y,
                                segment_end_position.z,
                            ],
                            color: final_line_segment_color,
                        });
                    }
                } else {
                    let backbone_segment_collection =
                        protein_locked_data.get_backbone_segments_for_rendering();
                    for (segment_start_position, segment_end_position, target_chain_identifier) in
                        backbone_segment_collection
                    {
                        let final_line_segment_color =
                            get_color_for_chain_identifier(&target_chain_identifier);

                        line_vertices_collection.push(LineVertex {
                            position: [
                                segment_start_position.x,
                                segment_start_position.y,
                                segment_start_position.z,
                            ],
                            color: final_line_segment_color,
                        });
                        line_vertices_collection.push(LineVertex {
                            position: [
                                segment_end_position.x,
                                segment_end_position.y,
                                segment_end_position.z,
                            ],
                            color: final_line_segment_color,
                        });
                    }
                }
            }
        }

        // Add measurement lines between selected atoms
        let measurement_line_color = [1.0, 1.0, 1.0]; // Pure white for measurement lines
        for &(first_atom_selection_index, second_atom_selection_index) in active_measurement_pairs {
            if let (Some(first_atom_selection_handle), Some(second_atom_selection_handle)) = (
                currently_selected_atoms.get(first_atom_selection_index),
                currently_selected_atoms.get(second_atom_selection_index),
            ) {
                if let (Some(first_atom_world_position), Some(second_atom_world_position)) = (
                    atom_world_positions_lookup_table.get(first_atom_selection_handle),
                    atom_world_positions_lookup_table.get(second_atom_selection_handle),
                ) {
                    line_vertices_collection.push(LineVertex {
                        position: [
                            first_atom_world_position.x,
                            first_atom_world_position.y,
                            first_atom_world_position.z,
                        ],
                        color: measurement_line_color,
                    });
                    line_vertices_collection.push(LineVertex {
                        position: [
                            second_atom_world_position.x,
                            second_atom_world_position.y,
                            second_atom_world_position.z,
                        ],
                        color: measurement_line_color,
                    });
                }
            }
        }

        // Update GPU cylinder instance buffer
        self.cylinder_count = cylinder_instances_collection.len() as u32;
        if !cylinder_instances_collection.is_empty() {
            let serialized_cylinder_instance_data =
                bytemuck::cast_slice(&cylinder_instances_collection);
            if serialized_cylinder_instance_data.len() as u64 > self.cylinder_instance_buffer.size()
            {
                self.cylinder_instance_buffer =
                    self.device
                        .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                            label: Some("Cylinder Instance Buffer"),
                            contents: serialized_cylinder_instance_data,
                            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                        });
            } else {
                self.queue.write_buffer(
                    &self.cylinder_instance_buffer,
                    0,
                    serialized_cylinder_instance_data,
                );
            }
        }

        // Update GPU sphere instance buffer
        self.sphere_count = sphere_instances_collection.len() as u32;
        if !sphere_instances_collection.is_empty() {
            let serialized_sphere_instance_data =
                bytemuck::cast_slice(&sphere_instances_collection);
            if serialized_sphere_instance_data.len() as u64 > self.sphere_instance_buffer.size() {
                self.sphere_instance_buffer =
                    self.device
                        .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                            label: Some("Sphere Instance Buffer"),
                            contents: serialized_sphere_instance_data,
                            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                        });
            } else {
                self.queue.write_buffer(
                    &self.sphere_instance_buffer,
                    0,
                    serialized_sphere_instance_data,
                );
            }
        }

        // Update GPU line vertex buffer
        self.line_vertex_count = line_vertices_collection.len() as u32;
        if !line_vertices_collection.is_empty() {
            let serialized_line_vertex_data = bytemuck::cast_slice(&line_vertices_collection);
            if serialized_line_vertex_data.len() as u64 > self.line_vertex_buffer.size() {
                self.line_vertex_buffer =
                    self.device
                        .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                            label: Some("Line Vertex Buffer"),
                            contents: serialized_line_vertex_data,
                            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                        });
            } else {
                self.queue
                    .write_buffer(&self.line_vertex_buffer, 0, serialized_line_vertex_data);
            }
        }
    }

    /// Encodes and returns a command buffer for rendering the current frame
    pub fn render(
        &mut self,
        target_texture_view: &wgpu::TextureView,
        camera_object_handle: &Camera,
        text_labels_collection: &[TextLabel],
        viewport_width: u32,
        viewport_height: u32,
    ) -> wgpu::CommandBuffer {
        let global_uniform_values_structure = Uniforms {
            view_projection_matrix: camera_object_handle
                .view_projection_matrix()
                .to_cols_array_2d(),
            camera_world_position: camera_object_handle.position().to_array(),
            _padding: 0.0,
        };
        self.queue.write_buffer(
            &self.global_uniform_buffer,
            0,
            bytemuck::bytes_of(&global_uniform_values_structure),
        );

        let mut graphics_command_encoder =
            self.device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Main Render Command Encoder"),
                });

        {
            let mut active_render_pass =
                graphics_command_encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("Protein Visualization Render Pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: target_texture_view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color {
                                r: 0.1,
                                g: 0.1,
                                b: 0.15,
                                a: 1.0,
                            }),
                            store: wgpu::StoreOp::Store,
                        },
                    })],
                    depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                        view: &self.depth_stencil_texture_view,
                        depth_ops: Some(wgpu::Operations {
                            load: wgpu::LoadOp::Clear(1.0),
                            store: wgpu::StoreOp::Store,
                        }),
                        stencil_ops: None,
                    }),
                    timestamp_writes: None,
                    occlusion_query_set: None,
                });

            // Draw backbone trace lines first
            if self.line_vertex_count > 0 {
                active_render_pass.set_pipeline(&self.line_render_pipeline);
                active_render_pass.set_bind_group(0, &self.global_uniform_bind_group, &[]);
                active_render_pass.set_vertex_buffer(0, self.line_vertex_buffer.slice(..));
                active_render_pass.draw(0..self.line_vertex_count, 0..1);
            }

            // Draw atom spheres
            if self.sphere_count > 0 {
                active_render_pass.set_pipeline(&self.sphere_render_pipeline);
                active_render_pass.set_bind_group(0, &self.global_uniform_bind_group, &[]);
                active_render_pass.set_vertex_buffer(0, self.sphere_instance_buffer.slice(..));
                active_render_pass.draw(0..4, 0..self.sphere_count);
            }

            // Draw cylinders (sticks/bonds)
            if self.cylinder_count > 0 {
                active_render_pass.set_pipeline(&self.cylinder_render_pipeline);
                active_render_pass.set_bind_group(0, &self.global_uniform_bind_group, &[]);
                active_render_pass.set_vertex_buffer(0, self.cylinder_instance_buffer.slice(..));
                active_render_pass.draw(0..4, 0..self.cylinder_count);
            }
        }

        // Prepare and render UI text labels
        let view_projection_matrix = camera_object_handle.view_projection_matrix();

        let mut combined_ui_text_string = String::new();
        for current_text_label in text_labels_collection {
            if current_text_label.position != glam::Vec3::ZERO {
                let clip_space_position =
                    view_projection_matrix * current_text_label.position.extend(1.0);
                if clip_space_position.w > 0.0 {
                    let normalized_device_coordinates =
                        clip_space_position.xyz() / clip_space_position.w;
                    if normalized_device_coordinates.x.abs() <= 1.0
                        && normalized_device_coordinates.y.abs() <= 1.0
                        && normalized_device_coordinates.z >= 0.0
                        && normalized_device_coordinates.z <= 1.0
                    {
                        // Position on screen (placeholder logic for combined string)
                        combined_ui_text_string.push_str(&current_text_label.text);
                        combined_ui_text_string.push(' ');
                    }
                }
            } else {
                // Fixed UI labels
                combined_ui_text_string.push_str(&current_text_label.text);
                combined_ui_text_string.push('\n');
            }
        }

        if !combined_ui_text_string.is_empty() {
            self.text_buffer.set_text(
                &mut self.font_system,
                &combined_ui_text_string,
                glyphon::Attrs::new().family(glyphon::Family::SansSerif),
                glyphon::Shaping::Advanced,
            );
            self.text_buffer
                .shape_until_scroll(&mut self.font_system, true);

            self.text_renderer
                .prepare(
                    &self.device,
                    &self.queue,
                    &mut self.font_system,
                    &mut self.text_atlas,
                    &self.viewport,
                    [glyphon::TextArea {
                        buffer: &self.text_buffer,
                        left: 20.0,
                        top: 20.0,
                        scale: 1.0,
                        bounds: glyphon::TextBounds {
                            left: 0,
                            top: 0,
                            right: viewport_width as i32,
                            bottom: viewport_height as i32,
                        },
                        default_color: glyphon::Color::rgb(255, 255, 255),
                        custom_glyphs: &[],
                    }],
                    &mut self.swash_cache,
                )
                .unwrap();

            {
                let mut text_rendering_pass =
                    graphics_command_encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                        label: Some("Text Overlay Rendering Pass"),
                        color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                            view: target_texture_view,
                            resolve_target: None,
                            ops: wgpu::Operations {
                                load: wgpu::LoadOp::Load,
                                store: wgpu::StoreOp::Store,
                            },
                        })],
                        depth_stencil_attachment: None,
                        timestamp_writes: None,
                        occlusion_query_set: None,
                    });

                self.text_renderer
                    .render(&self.text_atlas, &self.viewport, &mut text_rendering_pass)
                    .unwrap();
            }
        }

        graphics_command_encoder.finish()
    }
}
