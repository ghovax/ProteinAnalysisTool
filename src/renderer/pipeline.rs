//! WGPU renderer implementation for protein visualization
//!
//! This module handles the creation of render pipelines, management of GPU buffers,
//! and the actual drawing of protein structures (spheres and lines)

use std::sync::Arc;
use wgpu::util::DeviceExt;

use super::Camera;
use crate::protein::{ProteinStore, Representation, ColorScheme};

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
    /// Explicit padding for alignment
    pub _padding: f32,
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
    [0.2, 0.6, 1.0],  // Blue
    [1.0, 0.4, 0.4],  // Red
    [0.4, 0.9, 0.4],  // Green
    [1.0, 0.8, 0.2],  // Yellow
    [0.9, 0.5, 0.9],  // Magenta
    [0.5, 0.9, 0.9],  // Cyan
    [1.0, 0.6, 0.3],  // Orange
    [0.7, 0.7, 0.9],  // Light purple
];

/// Colors assigned to chemical elements
const ELEMENT_COLORS: &[(&str, [f32; 3])] = &[
    ("C", [0.5, 0.5, 0.5]),   // Carbon - gray
    ("N", [0.2, 0.2, 1.0]),   // Nitrogen - blue
    ("O", [1.0, 0.2, 0.2]),   // Oxygen - red
    ("S", [1.0, 1.0, 0.2]),   // Sulfur - yellow
    ("H", [1.0, 1.0, 1.0]),   // Hydrogen - white
    ("P", [1.0, 0.5, 0.0]),   // Phosphorus - orange
];

/// Returns the color associated with a given element symbol
fn get_element_color(element_symbol: &str) -> [f32; 3] {
    for (current_element_symbol, element_color) in ELEMENT_COLORS {
        if *current_element_symbol == element_symbol {
            return *element_color;
        }
    }
    [0.8, 0.8, 0.8] // Default gray
}

/// Maps a B-factor value to a color using a Blue-White-Red gradient
fn bfactor_to_color(temperature_factor: f32, minimum_bfactor: f32, maximum_bfactor: f32) -> [f32; 3] {
    let normalized_factor = if maximum_bfactor > minimum_bfactor {
        ((temperature_factor - minimum_bfactor) / (maximum_bfactor - minimum_bfactor)).clamp(0.0, 1.0)
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
    line_render_pipeline: wgpu::RenderPipeline,
    global_uniform_buffer: wgpu::Buffer,
    global_uniform_bind_group: wgpu::BindGroup,
    sphere_instance_buffer: wgpu::Buffer,
    sphere_count: u32,
    line_vertex_buffer: wgpu::Buffer,
    line_vertex_count: u32,
    depth_stencil_texture_view: wgpu::TextureView,
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

        let render_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Render Pipeline Layout"),
            bind_group_layouts: &[&global_uniform_bind_group_layout],
            push_constant_ranges: &[],
        });

        // Sphere pipeline (billboard quads)
        let sphere_render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
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

        let line_vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Line Buffer"),
            size: 1024 * 1024,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let depth_stencil_texture_view = Self::create_depth_texture(&device, width, height);

        Self {
            device,
            queue,
            sphere_render_pipeline,
            line_render_pipeline,
            global_uniform_buffer,
            global_uniform_bind_group,
            sphere_instance_buffer,
            sphere_count: 0,
            line_vertex_buffer,
            line_vertex_count: 0,
            depth_stencil_texture_view,
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
    }

    /// Updates the GPU buffers with the latest protein instance data
    pub fn update_instances(&mut self, store: &ProteinStore) {
        let mut sphere_instances = Vec::new();
        let mut line_vertices = Vec::new();

        for protein_shared_reference in store.iter() {
            let protein_locked_data = protein_shared_reference.read().unwrap();
            if !protein_locked_data.visible {
                continue;
            }

            let available_chain_identifiers: Vec<_> = protein_locked_data.chain_ids();
            let active_representation_mode = protein_locked_data.representation;
            let active_color_scheme = protein_locked_data.color_scheme;

            // Get B-factor range for coloring
            let (minimum_bfactor_value, maximum_bfactor_value) = protein_locked_data.bfactor_range();

            // Helper to get color for a chain
            let get_color_for_chain_identifier = |chain_id: &str| -> [f32; 3] {
                match active_color_scheme {
                    ColorScheme::ByChain => {
                        let chain_idx = available_chain_identifiers.iter().position(|c| c == chain_id).unwrap_or(0);
                        CHAIN_COLORS[chain_idx % CHAIN_COLORS.len()]
                    }
                    ColorScheme::ByElement => [0.5, 0.5, 0.5], // CA is always carbon
                    ColorScheme::ByBFactor => [0.5, 0.5, 0.5], // Will be overridden per-atom
                    ColorScheme::BySecondary => [0.7, 0.7, 0.7], // TODO: Implement
                    ColorScheme::Uniform(c) => c,
                }
            };

            // Generate spheres if needed
            if active_representation_mode == Representation::Spheres || active_representation_mode == Representation::BackboneAndSpheres {
                let alpha_carbon_data_collection = protein_locked_data.ca_with_bfactor();
                for (atom_position_vector, target_chain_identifier, atom_temperature_factor) in alpha_carbon_data_collection {
                    let final_atom_color = match active_color_scheme {
                        ColorScheme::ByBFactor => bfactor_to_color(atom_temperature_factor, minimum_bfactor_value, maximum_bfactor_value),
                        _ => get_color_for_chain_identifier(&target_chain_identifier),
                    };

                    sphere_instances.push(SphereInstance {
                        position: [atom_position_vector.x, atom_position_vector.y, atom_position_vector.z],
                        radius: 1.5,
                        color: final_atom_color,
                        _padding: 0.0,
                    });
                }
            }

            // Generate backbone lines if needed
            if active_representation_mode == Representation::Backbone || active_representation_mode == Representation::BackboneAndSpheres {
                let backbone_segment_collection = protein_locked_data.backbone_segments();
                for (segment_start_position, segment_end_position, target_chain_identifier) in backbone_segment_collection {
                    let final_segment_color = get_color_for_chain_identifier(&target_chain_identifier);

                    line_vertices.push(LineVertex {
                        position: [segment_start_position.x, segment_start_position.y, segment_start_position.z],
                        color: final_segment_color,
                    });
                    line_vertices.push(LineVertex {
                        position: [segment_end_position.x, segment_end_position.y, segment_end_position.z],
                        color: final_segment_color,
                    });
                }
            }
        }

        // Update sphere buffer
        self.sphere_count = sphere_instances.len() as u32;
        if !sphere_instances.is_empty() {
            let serialized_buffer_data = bytemuck::cast_slice(&sphere_instances);
            if serialized_buffer_data.len() as u64 > self.sphere_instance_buffer.size() {
                self.sphere_instance_buffer = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Sphere Buffer"),
                    contents: serialized_buffer_data,
                    usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                });
            } else {
                self.queue.write_buffer(&self.sphere_instance_buffer, 0, serialized_buffer_data);
            }
        }

        // Update line buffer
        self.line_vertex_count = line_vertices.len() as u32;
        if !line_vertices.is_empty() {
            let serialized_buffer_data = bytemuck::cast_slice(&line_vertices);
            if serialized_buffer_data.len() as u64 > self.line_vertex_buffer.size() {
                self.line_vertex_buffer = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Line Buffer"),
                    contents: serialized_buffer_data,
                    usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                });
            } else {
                self.queue.write_buffer(&self.line_vertex_buffer, 0, serialized_buffer_data);
            }
        }
    }

    /// Encodes and returns a command buffer for rendering the current frame
    pub fn render(
        &self,
        target_texture_view: &wgpu::TextureView,
        camera_object: &Camera,
    ) -> wgpu::CommandBuffer {
        let global_uniform_values = Uniforms {
            view_projection_matrix: camera_object.view_projection_matrix().to_cols_array_2d(),
            camera_world_position: camera_object.position().to_array(),
            _padding: 0.0,
        };
        self.queue.write_buffer(&self.global_uniform_buffer, 0, bytemuck::bytes_of(&global_uniform_values));

        let mut command_encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Render Encoder"),
        });

        {
            let mut active_render_pass = command_encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
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

            // Draw lines first (behind spheres)
            if self.line_vertex_count > 0 {
                active_render_pass.set_pipeline(&self.line_render_pipeline);
                active_render_pass.set_bind_group(0, &self.global_uniform_bind_group, &[]);
                active_render_pass.set_vertex_buffer(0, self.line_vertex_buffer.slice(..));
                active_render_pass.draw(0..self.line_vertex_count, 0..1);
            }

            // Draw spheres
            if self.sphere_count > 0 {
                active_render_pass.set_pipeline(&self.sphere_render_pipeline);
                active_render_pass.set_bind_group(0, &self.global_uniform_bind_group, &[]);
                active_render_pass.set_vertex_buffer(0, self.sphere_instance_buffer.slice(..));
                active_render_pass.draw(0..4, 0..self.sphere_count);
            }
        }

        command_encoder.finish()
    }
}