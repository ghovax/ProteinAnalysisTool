use std::sync::Arc;
use wgpu::util::DeviceExt;

use super::Camera;
use crate::protein::ProteinStore;

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct SphereInstance {
    pub position: [f32; 3],
    pub radius: f32,
    pub color: [f32; 3],
    pub _padding: f32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Uniforms {
    view_proj: [[f32; 4]; 4],
    camera_pos: [f32; 3],
    _padding: f32,
}

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

pub struct Renderer {
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    pipeline: wgpu::RenderPipeline,
    uniform_buffer: wgpu::Buffer,
    uniform_bind_group: wgpu::BindGroup,
    instance_buffer: wgpu::Buffer,
    instance_count: u32,
    depth_texture: wgpu::TextureView,
}

impl Renderer {
    pub fn new(
        device: Arc<wgpu::Device>,
        queue: Arc<wgpu::Queue>,
        surface_format: wgpu::TextureFormat,
        width: u32,
        height: u32,
    ) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Sphere Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("../shader.wgsl").into()),
        });

        let uniform_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Uniform Buffer"),
            size: std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let uniform_bind_group_layout =
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

        let uniform_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Uniform Bind Group"),
            layout: &uniform_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Render Pipeline Layout"),
            bind_group_layouts: &[&uniform_bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Sphere Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
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
                module: &shader,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleStrip,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
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

        let instance_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Instance Buffer"),
            size: 1024 * 1024, // 1MB initial size
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let depth_texture = Self::create_depth_texture(&device, width, height);

        Self {
            device,
            queue,
            pipeline,
            uniform_buffer,
            uniform_bind_group,
            instance_buffer,
            instance_count: 0,
            depth_texture,
        }
    }

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

    pub fn resize(&mut self, width: u32, height: u32) {
        self.depth_texture = Self::create_depth_texture(&self.device, width, height);
    }

    pub fn update_instances(&mut self, store: &ProteinStore) {
        let mut instances = Vec::new();

        for protein_arc in store.iter() {
            let protein = protein_arc.read().unwrap();
            if !protein.visible {
                continue;
            }

            let ca_positions = protein.ca_positions();
            let chain_ids: Vec<_> = protein.chain_ids();

            for (pos, chain_id) in ca_positions {
                let chain_idx = chain_ids.iter().position(|c| c == &chain_id).unwrap_or(0);
                let color = CHAIN_COLORS[chain_idx % CHAIN_COLORS.len()];

                instances.push(SphereInstance {
                    position: [pos.x, pos.y, pos.z],
                    radius: 1.5, // CA sphere radius in Angstroms
                    color,
                    _padding: 0.0,
                });
            }
        }

        self.instance_count = instances.len() as u32;

        if !instances.is_empty() {
            let data = bytemuck::cast_slice(&instances);
            if data.len() as u64 > self.instance_buffer.size() {
                self.instance_buffer = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Instance Buffer"),
                    contents: data,
                    usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                });
            } else {
                self.queue.write_buffer(&self.instance_buffer, 0, data);
            }
        }
    }

    pub fn render(
        &self,
        view: &wgpu::TextureView,
        camera: &Camera,
    ) -> wgpu::CommandBuffer {
        let uniforms = Uniforms {
            view_proj: camera.view_projection_matrix().to_cols_array_2d(),
            camera_pos: camera.position().to_array(),
            _padding: 0.0,
        };
        self.queue.write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(&uniforms));

        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Render Encoder"),
        });

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view,
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
                    view: &self.depth_texture,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            render_pass.set_pipeline(&self.pipeline);
            render_pass.set_bind_group(0, &self.uniform_bind_group, &[]);
            render_pass.set_vertex_buffer(0, self.instance_buffer.slice(..));
            render_pass.draw(0..4, 0..self.instance_count);
        }

        encoder.finish()
    }
}
