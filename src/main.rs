use glam::Vec3;
use mlua::{FromLua, Lua, UserData, UserDataMethods, Value};
use notify::Watcher;
use std::{path::PathBuf, sync::{Arc, Mutex}};
use winit::{event::*, event_loop::{ControlFlow, EventLoop}, window::WindowBuilder};

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Vertex {
    position: [f32; 3],
    color: [f32; 3],
}

#[derive(Clone, Copy, Debug)]
struct Transform {
    position: Vec3,
}

impl UserData for Transform {
    fn add_methods<'lua, M: UserDataMethods<'lua, Self>>(methods: &mut M) {
        methods.add_method_mut("set_pos", |_, this, (x, y, z): (f32, f32, f32)| {
            this.position = Vec3::new(x, y, z);
            Ok(())
        });
    }
}

impl<'lua> FromLua<'lua> for Transform {
    fn from_lua(value: Value<'lua>, _: &'lua Lua) -> mlua::Result<Self> {
        match value {
            Value::UserData(ud) => Ok(*ud.borrow::<Self>()?),
            _ => Err(mlua::Error::FromLuaConversionError { from: value.type_name(), to: "Transform", message: None }),
        }
    }
}

struct ScriptEngine {
    lua: Lua,
    transform: Arc<Mutex<Transform>>,
    rx: crossbeam_channel::Receiver<PathBuf>,
    _watcher: notify::RecommendedWatcher,
}

impl ScriptEngine {
    fn init(initial_transform: Arc<Mutex<Transform>>) -> Self {
        let lua = Lua::new();
        let (tx, rx) = crossbeam_channel::unbounded::<PathBuf>();
        let mut watcher = notify::recommended_watcher(move |res: notify::Result<notify::Event>| {
            if let Ok(event) = res {
                if event.kind.is_modify() {
                    for path in event.paths { let _ = tx.send(path); }
                }
            }
        }).unwrap();
        watcher.watch(std::path::Path::new("."), notify::RecursiveMode::NonRecursive).unwrap();
        
        let initial_val = *initial_transform.lock().unwrap();
        lua.globals().set("transform", initial_val).unwrap();

        Self { lua, transform: initial_transform, rx, _watcher: watcher }
    }

    fn update(&mut self) {
    match std::fs::read_to_string("player.lua") {
        Ok(code) => {
            if let Err(e) = self.lua.load(&code).exec() {
                eprintln!("LUA RUNTIME ERROR: {}", e);
            }
        }
        Err(e) => {
            eprintln!("FILE ERROR: Could not find 'player.lua' in the root directory. Error: {}", e);
        }
    }

    if let Ok(lua_t) = self.lua.globals().get::<_, Transform>("transform") {
        let mut t = self.transform.lock().unwrap();
        *t = lua_t;
        // THIS MUST PRINT IN YOUR TERMINAL
        println!("RUST: Position is now {}, {}", t.position.x, t.position.y);
    }
}
}

const VERTICES: &[Vertex] = &[
    Vertex { position: [-0.2, -0.2, 0.0], color: [1.0, 0.0, 0.0] },
    Vertex { position: [0.2, -0.2, 0.0], color: [0.0, 1.0, 0.0] },
    Vertex { position: [0.0, 0.2, 0.0], color: [0.0, 0.0, 1.0] },
];

async fn run() {
    let event_loop = EventLoop::new().unwrap();
    let window = Arc::new(WindowBuilder::new().with_title("Unity-like Rust Engine").build(&event_loop).unwrap());
    let instance = wgpu::Instance::default();
    let surface = instance.create_surface(Arc::clone(&window)).unwrap();
    let adapter = instance.request_adapter(&wgpu::RequestAdapterOptions::default()).await.unwrap();
    let (device, queue) = adapter.request_device(&wgpu::DeviceDescriptor::default(), None).await.unwrap();

    let size = window.inner_size();
    let config = surface.get_default_config(&adapter, size.width, size.height).unwrap();
    surface.configure(&device, &config);

    // Buffer for u_offset (Must be 16 bytes for vec4 alignment)
    let uniform_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Uniform Buffer"),
        size: 16,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        entries: &[wgpu::BindGroupLayoutEntry {
            binding: 0,
            visibility: wgpu::ShaderStages::VERTEX,
            ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None },
            count: None,
        }],
        label: None,
    });

    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        layout: &bind_group_layout,
        entries: &[wgpu::BindGroupEntry { binding: 0, resource: uniform_buffer.as_entire_binding() }],
        label: None,
    });

    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: None,
        source: wgpu::ShaderSource::Wgsl(std::borrow::Cow::Borrowed(include_str!("shader.wgsl"))),
    });

    let render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: None,
        layout: Some(&device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { bind_group_layouts: &[&bind_group_layout], ..Default::default() })),
        vertex: wgpu::VertexState { 
            module: &shader, entry_point: "vs_main", 
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<Vertex>() as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3],
            }],
            compilation_options: Default::default(),
        },
        fragment: Some(wgpu::FragmentState { module: &shader, entry_point: "fs_main", targets: &[Some(config.format.into())], compilation_options: Default::default() }),
        primitive: wgpu::PrimitiveState::default(),
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        multiview: None,
        cache: None,
    });

    let vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: None,
        size: (VERTICES.len() * std::mem::size_of::<Vertex>()) as u64,
        usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&vertex_buffer, 0, bytemuck::cast_slice(VERTICES));

    let transform = Arc::new(Mutex::new(Transform { position: Vec3::ZERO }));
    let mut script_engine = ScriptEngine::init(transform.clone());

    event_loop.set_control_flow(ControlFlow::Poll);
    let _ = event_loop.run(move |event, elwt| match event {
        Event::WindowEvent { event: WindowEvent::CloseRequested, .. } => elwt.exit(),
        Event::AboutToWait => {
            script_engine.update();
            let pos = transform.lock().unwrap().position;
            // Write 16 bytes: x, y, z, padding
            queue.write_buffer(&uniform_buffer, 0, bytemuck::cast_slice(&[pos.x, pos.y, pos.z, 0.0f32]));
            window.request_redraw();
        }
        Event::WindowEvent { event: WindowEvent::RedrawRequested, .. } => {
            let frame = surface.get_current_texture().unwrap();
            let view = frame.texture.create_view(&wgpu::TextureViewDescriptor::default());
            let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
            {
                let mut rp = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment { view: &view, resolve_target: None, ops: wgpu::Operations { load: wgpu::LoadOp::Clear(wgpu::Color::BLACK), store: wgpu::StoreOp::Store } })],
                    ..Default::default()
                });
                rp.set_pipeline(&render_pipeline);
                rp.set_bind_group(0, &bind_group, &[]);
                rp.set_vertex_buffer(0, vertex_buffer.slice(..));
                rp.draw(0..VERTICES.len() as u32, 0..1);
            }
            queue.submit(std::iter::once(encoder.finish()));
            frame.present();
        }
        _ => {}
    });
}

fn main() { pollster::block_on(run()); }