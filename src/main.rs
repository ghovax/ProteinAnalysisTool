mod lua_api;
mod protein;
mod renderer;

use notify::Watcher;
use std::path::PathBuf;
use std::sync::{Arc, RwLock};
use winit::{
    event::*,
    event_loop::{ControlFlow, EventLoop},
    keyboard::{KeyCode, PhysicalKey},
    window::WindowBuilder,
};

use lua_api::ScriptEngine;
use protein::ProteinStore;
use renderer::{Camera, Renderer};

struct App {
    surface: wgpu::Surface<'static>,
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    config: wgpu::SurfaceConfiguration,
    renderer: Renderer,
    camera: Camera,
    script_engine: ScriptEngine,
    store: Arc<RwLock<ProteinStore>>,

    // Input state
    mouse_pressed: bool,
    last_mouse_pos: Option<(f64, f64)>,

    // Script hot-reload
    script_rx: crossbeam_channel::Receiver<PathBuf>,
    _watcher: notify::RecommendedWatcher,
}

impl App {
    async fn new(window: Arc<winit::window::Window>) -> Self {
        let size = window.inner_size();
        let instance = wgpu::Instance::default();
        let surface = instance.create_surface(window.clone()).unwrap();
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .expect("Failed to find adapter");

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .expect("Failed to create device");

        let device = Arc::new(device);
        let queue = Arc::new(queue);

        let surface_caps = surface.get_capabilities(&adapter);
        let surface_format = surface_caps
            .formats
            .iter()
            .find(|f| f.is_srgb())
            .copied()
            .unwrap_or(surface_caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let renderer = Renderer::new(
            device.clone(),
            queue.clone(),
            surface_format,
            size.width,
            size.height,
        );

        let camera = Camera::new(size.width as f32 / size.height as f32);

        let store = Arc::new(RwLock::new(ProteinStore::new()));
        let script_engine = ScriptEngine::new(store.clone()).expect("Failed to create Lua engine");

        // Set up file watcher for hot-reload
        let (script_tx, script_rx) = crossbeam_channel::unbounded::<PathBuf>();
        let mut watcher =
            notify::recommended_watcher(move |res: notify::Result<notify::Event>| {
                if let Ok(event) = res {
                    if event.kind.is_modify() {
                        for path in event.paths {
                            if path.extension().map(|e| e == "lua").unwrap_or(false) {
                                let _ = script_tx.send(path);
                            }
                        }
                    }
                }
            })
            .expect("Failed to create watcher");

        // Watch scripts directory and current directory
        let _ = watcher.watch(
            std::path::Path::new("scripts"),
            notify::RecursiveMode::Recursive,
        );
        let _ = watcher.watch(
            std::path::Path::new("."),
            notify::RecursiveMode::NonRecursive,
        );

        // Run initial script
        if std::path::Path::new("scripts/init.lua").exists() {
            if let Err(e) = script_engine.run_file("scripts/init.lua") {
                eprintln!("Script error: {}", e);
            }
        } else if std::path::Path::new("init.lua").exists() {
            if let Err(e) = script_engine.run_file("init.lua") {
                eprintln!("Script error: {}", e);
            }
        }

        Self {
            surface,
            device,
            queue,
            config,
            renderer,
            camera,
            script_engine,
            store,
            mouse_pressed: false,
            last_mouse_pos: None,
            script_rx,
            _watcher: watcher,
        }
    }

    fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);
            self.renderer.resize(new_size.width, new_size.height);
            self.camera
                .set_aspect(new_size.width as f32 / new_size.height as f32);
        }
    }

    fn update(&mut self) {
        // Check for script hot-reload
        while let Ok(path) = self.script_rx.try_recv() {
            println!("Reloading: {:?}", path);
            if let Some(path_str) = path.to_str() {
                if let Err(e) = self.script_engine.run_file(path_str) {
                    eprintln!("Script error: {}", e);
                }
            }
        }

        // Update renderer with current protein data
        {
            let store = self.store.read().unwrap();
            self.renderer.update_instances(&store);
        }

        // Auto-focus camera on first protein if we haven't moved it
        if self.camera.target == glam::Vec3::ZERO {
            let store = self.store.read().unwrap();
            let first_protein = store.iter().next().cloned();
            drop(store);

            if let Some(protein_arc) = first_protein {
                let protein = protein_arc.read().unwrap();
                let center = protein.center_of_mass();
                let (min, max) = protein.bounding_box();
                let radius = (max - min).length() / 2.0;
                drop(protein);
                self.camera.focus_on(center, radius);
            }
        }
    }

    fn render(&self) -> Result<(), wgpu::SurfaceError> {
        let output = self.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        let command_buffer = self.renderer.render(&view, &self.camera);
        self.queue.submit(std::iter::once(command_buffer));
        output.present();

        Ok(())
    }

    fn handle_mouse_input(&mut self, state: ElementState, button: MouseButton) {
        if button == MouseButton::Left {
            self.mouse_pressed = state == ElementState::Pressed;
            if !self.mouse_pressed {
                self.last_mouse_pos = None;
            }
        }
    }

    fn handle_mouse_move(&mut self, position: (f64, f64)) {
        if self.mouse_pressed {
            if let Some((last_x, last_y)) = self.last_mouse_pos {
                let delta_x = (position.0 - last_x) as f32 * 0.01;
                let delta_y = (position.1 - last_y) as f32 * 0.01;
                self.camera.rotate(-delta_x, -delta_y);
            }
        }
        self.last_mouse_pos = Some(position);
    }

    fn handle_scroll(&mut self, delta: f32) {
        self.camera.zoom(delta);
    }

    fn handle_key(&mut self, key: KeyCode) {
        match key {
            KeyCode::KeyR => {
                // Reset camera
                self.camera = Camera::new(self.config.width as f32 / self.config.height as f32);
                let store = self.store.read().unwrap();
                let first_protein = store.iter().next().cloned();
                drop(store);

                if let Some(protein_arc) = first_protein {
                    let protein = protein_arc.read().unwrap();
                    let center = protein.center_of_mass();
                    let (min, max) = protein.bounding_box();
                    let radius = (max - min).length() / 2.0;
                    drop(protein);
                    self.camera.focus_on(center, radius);
                }
            }
            KeyCode::Escape => {
                std::process::exit(0);
            }
            _ => {}
        }
    }
}

async fn run() {
    let event_loop = EventLoop::new().unwrap();
    let window = Arc::new(
        WindowBuilder::new()
            .with_title("Protein Viewer")
            .with_inner_size(winit::dpi::LogicalSize::new(1280, 800))
            .build(&event_loop)
            .unwrap(),
    );

    let mut app = App::new(window.clone()).await;

    event_loop.set_control_flow(ControlFlow::Poll);
    let _ = event_loop.run(move |event, elwt| match event {
        Event::WindowEvent { event, .. } => match event {
            WindowEvent::CloseRequested => elwt.exit(),
            WindowEvent::Resized(size) => app.resize(size),
            WindowEvent::MouseInput { state, button, .. } => app.handle_mouse_input(state, button),
            WindowEvent::CursorMoved { position, .. } => {
                app.handle_mouse_move((position.x, position.y))
            }
            WindowEvent::MouseWheel { delta, .. } => {
                let scroll = match delta {
                    MouseScrollDelta::LineDelta(_, y) => y,
                    MouseScrollDelta::PixelDelta(pos) => pos.y as f32 / 100.0,
                };
                app.handle_scroll(scroll);
            }
            WindowEvent::KeyboardInput {
                event:
                    KeyEvent {
                        physical_key: PhysicalKey::Code(keycode),
                        state: ElementState::Pressed,
                        ..
                    },
                ..
            } => app.handle_key(keycode),
            WindowEvent::RedrawRequested => {
                app.update();
                match app.render() {
                    Ok(_) => {}
                    Err(wgpu::SurfaceError::Lost) => app.resize(winit::dpi::PhysicalSize {
                        width: app.config.width,
                        height: app.config.height,
                    }),
                    Err(wgpu::SurfaceError::OutOfMemory) => elwt.exit(),
                    Err(e) => eprintln!("Render error: {:?}", e),
                }
            }
            _ => {}
        },
        Event::AboutToWait => {
            window.request_redraw();
        }
        _ => {}
    });
}

fn main() {
    println!("Protein Viewer");
    println!("==============");
    println!("Controls:");
    println!("  Mouse drag: Rotate");
    println!("  Scroll: Zoom");
    println!("  R: Reset camera");
    println!("  Esc: Quit");
    println!();
    println!("Looking for scripts/init.lua or init.lua...");
    println!();

    pollster::block_on(run());
}
