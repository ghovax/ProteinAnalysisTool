//! Main entry point for the Protein Viewer application
//!
//! This module coordinates the interaction between the WGPU renderer,
//! the Lua script engine, and the protein data store

mod lua_api;
mod protein;
mod renderer;

use clap::Parser;
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
use protein::{ProteinStore, Representation, ColorScheme};
use renderer::{Camera, Renderer};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the initial Lua script to run
    #[arg(short, long)]
    script: Option<String>,
}

/// The main application state
struct App {
    /// WGPU surface for rendering
    surface: wgpu::Surface<'static>,
    /// The GPU device handle
    device: Arc<wgpu::Device>,
    /// The GPU command queue
    queue: Arc<wgpu::Queue>,
    /// Current surface configuration
    config: wgpu::SurfaceConfiguration,
    /// The protein renderer
    renderer: Renderer,
    /// The interactive camera
    camera: Arc<RwLock<Camera>>,
    /// The Lua script engine for automation
    script_engine: ScriptEngine,
    /// The shared store containing all loaded proteins
    store: Arc<RwLock<ProteinStore>>,

    /// Selected atoms: (protein_name, atom_index_in_ca_positions)
    selected_atoms: Vec<(String, usize)>,
    /// Measurement pairs: (atom1, atom2) where atomX is index in selected_atoms
    measurements: Vec<(usize, usize)>,

    /// Whether the left mouse button is currently pressed
    mouse_pressed: bool,
    /// The last recorded mouse position for rotation calculations
    last_mouse_pos: Option<(f64, f64)>,

    /// Receiver for script hot-reload notifications
    script_rx: crossbeam_channel::Receiver<PathBuf>,
    /// The file system watcher for scripts
    _watcher: notify::RecommendedWatcher,
}

impl App {
    /// Initializes a new instance of the application
    async fn new(main_window: Arc<winit::window::Window>, initial_script: Option<String>) -> Self {
        let window_inner_size = main_window.inner_size();
        let wgpu_instance = wgpu::Instance::default();
        let rendering_surface = wgpu_instance.create_surface(main_window.clone()).unwrap();
        let graphics_adapter = wgpu_instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&rendering_surface),
                force_fallback_adapter: false,
            })
            .await
            .expect("Failed to find adapter");

        let (graphics_device, command_queue) = graphics_adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .expect("Failed to create device");

        let graphics_device = Arc::new(graphics_device);
        let command_queue = Arc::new(command_queue);

        let surface_capabilities = rendering_surface.get_capabilities(&graphics_adapter);
        let preferred_surface_format = surface_capabilities
            .formats
            .iter()
            .find(|f| f.is_srgb())
            .copied()
            .unwrap_or(surface_capabilities.formats[0]);

        let surface_configuration = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: preferred_surface_format,
            width: window_inner_size.width.max(1),
            height: window_inner_size.height.max(1),
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: surface_capabilities.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        rendering_surface.configure(&graphics_device, &surface_configuration);

        let main_renderer = Renderer::new(
            graphics_device.clone(),
            command_queue.clone(),
            preferred_surface_format,
            window_inner_size.width,
            window_inner_size.height,
        );

        let main_camera = Arc::new(RwLock::new(Camera::new(window_inner_size.width as f32 / window_inner_size.height as f32)));

        let protein_data_store = Arc::new(RwLock::new(ProteinStore::new()));
        let lua_script_engine = ScriptEngine::new(protein_data_store.clone(), main_camera.clone()).expect("Failed to create Lua engine");

        // Set up file watcher for hot-reload
        let (script_path_sender, script_path_receiver) = crossbeam_channel::unbounded::<PathBuf>();
        let mut file_system_watcher =
            notify::recommended_watcher(move |watcher_result: notify::Result<notify::Event>| {
                if let Ok(notify_event) = watcher_result {
                    if notify_event.kind.is_modify() {
                        for modified_file_path in notify_event.paths {
                            if modified_file_path.extension().map(|extension| extension == "lua").unwrap_or(false) {
                                let _ = script_path_sender.send(modified_file_path);
                            }
                        }
                    }
                }
            })
            .expect("Failed to create watcher");

        // Watch scripts directory and current directory
        let _ = file_system_watcher.watch(
            std::path::Path::new("scripts"),
            notify::RecursiveMode::Recursive,
        );
        let _ = file_system_watcher.watch(
            std::path::Path::new("."),
            notify::RecursiveMode::NonRecursive,
        );

        // Run initial script if provided
        if let Some(script_path) = initial_script {
            if let Err(script_error_message) = lua_script_engine.run_file(&script_path) {
                eprintln!("Script error: {}", script_error_message);
            }
        }

        Self {
            surface: rendering_surface,
            device: graphics_device,
            queue: command_queue,
            config: surface_configuration,
            renderer: main_renderer,
            camera: main_camera,
            script_engine: lua_script_engine,
            store: protein_data_store,
            selected_atoms: Vec::new(),
            measurements: Vec::new(),
            mouse_pressed: false,
            last_mouse_pos: None,
            script_rx: script_path_receiver,
            _watcher: file_system_watcher,
        }
    }

    /// Resizes the application surface and updates the renderer and camera
    fn resize(&mut self, physical_window_size: winit::dpi::PhysicalSize<u32>) {
        if physical_window_size.width > 0 && physical_window_size.height > 0 {
            self.config.width = physical_window_size.width;
            self.config.height = physical_window_size.height;
            self.surface.configure(&self.device, &self.config);
            self.renderer.resize(physical_window_size.width, physical_window_size.height);
            self.camera.write().unwrap()
                .set_aspect(physical_window_size.width as f32 / physical_window_size.height as f32);
        }
    }

    /// Updates the application state, including script hot-reloads and camera focus
    fn update(&mut self) {
        // Check for script hot-reload
        while let Ok(reloaded_script_path) = self.script_rx.try_recv() {
            println!("Reloading: {:?}", reloaded_script_path);
            if let Some(reloaded_script_path_string) = reloaded_script_path.to_str() {
                if let Err(hot_reload_error) = self.script_engine.run_file(reloaded_script_path_string) {
                    eprintln!("Script error: {}", hot_reload_error);
                }
            }
        }

        // Update renderer with current protein data
        {
            let locked_protein_store = self.store.read().unwrap();
            self.renderer.update_instances(&locked_protein_store, &self.selected_atoms, &self.measurements);
        }

        // Auto-focus camera on first protein if we haven't moved it
        let mut camera = self.camera.write().unwrap();
        if camera.target == glam::Vec3::ZERO {
            let locked_protein_store = self.store.read().unwrap();
            let optional_first_protein_entry = locked_protein_store.iter().next().cloned();
            drop(locked_protein_store);

            if let Some(protein_shared_reference) = optional_first_protein_entry {
                let protein_locked_data = protein_shared_reference.read().unwrap();
                let protein_center_of_mass = protein_locked_data.center_of_mass();
                let (bounding_box_minimum, bounding_box_maximum) = protein_locked_data.bounding_box();
                let bounding_sphere_radius = (bounding_box_maximum - bounding_box_minimum).length() / 2.0;
                drop(protein_locked_data);
                camera.focus_on(protein_center_of_mass, bounding_sphere_radius);
            }
        }
    }

    /// Renders the current frame to the surface
    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let surface_texture_output = self.surface.get_current_texture()?;
        let texture_view_for_rendering = surface_texture_output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        let mut text_labels_collection = Vec::new();
        
        // Add distance measurement labels
        {
            let locked_protein_store = self.store.read().unwrap();
            let mut global_atom_positions_lookup_table = std::collections::HashMap::new();
            for protein_shared_handle in locked_protein_store.iter() {
                let protein_locked_data = protein_shared_handle.read().unwrap();
                let alpha_carbon_positions_and_chain_identifiers_collection = protein_locked_data.get_alpha_carbon_positions_and_chain_identifiers();
                for (atom_indexing_counter, (atom_world_position, _chain_identifier)) in alpha_carbon_positions_and_chain_identifiers_collection.into_iter().enumerate() {
                    global_atom_positions_lookup_table.insert((protein_locked_data.name.clone(), atom_indexing_counter), atom_world_position);
                }
            }

            for &(first_selected_atom_index, second_selected_atom_index) in &self.measurements {
                if let (Some(first_atom_selection_data), Some(second_atom_selection_data)) = (self.selected_atoms.get(first_selected_atom_index), self.selected_atoms.get(second_selected_atom_index)) {
                    if let (Some(first_atom_world_position), Some(second_atom_world_position)) = (global_atom_positions_lookup_table.get(first_atom_selection_data), global_atom_positions_lookup_table.get(second_atom_selection_data)) {
                        let calculated_euclidean_distance = first_atom_world_position.distance(*second_atom_world_position);
                        let midpoint_between_atoms = (*first_atom_world_position + *second_atom_world_position) * 0.5;
                        text_labels_collection.push(renderer::pipeline::TextLabel {
                            position: midpoint_between_atoms,
                            text: format!("{:.2} A", calculated_euclidean_distance),
                            color: [1.0, 1.0, 1.0, 1.0],
                        });
                    }
                }
            }
        }

        // Add protein information labels for the UI overlay
        {
            let locked_protein_store = self.store.read().unwrap();
            for protein_shared_handle in locked_protein_store.iter() {
                let protein_locked_data = protein_shared_handle.read().unwrap();
                text_labels_collection.push(renderer::pipeline::TextLabel {
                    position: glam::Vec3::ZERO, // UI labels use special screen-space positioning
                    text: format!("Protein: {} ({} atoms)", protein_locked_data.name, protein_locked_data.atom_count()),
                    color: [1.0, 1.0, 1.0, 1.0],
                });
            }
        }

        let locked_camera_handle = self.camera.read().unwrap();
        let rendered_graphics_command_buffer = self.renderer.render(
            &texture_view_for_rendering, 
            &locked_camera_handle, 
            &text_labels_collection,
            self.config.width,
            self.config.height
        );
        self.queue.submit(std::iter::once(rendered_graphics_command_buffer));
        surface_texture_output.present();

        Ok(())
    }

    /// Handles mouse button input events
    fn handle_mouse_input(&mut self, element_press_state: ElementState, mouse_button_identifier: MouseButton) {
        if mouse_button_identifier == MouseButton::Left {
            self.mouse_pressed = element_press_state == ElementState::Pressed;
            
            if self.mouse_pressed {
                // Ray casting for selection
                if let Some((current_mouse_x_coordinate, current_mouse_y_coordinate)) = self.last_mouse_pos {
                    let (ray_origin_point, ray_direction_unit_vector) = self.camera.read().unwrap().ray_from_screen(
                        glam::Vec2::new(current_mouse_x_coordinate as f32, current_mouse_y_coordinate as f32),
                        glam::Vec2::new(self.config.width as f32, self.config.height as f32)
                    );

                    let mut closest_intersection_hit_data: Option<(String, usize, f32)> = None;

                    let locked_protein_store = self.store.read().unwrap();
                    for protein_shared_handle in locked_protein_store.iter() {
                        let protein_locked_data = protein_shared_handle.read().unwrap();
                        let protein_identifier_name = protein_locked_data.name.clone();
                        
                        // Check Alpha Carbon (CA) atoms for intersection
                        let alpha_carbon_positions_and_chain_identifiers_collection = protein_locked_data.get_alpha_carbon_positions_and_chain_identifiers();
                        for (atom_indexing_counter, (atom_world_position, _chain_identifier)) in alpha_carbon_positions_and_chain_identifiers_collection.into_iter().enumerate() {
                            let sphere_intersection_radius = 1.5; // Same radius as used in the renderer
                            let vector_from_ray_origin_to_atom_center = ray_origin_point - atom_world_position;
                            let quadratic_coefficient_b = vector_from_ray_origin_to_atom_center.dot(ray_direction_unit_vector);
                            let quadratic_coefficient_c = vector_from_ray_origin_to_atom_center.dot(vector_from_ray_origin_to_atom_center) - sphere_intersection_radius * sphere_intersection_radius;
                            let intersection_discriminant_value = quadratic_coefficient_b * quadratic_coefficient_b - quadratic_coefficient_c;
                            
                            if intersection_discriminant_value >= 0.0 {
                                let distance_to_intersection_point = -quadratic_coefficient_b - intersection_discriminant_value.sqrt();
                                if distance_to_intersection_point > 0.0 {
                                    if closest_intersection_hit_data.is_none() || distance_to_intersection_point < closest_intersection_hit_data.as_ref().unwrap().2 {
                                        closest_intersection_hit_data = Some((protein_identifier_name.clone(), atom_indexing_counter, distance_to_intersection_point));
                                    }
                                }
                            }
                        }
                    }

                    if let Some((hit_protein_name, hit_atom_index, _)) = closest_intersection_hit_data {
                        let target_atom_selection_identifier = (hit_protein_name, hit_atom_index);
                        if !self.selected_atoms.contains(&target_atom_selection_identifier) {
                            self.selected_atoms.push(target_atom_selection_identifier);
                        }
                    }
                }
            }

            if !self.mouse_pressed {
                self.last_mouse_pos = None;
            }
        }
    }

    /// Handles cursor movement events for camera rotation
    fn handle_mouse_move(&mut self, cursor_screen_position: (f64, f64)) {
        if self.mouse_pressed {
            if let Some((previous_mouse_x, previous_mouse_y)) = self.last_mouse_pos {
                let mouse_movement_delta_x = (cursor_screen_position.0 - previous_mouse_x) as f32 * 0.01;
                let mouse_movement_delta_y = (cursor_screen_position.1 - previous_mouse_y) as f32 * 0.01;
                self.camera.write().unwrap().rotate(-mouse_movement_delta_x, mouse_movement_delta_y);
            }
        }
        self.last_mouse_pos = Some(cursor_screen_position);
    }

    /// Handles mouse scroll events for camera zoom
    fn handle_scroll(&mut self, scroll_delta_value: f32) {
        self.camera.write().unwrap().zoom(scroll_delta_value);
    }

    /// Handles keyboard input events for various shortcuts
    fn handle_key(&mut self, pressed_physical_keycode: KeyCode) {
        match pressed_physical_keycode {
            KeyCode::KeyR => {
                // Reset camera and selection
                self.selected_atoms.clear();
                self.measurements.clear();
                let mut camera = self.camera.write().unwrap();
                *camera = Camera::new(self.config.width as f32 / self.config.height as f32);
                let locked_protein_store = self.store.read().unwrap();
                let optional_first_protein_entry = locked_protein_store.iter().next().cloned();
                drop(locked_protein_store);

                if let Some(protein_shared_reference) = optional_first_protein_entry {
                    let protein_locked_data = protein_shared_reference.read().unwrap();
                    let protein_center_of_mass = protein_locked_data.center_of_mass();
                    let (bounding_box_minimum, bounding_box_maximum) = protein_locked_data.bounding_box();
                    let bounding_sphere_radius = (bounding_box_maximum - bounding_box_minimum).length() / 2.0;
                    drop(protein_locked_data);
                    camera.focus_on(protein_center_of_mass, bounding_sphere_radius);
                }
            }
            // Representation modes
            KeyCode::Digit1 => {
                self.set_all_representation(Representation::Spheres);
            }
            KeyCode::Digit2 => {
                self.set_all_representation(Representation::Backbone);
            }
            KeyCode::Digit3 => {
                self.set_all_representation(Representation::BackboneAndSpheres);
            }
            // Color schemes
            KeyCode::KeyC => {
                self.set_all_color_scheme(ColorScheme::ByChain);
            }
            KeyCode::KeyB => {
                self.set_all_color_scheme(ColorScheme::ByBFactor);
            }
            KeyCode::KeyE => {
                self.set_all_color_scheme(ColorScheme::ByElement);
            }
            KeyCode::KeyS => {
                self.set_all_color_scheme(ColorScheme::BySecondary);
            }
            KeyCode::KeyM => {
                // Calculate distance between the last two selected atoms
                if self.selected_atoms.len() >= 2 {
                    let second_atom_selection_index = self.selected_atoms.len() - 1;
                    let first_atom_selection_index = self.selected_atoms.len() - 2;
                    self.measurements.push((first_atom_selection_index, second_atom_selection_index));
                }
            }
            KeyCode::Escape => {
                std::process::exit(0);
            }
            _ => {}
        }
    }

    /// Sets the representation mode for all loaded proteins
    fn set_all_representation(&self, target_representation_mode: Representation) {
        let locked_protein_store = self.store.read().unwrap();
        for protein_shared_handle in locked_protein_store.iter() {
            let mut protein_mutable_data = protein_shared_handle.write().unwrap();
            protein_mutable_data.representation = target_representation_mode;
        }
    }

    /// Sets the color scheme for all loaded proteins
    fn set_all_color_scheme(&self, target_color_scheme: ColorScheme) {
        let locked_protein_store = self.store.read().unwrap();
        for protein_shared_handle in locked_protein_store.iter() {
            let mut protein_mutable_data = protein_shared_handle.write().unwrap();
            protein_mutable_data.color_scheme = target_color_scheme;
        }
    }
}



/// The main application event loop
async fn run(initial_script: Option<String>) {
    let main_event_loop = EventLoop::new().unwrap();
    let application_window = Arc::new(
        WindowBuilder::new()
            .with_title("Protein Viewer")
            .with_inner_size(winit::dpi::LogicalSize::new(1280, 800))
            .build(&main_event_loop)
            .unwrap(),
    );

    let mut application_state = App::new(application_window.clone(), initial_script).await;

    main_event_loop.set_control_flow(ControlFlow::Poll);
    let _ = main_event_loop.run(move |winit_event, event_loop_window_target| match winit_event {
        Event::WindowEvent { event, .. } => match event {
            WindowEvent::CloseRequested => event_loop_window_target.exit(),
            WindowEvent::Resized(new_window_size) => application_state.resize(new_window_size),
            WindowEvent::MouseInput { state: element_press_state, button: mouse_button_identifier, .. } => application_state.handle_mouse_input(element_press_state, mouse_button_identifier),
            WindowEvent::CursorMoved { position: cursor_screen_position, .. } => {
                application_state.handle_mouse_move((cursor_screen_position.x, cursor_screen_position.y))
            }
            WindowEvent::MouseWheel { delta: scroll_delta_value, .. } => {
                let calculated_scroll_amount = match scroll_delta_value {
                    MouseScrollDelta::LineDelta(_, y) => y,
                    MouseScrollDelta::PixelDelta(pos) => pos.y as f32 / 100.0,
                };
                application_state.handle_scroll(calculated_scroll_amount);
            }
            WindowEvent::KeyboardInput {
                event:
                    KeyEvent {
                        physical_key: PhysicalKey::Code(pressed_physical_keycode),
                        state: ElementState::Pressed,
                        ..
                    },
                ..
            } => application_state.handle_key(pressed_physical_keycode),
            WindowEvent::RedrawRequested => {
                application_state.update();
                match application_state.render() {
                    Ok(_) => {}
                    Err(wgpu::SurfaceError::Lost) => application_state.resize(winit::dpi::PhysicalSize {
                        width: application_state.config.width,
                        height: application_state.config.height,
                    }),
                    Err(wgpu::SurfaceError::OutOfMemory) => event_loop_window_target.exit(),
                    Err(wgpu_render_error) => eprintln!("Render error: {:?}", wgpu_render_error),
                }
            }
            _ => {}
        },
        Event::AboutToWait => {
            application_window.request_redraw();
        }
        _ => {}
    });
}

/// Application entry point
fn main() {
    let args = Args::parse();
    pollster::block_on(run(args.script));
}