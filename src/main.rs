//! Main entry point for the Protein Viewer application
//!
//! This module coordinates the interaction between the WGPU renderer,
//! the Lua script engine, and the protein data store

mod lua_api;
mod protein;
mod renderer;
mod selection;
mod analysis;
mod surface;

use clap::Parser;
use notify::Watcher;
use pdbtbx::ContainsAtomConformer;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, RwLock};
use winit::{
    event::*,
    event_loop::{ControlFlow, EventLoop},
    keyboard::{KeyCode, PhysicalKey},
    window::{WindowBuilder, WindowId},
};

use lua_api::{ScriptEngine, ScriptEvent};
use protein::{ColorScheme, ProteinStore, Representation};
use renderer::{Camera, Renderer};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the initial Lua script to run
    #[arg(short, long)]
    script: Option<String>,
}

/// State associated with a single window
struct WindowState {
    window: Arc<winit::window::Window>,
    surface: wgpu::Surface<'static>,
    config: wgpu::SurfaceConfiguration,
    renderer: Renderer,
    camera: Arc<RwLock<Camera>>,
    store: Arc<RwLock<ProteinStore>>,

    /// Selected atoms: (protein_name, atom_index_in_ca_positions)
    selected_atoms: Vec<(String, usize)>,
    /// Measurement pairs: (atom1, atom2) where atomX is index in selected_atoms
    measurements: Vec<(usize, usize)>,

    /// Whether the left mouse button is currently pressed
    mouse_pressed: bool,
    /// The last recorded mouse position for rotation calculations
    last_recorded_mouse_cursor_position: Option<(f64, f64)>,
}

impl WindowState {
    fn new(
        window: Arc<winit::window::Window>,
        instance: &wgpu::Instance,
        adapter: &wgpu::Adapter,
        device: Arc<wgpu::Device>,
        queue: Arc<wgpu::Queue>,
        store: Arc<RwLock<ProteinStore>>,
        camera: Arc<RwLock<Camera>>,
    ) -> Self {
        let window_inner_size = window.inner_size();
        let surface = instance.create_surface(window.clone()).unwrap();

        let surface_capabilities = surface.get_capabilities(adapter);
        let preferred_surface_format = surface_capabilities
            .formats
            .iter()
            .find(|f| f.is_srgb())
            .copied()
            .unwrap_or(surface_capabilities.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: preferred_surface_format,
            width: window_inner_size.width.max(1),
            height: window_inner_size.height.max(1),
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: surface_capabilities.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let renderer = Renderer::new(
            device.clone(),
            queue.clone(),
            preferred_surface_format,
            window_inner_size.width,
            window_inner_size.height,
        );

        // Ensure camera aspect ratio is correct for the window size
        camera.write().unwrap().set_aspect(window_inner_size.width as f32 / window_inner_size.height as f32);

        Self {
            window,
            surface,
            config,
            renderer,
            camera,
            store,
            selected_atoms: Vec::new(),
            measurements: Vec::new(),
            mouse_pressed: false,
            last_recorded_mouse_cursor_position: None,
        }
    }

    fn resize(&mut self, device: &wgpu::Device, physical_window_size: winit::dpi::PhysicalSize<u32>) {
        if physical_window_size.width > 0 && physical_window_size.height > 0 {
            self.config.width = physical_window_size.width;
            self.config.height = physical_window_size.height;
            self.surface.configure(device, &self.config);
            self.renderer
                .resize(physical_window_size.width, physical_window_size.height);
            self.camera
                .write()
                .unwrap()
                .set_aspect(physical_window_size.width as f32 / physical_window_size.height as f32);
        }
    }

    fn update(&mut self) {
        // Update renderer with current protein data
        {
            let locked_protein_store = self.store.read().unwrap();
            self.renderer.update_instances(
                &locked_protein_store,
                &self.selected_atoms,
                &self.measurements,
            );
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
                let (bounding_box_minimum, bounding_box_maximum) =
                    protein_locked_data.bounding_box();
                let bounding_sphere_radius =
                    (bounding_box_maximum - bounding_box_minimum).length() / 2.0;
                drop(protein_locked_data);
                camera.focus_on(protein_center_of_mass, bounding_sphere_radius);
            }
        }
    }

    fn render(&mut self, queue: &wgpu::Queue) -> Result<(), wgpu::SurfaceError> {
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
                let alpha_carbon_positions_and_chain_identifiers_collection =
                    protein_locked_data.get_alpha_carbon_positions_and_chain_identifiers();
                for (atom_indexing_counter, (atom_world_position, _chain_identifier)) in
                    alpha_carbon_positions_and_chain_identifiers_collection
                        .into_iter()
                        .enumerate()
                {
                    global_atom_positions_lookup_table.insert(
                        (protein_locked_data.name.clone(), atom_indexing_counter),
                        atom_world_position,
                    );
                }
            }

            for &(first_selected_atom_index, second_selected_atom_index) in &self.measurements {
                if let (Some(first_atom_selection_data), Some(second_atom_selection_data)) = (
                    self.selected_atoms.get(first_selected_atom_index),
                    self.selected_atoms.get(second_selected_atom_index),
                ) {
                    if let (Some(first_atom_world_position), Some(second_atom_world_position)) = (
                        global_atom_positions_lookup_table.get(first_atom_selection_data),
                        global_atom_positions_lookup_table.get(second_atom_selection_data),
                    ) {
                        let calculated_euclidean_distance =
                            first_atom_world_position.distance(*second_atom_world_position);
                        let midpoint_between_atoms =
                            (*first_atom_world_position + *second_atom_world_position) * 0.5;
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
                    text: format!(
                        "Protein: {} ({} atoms)",
                        protein_locked_data.name,
                        protein_locked_data.atom_count()
                    ),
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
            self.config.height,
        );
        queue.submit(std::iter::once(rendered_graphics_command_buffer));
        surface_texture_output.present();

        Ok(())
    }

    fn handle_mouse_input(
        &mut self,
        element_press_state: ElementState,
        mouse_button_identifier: MouseButton,
    ) {
        if mouse_button_identifier == MouseButton::Left {
            let was_previously_pressed = self.mouse_pressed;
            self.mouse_pressed = element_press_state == ElementState::Pressed;

            // Only perform ray casting on the initial press (click)
            if self.mouse_pressed && !was_previously_pressed {
                if let Some((current_mouse_x_coordinate, current_mouse_y_coordinate)) =
                    self.last_recorded_mouse_cursor_position
                {
                    let (ray_origin_point, ray_direction_unit_vector) =
                        self.camera.read().unwrap().calculate_ray_from_screen_coordinates(
                            glam::Vec2::new(
                                current_mouse_x_coordinate as f32,
                                current_mouse_y_coordinate as f32,
                            ),
                            glam::Vec2::new(self.config.width as f32, self.config.height as f32),
                        );

                    let mut closest_intersection_hit_data: Option<(String, usize, f32)> = None;

                    let locked_protein_store = self.store.read().unwrap();
                    for protein_shared_handle in locked_protein_store.iter() {
                        let protein_locked_data = protein_shared_handle.read().unwrap();
                        let protein_identifier_name = protein_locked_data.name.clone();
                        
                        let current_representation_mode = protein_locked_data.representation;
                        
                        // Use appropriate hit-box radius based on visual representation
                        let base_sphere_hitbox_radius = match current_representation_mode {
                            Representation::Spheres | Representation::BackboneAndSpheres => 1.5,
                            Representation::BallAndStick => 0.4,
                            Representation::SpaceFilling => 1.7,
                            _ => 1.0, // Fallback for modes without primary spheres
                        };

                        // Check ALL atoms for intersection, matching the renderer's indexing
                        for (atom_global_index, current_atom_hierarchy) in
                            protein_locked_data.pdb.atoms_with_hierarchy().into_iter().enumerate()
                        {
                            let current_atom_reference = current_atom_hierarchy.atom();
                            
                            // If in CA-only mode, skip non-CA atoms for consistency with visualization
                            if matches!(current_representation_mode, Representation::Spheres | Representation::BackboneAndSpheres) 
                               && current_atom_reference.name() != "CA" {
                                continue;
                            }

                            let atom_position_tuple = current_atom_reference.pos();
                            let atom_world_position_vector = glam::Vec3::new(
                                atom_position_tuple.0 as f32,
                                atom_position_tuple.1 as f32,
                                atom_position_tuple.2 as f32
                            );

                            let vector_from_ray_origin_to_atom_center =
                                ray_origin_point - atom_world_position_vector;
                            
                            let quadratic_coefficient_b = vector_from_ray_origin_to_atom_center
                                .dot(ray_direction_unit_vector);
                            let quadratic_coefficient_c = vector_from_ray_origin_to_atom_center
                                .dot(vector_from_ray_origin_to_atom_center)
                                - base_sphere_hitbox_radius * base_sphere_hitbox_radius;
                            
                            let intersection_discriminant_value = quadratic_coefficient_b
                                * quadratic_coefficient_b
                                - quadratic_coefficient_c;

                            if intersection_discriminant_value >= 0.0 {
                                let distance_to_intersection_point = -quadratic_coefficient_b
                                    - intersection_discriminant_value.sqrt();
                                if distance_to_intersection_point > 0.0 {
                                    if closest_intersection_hit_data.is_none()
                                        || distance_to_intersection_point
                                            < closest_intersection_hit_data.as_ref().unwrap().2
                                    {
                                        closest_intersection_hit_data = Some((
                                            protein_identifier_name.clone(),
                                            atom_global_index,
                                            distance_to_intersection_point,
                                        ));
                                    }
                                }
                            }
                        }
                    }

                    if let Some((hit_protein_name, hit_atom_index, _)) =
                        closest_intersection_hit_data
                    {
                        let target_atom_selection_identifier = (hit_protein_name, hit_atom_index);
                        if self.selected_atoms.contains(&target_atom_selection_identifier) {
                            // Toggle selection off if already selected
                            self.selected_atoms.retain(|item| item != &target_atom_selection_identifier);
                        } else {
                            self.selected_atoms.push(target_atom_selection_identifier);
                        }
                    }
                }
            }
        }
    }

    fn handle_mouse_move(&mut self, cursor_screen_position: (f64, f64)) {
        if self.mouse_pressed {
            if let Some((previous_mouse_x, previous_mouse_y)) = self.last_recorded_mouse_cursor_position {
                let mouse_movement_delta_x =
                    (cursor_screen_position.0 - previous_mouse_x) as f32 * 0.01;
                let mouse_movement_delta_y =
                    (cursor_screen_position.1 - previous_mouse_y) as f32 * 0.01;
                self.camera
                    .write()
                    .unwrap()
                    .rotate(-mouse_movement_delta_x, mouse_movement_delta_y);
            }
        }
        self.last_recorded_mouse_cursor_position = Some(cursor_screen_position);
    }

    fn handle_scroll(&mut self, scroll_delta_value: f32) {
        self.camera.write().unwrap().zoom(scroll_delta_value);
    }

    fn handle_key(&mut self, pressed_physical_keycode: KeyCode) {
        match pressed_physical_keycode {
            // Camera and Selection Control
            KeyCode::KeyR => {
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
                    let (bounding_box_minimum, bounding_box_maximum) =
                        protein_locked_data.bounding_box();
                    let bounding_sphere_radius =
                        (bounding_box_maximum - bounding_box_minimum).length() / 2.0;
                    drop(protein_locked_data);
                    camera.focus_on(protein_center_of_mass, bounding_sphere_radius);
                }
            }
            
            // Representation Modes
            KeyCode::Digit1 => self.set_all_representation(Representation::Spheres),
            KeyCode::Digit2 => self.set_all_representation(Representation::Backbone),
            KeyCode::Digit3 => self.set_all_representation(Representation::BackboneAndSpheres),
            KeyCode::Digit4 => self.set_all_representation(Representation::Sticks),
            KeyCode::Digit5 => self.set_all_representation(Representation::BallAndStick),
            KeyCode::Digit6 => self.set_all_representation(Representation::SpaceFilling),
            KeyCode::Digit7 => self.set_all_representation(Representation::Lines),
            
            // Color Schemes
            KeyCode::KeyC => self.set_all_color_scheme(ColorScheme::ByChain),
            KeyCode::KeyB => self.set_all_color_scheme(ColorScheme::ByBFactor),
            KeyCode::KeyE => self.set_all_color_scheme(ColorScheme::ByElement),
            KeyCode::KeyS => self.set_all_color_scheme(ColorScheme::BySecondary),
            
            // Analysis and Surface
            KeyCode::KeyM => {
                if self.selected_atoms.len() >= 2 {
                    let second_atom_selection_index = self.selected_atoms.len() - 1;
                    let first_atom_selection_index = self.selected_atoms.len() - 2;
                    self.measurements
                        .push((first_atom_selection_index, second_atom_selection_index));
                }
            }
            KeyCode::KeyF => {
                // Toggle surface for all proteins (placeholder for a more complex toggle)
                let locked_protein_store = self.store.read().unwrap();
                for protein_shared_handle in locked_protein_store.iter() {
                    let mut protein_mutable_data = protein_shared_handle.write().unwrap();
                    let current_visibility = protein_mutable_data.molecular_surface_mesh.is_surface_visible;
                    protein_mutable_data.molecular_surface_mesh.is_surface_visible = !current_visibility;
                }
            }
            
            KeyCode::Escape => std::process::exit(0),
            _ => {}
        }
    }

    fn set_all_representation(&self, target_representation_mode: Representation) {
        let locked_protein_store = self.store.read().unwrap();
        for protein_shared_handle in locked_protein_store.iter() {
            let mut protein_mutable_data = protein_shared_handle.write().unwrap();
            protein_mutable_data.representation = target_representation_mode;
        }
    }

    fn set_all_color_scheme(&self, target_color_scheme: ColorScheme) {
        let locked_protein_store = self.store.read().unwrap();
        for protein_shared_handle in locked_protein_store.iter() {
            let mut protein_mutable_data = protein_shared_handle.write().unwrap();
            protein_mutable_data.color_scheme = target_color_scheme;
        }
    }
}

/// The main application state managing multiple windows
struct App {
    instance: wgpu::Instance,
    adapter: wgpu::Adapter,
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    windows: HashMap<WindowId, WindowState>,
    script_engine: ScriptEngine,
    global_store: Arc<RwLock<ProteinStore>>,
    global_camera: Arc<RwLock<Camera>>,
    script_rx: crossbeam_channel::Receiver<PathBuf>,
    script_event_rx: crossbeam_channel::Receiver<ScriptEvent>,
    _watcher: notify::RecommendedWatcher,
}

impl App {
    async fn new() -> Self {
        let wgpu_instance = wgpu::Instance::default();
        let graphics_adapter = wgpu_instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
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

        let (script_event_tx, script_event_rx) = crossbeam_channel::unbounded();
        
        let global_store = Arc::new(RwLock::new(ProteinStore::new()));
        let global_camera = Arc::new(RwLock::new(Camera::new(1.0)));
        
        let lua_script_engine = ScriptEngine::new(global_store.clone(), global_camera.clone(), script_event_tx)
            .expect("Failed to create Lua engine");

        let (script_path_sender, script_path_receiver) = crossbeam_channel::unbounded::<PathBuf>();
        let mut file_system_watcher =
            notify::recommended_watcher(move |watcher_result: notify::Result<notify::Event>| {
                if let Ok(notify_event) = watcher_result {
                    if notify_event.kind.is_modify() {
                        for modified_file_path in notify_event.paths {
                            if modified_file_path
                                .extension()
                                .map(|extension| extension == "lua")
                                .unwrap_or(false)
                            {
                                let _ = script_path_sender.send(modified_file_path);
                            }
                        }
                    }
                }
            })
            .expect("Failed to create watcher");

        let _ = file_system_watcher.watch(
            std::path::Path::new("scripts"),
            notify::RecursiveMode::Recursive,
        );
        let _ = file_system_watcher.watch(
            std::path::Path::new("."),
            notify::RecursiveMode::NonRecursive,
        );

        Self {
            instance: wgpu_instance,
            adapter: graphics_adapter,
            device: graphics_device,
            queue: command_queue,
            windows: HashMap::new(),
            script_engine: lua_script_engine,
            global_store,
            global_camera,
            script_rx: script_path_receiver,
            script_event_rx,
            _watcher: file_system_watcher,
        }
    }

    fn handle_script_events(&mut self, event_loop: &winit::event_loop::EventLoopWindowTarget<()>) {
        while let Ok(event) = self.script_event_rx.try_recv() {
            match event {
                ScriptEvent::NewWindow(store, camera) => {
                    let window = Arc::new(
                        WindowBuilder::new()
                            .with_title("Protein Viewer")
                            .with_inner_size(winit::dpi::LogicalSize::new(1280, 800))
                            .build(event_loop)
                            .unwrap(),
                    );
                    let window_id = window.id();
                    let window_state = WindowState::new(
                        window,
                        &self.instance,
                        &self.adapter,
                        self.device.clone(),
                        self.queue.clone(),
                        store,
                        camera,
                    );
                    self.windows.insert(window_id, window_state);
                }
            }
        }
    }

    fn update(&mut self, event_loop: &winit::event_loop::EventLoopWindowTarget<()>) {
        // Debounce script hot-reload: drain the channel and only keep the last path
        let mut last_reloaded_path = None;
        while let Ok(reloaded_script_path) = self.script_rx.try_recv() {
            last_reloaded_path = Some(reloaded_script_path);
        }

        if let Some(reloaded_script_path) = last_reloaded_path {
            println!("Reloading: {:?}", reloaded_script_path);
            if let Some(reloaded_script_path_string) = reloaded_script_path.to_str() {
                // DISCARD any pending window creation events from the previous run
                while self.script_event_rx.try_recv().is_ok() {}

                // RESET STATE: Close all windows and clear global store on reload
                self.windows.clear();
                self.global_store.write().unwrap().clear();

                if let Err(hot_reload_error) =
                    self.script_engine.run_file(reloaded_script_path_string)
                {
                    eprintln!("Script error: {}", hot_reload_error);
                }
            }
        }

        self.handle_script_events(event_loop);

        for window_state in self.windows.values_mut() {
            window_state.update();
        }
    }
}

async fn run(initial_script: Option<String>) {
    let main_event_loop = EventLoop::new().unwrap();
    let mut app = App::new().await;

    // Run the initial script if provided by the user
    if let Some(provided_initial_script_path) = initial_script {
        if let Err(script_execution_error) = app.script_engine.run_file(&provided_initial_script_path) {
            eprintln!("Error running explicitly provided script {}: {}", provided_initial_script_path, script_execution_error);
        }
    }

    main_event_loop.set_control_flow(ControlFlow::Poll);
    let _ = main_event_loop.run(
        move |winit_event, event_loop_window_target| {
            match winit_event {
                Event::WindowEvent { window_id, event } => {
                    if let Some(window_state) = app.windows.get_mut(&window_id) {
                        match event {
                            WindowEvent::CloseRequested => {
                                app.windows.remove(&window_id);
                                if app.windows.is_empty() {
                                    // Optionally exit if last window closed, 
                                    // but we might want to keep the script engine running for new windows
                                }
                            },
                            WindowEvent::Resized(new_window_size) => window_state.resize(&app.device, new_window_size),
                            WindowEvent::MouseInput {
                                state,
                                button,
                                ..
                            } => window_state.handle_mouse_input(state, button),
                            WindowEvent::CursorMoved {
                                position,
                                ..
                            } => window_state.handle_mouse_move((position.x, position.y)),
                            WindowEvent::MouseWheel {
                                delta,
                                ..
                            } => {
                                let scroll = match delta {
                                    MouseScrollDelta::LineDelta(_, y) => y,
                                    MouseScrollDelta::PixelDelta(pos) => pos.y as f32 / 100.0,
                                };
                                window_state.handle_scroll(scroll);
                            }
                            WindowEvent::KeyboardInput {
                                event:
                                    KeyEvent {
                                        physical_key: PhysicalKey::Code(code),
                                        state: ElementState::Pressed,
                                        ..
                                    },
                                ..
                            } => window_state.handle_key(code),
                            WindowEvent::RedrawRequested => {
                                match window_state.render(&app.queue) {
                                    Ok(_) => {}
                                    Err(wgpu::SurfaceError::Lost) => {
                                        window_state.resize(&app.device, window_state.window.inner_size())
                                    }
                                    Err(wgpu::SurfaceError::OutOfMemory) => {
                                        event_loop_window_target.exit();
                                    }
                                    Err(e) => eprintln!("Render error: {:?}", e),
                                }
                            }
                            _ => {}
                        }
                    }
                },
                Event::AboutToWait => {
                    app.update(event_loop_window_target);
                    for window_state in app.windows.values() {
                        window_state.window.request_redraw();
                    }
                }
                _ => {}
            }
        },
    );
}

fn main() {

    let args = Args::parse();

    pollster::block_on(run(args.script));

}
