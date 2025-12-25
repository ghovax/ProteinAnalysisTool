//! Main entry point for the Protein Viewer application
//!
//! This module coordinates the interaction between the WGPU renderer,
//! the Lua script engine, and the protein data store

mod analysis;
mod lua_api;
mod protein;
mod renderer;
mod selection;
mod surface;

use clap::Parser;
use notify::Watcher;
use pdbtbx::ContainsAtomConformer;
use std::path::PathBuf;
use std::sync::{Arc, RwLock};
use winit::{
    event::*,
    event_loop::{ControlFlow, EventLoop},
    keyboard::{KeyCode, PhysicalKey},
    window::WindowBuilder,
};

use lua_api::{ScriptEngine, ScriptEvent};
use protein::{ColorScheme, ProteinStore, Representation};
use renderer::{Camera, Renderer};
use tracing::{info, error};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the initial Lua script to run
    #[arg(short, long)]
    script: Option<String>,
}

/// State for a single rendering viewport
struct ViewportState {
    camera: Arc<RwLock<Camera>>,
    store: Arc<RwLock<ProteinStore>>,
    /// Selected atoms: (protein_name, atom_index_in_pdb_hierarchy)
    selected_atoms: Vec<(String, usize)>,
    /// Measurement pairs: (atom1, atom2) where atomX is index in selected_atoms
    measurements: Vec<(usize, usize)>,
    /// Whether the distance for the current selection has already been printed
    distance_has_been_printed_for_current_selection: bool,
}

impl ViewportState {
    fn new(store: Arc<RwLock<ProteinStore>>, camera: Arc<RwLock<Camera>>) -> Self {
        Self {
            store,
            camera,
            selected_atoms: Vec::new(),
            measurements: Vec::new(),
            distance_has_been_printed_for_current_selection: false,
        }
    }
}

/// State associated with a single window
struct WindowState {
    window: Arc<winit::window::Window>,
    surface: wgpu::Surface<'static>,
    config: wgpu::SurfaceConfiguration,
    renderer: Renderer,
    
    /// Multiple viewports displayed in this window
    viewports: Vec<ViewportState>,

    /// The index of the viewport currently being interacted with (captured during drag)
    captured_viewport_index: Option<usize>,

    /// Whether the left mouse button is currently pressed
    mouse_pressed: bool,
    /// Whether the right mouse button is currently pressed
    right_mouse_button_is_pressed: bool,
    /// The last recorded mouse position for rotation calculations
    last_recorded_mouse_cursor_position: Option<(f64, f64)>,

    /// The path of the currently active Lua script
    active_script_path_identifier: String,
}

impl WindowState {
    fn new(
        window: Arc<winit::window::Window>,
        instance: &wgpu::Instance,
        adapter: &wgpu::Adapter,
        device: Arc<wgpu::Device>,
        queue: Arc<wgpu::Queue>,
        script_path: String,
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

        let mut window_state_instance = Self {
            window,
            surface,
            config,
            renderer,
            viewports: Vec::new(),
            captured_viewport_index: None,
            mouse_pressed: false,
            right_mouse_button_is_pressed: false,
            last_recorded_mouse_cursor_position: None,
            active_script_path_identifier: script_path,
        };
        
        window_state_instance.update_window_title_with_current_state();
        window_state_instance
    }

    fn update_window_title_with_current_state(&mut self) {
        let mut loaded_protein_identifiers_collection = Vec::new();
        for viewport in &self.viewports {
            let locked_protein_store = viewport.store.read().unwrap();
            for protein_name in locked_protein_store.list() {
                if !loaded_protein_identifiers_collection.contains(&protein_name) {
                    loaded_protein_identifiers_collection.push(protein_name);
                }
            }
        }

        let formatted_window_title = if loaded_protein_identifiers_collection.is_empty() {
            format!("Protein Viewer (Script: {})", self.active_script_path_identifier)
        } else {
            format!(
                "Protein Viewer (Script: {}, Proteins: {})",
                self.active_script_path_identifier,
                loaded_protein_identifiers_collection.join(", ")
            )
        };
        
        self.window.set_title(&formatted_window_title);
    }

    fn add_viewport(&mut self, store: Arc<RwLock<ProteinStore>>, camera: Arc<RwLock<Camera>>) {
        self.viewports.push(ViewportState::new(store, camera));
        self.recalculate_viewport_aspect_ratios();
        self.update_window_title_with_current_state();
    }

    fn calculate_grid_dimensions(&self) -> (u32, u32) {
        let total_viewport_count = self.viewports.len() as u32;
        if total_viewport_count == 0 {
            return (0, 0);
        }
        
        let column_count = (total_viewport_count as f32).sqrt().ceil() as u32;
        let row_count = (total_viewport_count as f32 / column_count as f32).ceil() as u32;
        
        (column_count, row_count)
    }

    fn recalculate_viewport_aspect_ratios(&mut self) {
        let (column_count, row_count) = self.calculate_grid_dimensions();
        if column_count == 0 || row_count == 0 {
            return;
        }

        let window_width = self.config.width as f32;
        let window_height = self.config.height as f32;
        
        let viewport_width = window_width / column_count as f32;
        let viewport_height = window_height / row_count as f32;
        let aspect_ratio = viewport_width / viewport_height;

        for viewport in &self.viewports {
            viewport.camera.write().unwrap().set_aspect(aspect_ratio);
        }
    }

    fn resize(
        &mut self,
        device: &wgpu::Device,
        physical_window_size: winit::dpi::PhysicalSize<u32>,
    ) {
        if physical_window_size.width > 0 && physical_window_size.height > 0 {
            self.config.width = physical_window_size.width;
            self.config.height = physical_window_size.height;
            self.surface.configure(device, &self.config);
            self.renderer
                .resize(physical_window_size.width, physical_window_size.height);
            self.recalculate_viewport_aspect_ratios();
        }
    }

    fn update(&mut self) {
        // Auto-focus logic for each viewport if needed
        for viewport in &self.viewports {
            let mut camera = viewport.camera.write().unwrap();
            if camera.target == glam::Vec3::ZERO {
                let locked_protein_store = viewport.store.read().unwrap();
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
    }

    fn render(&mut self, queue: &wgpu::Queue) -> Result<(), wgpu::SurfaceError> {
        let surface_texture_output = self.surface.get_current_texture()?;
        let texture_view_for_rendering = surface_texture_output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        // Prepare data for all viewports in a grid
        let mut viewport_render_data_collection = Vec::new();
        let (column_count, row_count) = self.calculate_grid_dimensions();
        
        let viewport_width = self.config.width as f32 / column_count as f32;
        let viewport_height = self.config.height as f32 / row_count as f32;

        for (viewport_index, viewport_state) in self.viewports.iter().enumerate() {
            let column_index = (viewport_index as u32 % column_count) as f32;
            let row_index = (viewport_index as u32 / column_count) as f32;
            
            let x_coordinate_offset = column_index * viewport_width;
            let y_coordinate_offset = row_index * viewport_height;
            
            viewport_render_data_collection.push((
                viewport_state.store.clone(),
                viewport_state.camera.clone(),
                &viewport_state.selected_atoms[..],
                &viewport_state.measurements[..],
                [x_coordinate_offset, y_coordinate_offset, viewport_width, viewport_height],
            ));
        }

        let rendered_graphics_command_buffer = self
            .renderer
            .render_viewports(&texture_view_for_rendering, &viewport_render_data_collection);
        
        queue.submit(std::iter::once(rendered_graphics_command_buffer));
        surface_texture_output.present();

        Ok(())
    }

    fn get_viewport_at_mouse(&self, mouse_x_coordinate: f64, mouse_y_coordinate: f64) -> Option<(usize, f32, f32, f32)> {
        let (column_count, row_count) = self.calculate_grid_dimensions();
        if column_count == 0 || row_count == 0 { return None; }
        
        let viewport_width = self.config.width as f64 / column_count as f64;
        let viewport_height = self.config.height as f64 / row_count as f64;
        
        let column_index = (mouse_x_coordinate / viewport_width).floor() as u32;
        let row_index = (mouse_y_coordinate / viewport_height).floor() as u32;
        
        let target_viewport_index = (row_index * column_count + column_index) as usize;
        
        if target_viewport_index < self.viewports.len() {
            let local_x_coordinate = mouse_x_coordinate - (column_index as f64 * viewport_width);
            let local_y_coordinate = mouse_y_coordinate - (row_index as f64 * viewport_height);
            Some((target_viewport_index, local_x_coordinate as f32, local_y_coordinate as f32, viewport_width as f32))
        } else {
            None
        }
    }

    fn handle_mouse_input(
        &mut self,
        element_press_state: ElementState,
        mouse_button_identifier: MouseButton,
    ) {
        match mouse_button_identifier {
            MouseButton::Left => {
                let was_previously_pressed = self.mouse_pressed;
                self.mouse_pressed = element_press_state == ElementState::Pressed;

                if self.mouse_pressed && !was_previously_pressed {
                    if let Some((current_mouse_x_coordinate, current_mouse_y_coordinate)) =
                        self.last_recorded_mouse_cursor_position
                    {
                        let (column_count, row_count) = self.calculate_grid_dimensions();
                        if let Some((viewport_index, local_x_coordinate, local_y_coordinate, viewport_width)) = self.get_viewport_at_mouse(current_mouse_x_coordinate, current_mouse_y_coordinate) {
                            // Capture the viewport where the interaction started
                            self.captured_viewport_index = Some(viewport_index);
                            
                            let target_viewport = &mut self.viewports[viewport_index];
                            let (ray_origin_point, ray_direction_unit_vector) = target_viewport
                                .camera
                                .read()
                                .unwrap()
                                .calculate_ray_from_screen_coordinates(
                                    glam::Vec2::new(
                                        local_x_coordinate,
                                        local_y_coordinate,
                                    ),
                                    glam::Vec2::new(
                                        viewport_width,
                                        self.config.height as f32 / row_count as f32,
                                    ),
                                );

                            let mut closest_intersection_hit_data: Option<(String, usize, f32)> = None;

                            let locked_protein_store = target_viewport.store.read().unwrap();
                            for protein_shared_handle in locked_protein_store.iter() {
                                let protein_locked_data = protein_shared_handle.read().unwrap();
                                let protein_identifier_name = protein_locked_data.name.clone();

                                let current_representation_mode = protein_locked_data.representation;

                                let base_sphere_hitbox_radius = match current_representation_mode {
                                    Representation::Spheres | Representation::BackboneAndSpheres => 1.5,
                                    Representation::BallAndStick => 0.4,
                                    Representation::SpaceFilling => 1.7,
                                    _ => 1.0,
                                };

                                for (atom_global_index, current_atom_hierarchy) in protein_locked_data
                                    .pdb
                                    .atoms_with_hierarchy()
                                    .into_iter()
                                    .enumerate()
                                {
                                    let current_atom_reference = current_atom_hierarchy.atom();

                                    if matches!(
                                        current_representation_mode,
                                        Representation::Spheres | Representation::BackboneAndSpheres
                                    ) && current_atom_reference.name() != "CA"
                                    {
                                        continue;
                                    }

                                    let atom_position_tuple = current_atom_reference.pos();
                                    let atom_world_position_vector = glam::Vec3::new(
                                        atom_position_tuple.0 as f32,
                                        atom_position_tuple.1 as f32,
                                        atom_position_tuple.2 as f32,
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
                                let target_atom_selection_identifier =
                                    (hit_protein_name, hit_atom_index);

                                if target_viewport
                                    .selected_atoms
                                    .contains(&target_atom_selection_identifier)
                                {
                                    target_viewport.selected_atoms
                                        .retain(|item| item != &target_atom_selection_identifier);
                                    target_viewport.distance_has_been_printed_for_current_selection = false;
                                } else {
                                    if target_viewport.selected_atoms.len() >= 2 {
                                        target_viewport.selected_atoms.clear();
                                        target_viewport.measurements.clear();
                                    }
                                    target_viewport.selected_atoms.push(target_atom_selection_identifier);
                                    target_viewport.distance_has_been_printed_for_current_selection = false;
                                }
                            }
                        }
                    }
                } else if !self.mouse_pressed {
                    // Release capture when mouse is released
                    self.captured_viewport_index = None;
                }
            }
            MouseButton::Right => {
                self.right_mouse_button_is_pressed = element_press_state == ElementState::Pressed;
                if !self.right_mouse_button_is_pressed {
                    self.captured_viewport_index = None;
                } else if let Some((current_mouse_x_coordinate, current_mouse_y_coordinate)) = self.last_recorded_mouse_cursor_position {
                    // Start capture for right-click drag if not already captured
                    if self.captured_viewport_index.is_none() {
                        if let Some((viewport_index, _, _, _)) = self.get_viewport_at_mouse(current_mouse_x_coordinate, current_mouse_y_coordinate) {
                            self.captured_viewport_index = Some(viewport_index);
                        }
                    }
                }
            }
            _ => {}
        }
    }

    fn handle_mouse_move(&mut self, cursor_screen_position: (f64, f64)) {
        if let Some((previous_mouse_x_coordinate, previous_mouse_y_coordinate)) = self.last_recorded_mouse_cursor_position
        {
            let mouse_movement_delta_x = (cursor_screen_position.0 - previous_mouse_x_coordinate) as f32;
            let mouse_movement_delta_y = (cursor_screen_position.1 - previous_mouse_y_coordinate) as f32;

            // Prioritize the captured viewport for rotation and translation
            let viewport_index_to_interact_with = self.captured_viewport_index.or_else(|| {
                self.get_viewport_at_mouse(cursor_screen_position.0, cursor_screen_position.1).map(|(index, _, _, _)| index)
            });

            if let Some(active_viewport_index) = viewport_index_to_interact_with {
                if active_viewport_index < self.viewports.len() {
                    let target_viewport = &self.viewports[active_viewport_index];
                    if self.mouse_pressed {
                        target_viewport.camera.write().unwrap().rotate(
                            -mouse_movement_delta_x * 0.01,
                            mouse_movement_delta_y * 0.01,
                        );
                    } else if self.right_mouse_button_is_pressed {
                        target_viewport.camera
                            .write()
                            .unwrap()
                            .translate_camera_target_position(
                                mouse_movement_delta_x,
                                mouse_movement_delta_y,
                            );
                    }
                }
            }
        }
        self.last_recorded_mouse_cursor_position = Some(cursor_screen_position);
    }

    fn handle_scroll(&mut self, scroll_delta_value: f32) {
        if let Some((current_mouse_x, current_mouse_y)) = self.last_recorded_mouse_cursor_position {
            if let Some((viewport_index, _, _, _)) = self.get_viewport_at_mouse(current_mouse_x, current_mouse_y) {
                self.viewports[viewport_index].camera.write().unwrap().zoom(scroll_delta_value);
            }
        }
    }

    fn handle_key(&mut self, pressed_physical_keycode: KeyCode) {
        // Apply key events to ALL viewports for now, or we could track active viewport
        for viewport in &mut self.viewports {
            match pressed_physical_keycode {
                KeyCode::KeyR => {
                    viewport.selected_atoms.clear();
                    viewport.measurements.clear();
                    viewport.distance_has_been_printed_for_current_selection = false;
                    let mut camera = viewport.camera.write().unwrap();
                    *camera = Camera::new(camera.aspect);
                    let locked_protein_store = viewport.store.read().unwrap();
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

                KeyCode::Digit1 => set_viewport_representation(viewport, Representation::Spheres),
                KeyCode::Digit2 => set_viewport_representation(viewport, Representation::Backbone),
                KeyCode::Digit3 => set_viewport_representation(viewport, Representation::BackboneAndSpheres),
                KeyCode::Digit4 => set_viewport_representation(viewport, Representation::Sticks),
                KeyCode::Digit5 => set_viewport_representation(viewport, Representation::BallAndStick),
                KeyCode::Digit6 => set_viewport_representation(viewport, Representation::SpaceFilling),
                KeyCode::Digit7 => set_viewport_representation(viewport, Representation::Lines),

                KeyCode::KeyC => set_viewport_color_scheme(viewport, ColorScheme::ByChain),
                KeyCode::KeyB => set_viewport_color_scheme(viewport, ColorScheme::ByBFactor),
                KeyCode::KeyE => set_viewport_color_scheme(viewport, ColorScheme::ByElement),
                KeyCode::KeyS => set_viewport_color_scheme(viewport, ColorScheme::BySecondary),

                KeyCode::KeyM => {
                    if viewport.selected_atoms.len() == 2
                        && !viewport.distance_has_been_printed_for_current_selection
                    {
                        let locked_protein_store = viewport.store.read().unwrap();
                        let mut global_atom_positions_lookup_table = std::collections::HashMap::new();
                        for protein_shared_handle in locked_protein_store.iter() {
                            let protein_locked_data = protein_shared_handle.read().unwrap();
                            for (atom_indexing_counter, current_atom_hierarchy) in protein_locked_data
                                .pdb
                                .atoms_with_hierarchy()
                                .into_iter()
                                .enumerate()
                            {
                                let atom_position_tuple = current_atom_hierarchy.atom().pos();
                                global_atom_positions_lookup_table.insert(
                                    (protein_locked_data.name.clone(), atom_indexing_counter),
                                    glam::Vec3::new(
                                        atom_position_tuple.0 as f32,
                                        atom_position_tuple.1 as f32,
                                        atom_position_tuple.2 as f32,
                                    ),
                                );
                            }
                        }

                        if let (
                            Some(first_atom_selection_identifier),
                            Some(second_atom_selection_identifier),
                        ) = (viewport.selected_atoms.get(0), viewport.selected_atoms.get(1))
                        {
                            if let (Some(first_atom_world_position), Some(second_atom_world_position)) = (
                                global_atom_positions_lookup_table.get(first_atom_selection_identifier),
                                global_atom_positions_lookup_table
                                    .get(second_atom_selection_identifier),
                            ) {
                                let calculated_distance_value =
                                    first_atom_world_position.distance(*second_atom_world_position);
                                info!(
                                    "Calculated distance: {:.2} Angstroms",
                                    calculated_distance_value
                                );

                                viewport.measurements.clear();
                                viewport.measurements.push((0, 1));
                                viewport.distance_has_been_printed_for_current_selection = true;
                            }
                        }
                    }
                }
                KeyCode::KeyF => {
                    let locked_protein_store = viewport.store.read().unwrap();
                    for protein_shared_handle in locked_protein_store.iter() {
                        let mut protein_mutable_data = protein_shared_handle.write().unwrap();
                        let current_visibility_status = protein_mutable_data
                            .molecular_surface_mesh
                            .is_surface_visible;
                        protein_mutable_data
                            .molecular_surface_mesh
                            .is_surface_visible = !current_visibility_status;
                    }
                }
                _ => {}
            }
        }
        
        if pressed_physical_keycode == KeyCode::Escape {
            std::process::exit(0);
        }
    }
}

fn set_viewport_representation(
    viewport: &ViewportState,
    target_representation_mode: Representation,
) {
    let locked_protein_store = viewport.store.read().unwrap();
    for protein_shared_handle in locked_protein_store.iter() {
        let mut protein_mutable_data = protein_shared_handle.write().unwrap();
        protein_mutable_data.representation = target_representation_mode;
    }
}

fn set_viewport_color_scheme(viewport: &ViewportState, target_color_scheme: ColorScheme) {
    let locked_protein_store = viewport.store.read().unwrap();
    for protein_shared_handle in locked_protein_store.iter() {
        let mut protein_mutable_data = protein_shared_handle.write().unwrap();
        protein_mutable_data.color_scheme = target_color_scheme;
    }
}

/// The main application state managing a single window with viewports
struct App {
    instance: wgpu::Instance,
    adapter: wgpu::Adapter,
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    window_state: Option<WindowState>,
    script_engine: ScriptEngine,
    global_store: Arc<RwLock<ProteinStore>>,
    global_camera: Arc<RwLock<Camera>>,
    script_rx: crossbeam_channel::Receiver<PathBuf>,
    script_event_rx: crossbeam_channel::Receiver<ScriptEvent>,
    _watcher: notify::RecommendedWatcher,
    
    /// The path of the currently executing script
    active_script_path_identifier: String,
}

impl App {
    async fn new(initial_script_path: Option<String>) -> Self {
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

        let lua_script_engine =
            ScriptEngine::new(global_store.clone(), global_camera.clone(), script_event_tx)
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
            window_state: None,
            script_engine: lua_script_engine,
            global_store,
            global_camera,
            script_rx: script_path_receiver,
            script_event_rx,
            _watcher: file_system_watcher,
            active_script_path_identifier: initial_script_path.unwrap_or_else(|| "none".to_string()),
        }
    }

    fn ensure_window(&mut self, event_loop: &winit::event_loop::EventLoopWindowTarget<()>) -> &mut WindowState {
        if self.window_state.is_none() {
            let window = Arc::new(
                WindowBuilder::new()
                    .with_title("Protein Viewer (Split View)")
                    .with_inner_size(winit::dpi::LogicalSize::new(1280, 800))
                    .build(event_loop)
                    .unwrap(),
            );
            self.window_state = Some(WindowState::new(
                window,
                &self.instance,
                &self.adapter,
                self.device.clone(),
                self.queue.clone(),
                self.active_script_path_identifier.clone(),
            ));
        }
        self.window_state.as_mut().unwrap()
    }

    fn handle_script_events(&mut self, event_loop: &winit::event_loop::EventLoopWindowTarget<()>) {
        while let Ok(event) = self.script_event_rx.try_recv() {
            match event {
                ScriptEvent::NewWindow(store, camera) => {
                    let window_state = self.ensure_window(event_loop);
                    window_state.add_viewport(store, camera);
                }
            }
        }
    }

    fn update(&mut self, event_loop: &winit::event_loop::EventLoopWindowTarget<()>) {
        let mut last_reloaded_path = None;
        while let Ok(reloaded_script_path) = self.script_rx.try_recv() {
            last_reloaded_path = Some(reloaded_script_path);
        }

        if let Some(reloaded_script_path) = last_reloaded_path {
            info!("Reloading script from path: {:?}", reloaded_script_path);
            if let Some(reloaded_script_path_string) = reloaded_script_path.to_str() {
                self.active_script_path_identifier = reloaded_script_path_string.to_string();
                
                while self.script_event_rx.try_recv().is_ok() {}

                if let Some(current_window_state) = &mut self.window_state {
                    current_window_state.viewports.clear();
                    current_window_state.active_script_path_identifier = self.active_script_path_identifier.clone();
                    current_window_state.update_window_title_with_current_state();
                }
                self.global_store.write().unwrap().clear();

                if let Err(hot_reload_error) =
                    self.script_engine.run_file(reloaded_script_path_string)
                {
                    error!("Standardized script error log: {}", hot_reload_error);
                }
            }
        }

        self.handle_script_events(event_loop);

        if let Some(current_window_state) = &mut self.window_state {
            current_window_state.update();
        }
    }
}

async fn run(initial_script: Option<String>) {
    let main_event_loop = EventLoop::new().unwrap();
    let mut app = App::new(initial_script.clone()).await;

    if let Some(provided_initial_script_path) = initial_script {
        if let Err(script_execution_error) =
            app.script_engine.run_file(&provided_initial_script_path)
        {
            error!(
                "Error running explicitly provided script {}: {}",
                provided_initial_script_path, script_execution_error
            );
        }
    }

    main_event_loop.set_control_flow(ControlFlow::Poll);
    let _ = main_event_loop.run(move |winit_event, event_loop_window_target| {
        match winit_event {
            Event::WindowEvent { window_id: _, event } => {
                if let Some(window_state) = app.window_state.as_mut() {
                    match event {
                        WindowEvent::CloseRequested => {
                            event_loop_window_target.exit();
                        }
                        WindowEvent::Resized(new_window_size) => {
                            window_state.resize(&app.device, new_window_size)
                        }
                        WindowEvent::MouseInput { state, button, .. } => {
                            window_state.handle_mouse_input(state, button)
                        }
                        WindowEvent::CursorMoved { position, .. } => {
                            window_state.handle_mouse_move((position.x, position.y))
                        }
                        WindowEvent::MouseWheel { delta, .. } => {
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
                        WindowEvent::RedrawRequested => match window_state.render(&app.queue) {
                            Ok(_) => {}
                            Err(wgpu::SurfaceError::Lost) => {
                                window_state.resize(&app.device, window_state.window.inner_size())
                            }
                            Err(wgpu::SurfaceError::OutOfMemory) => {
                                event_loop_window_target.exit();
                            }
                            Err(e) => error!("Standardized render error log: {:?}", e),
                        },
                        _ => {}
                    }
                }
            }
            Event::AboutToWait => {
                app.update(event_loop_window_target);
                if let Some(ws) = &app.window_state {
                    ws.window.request_redraw();
                }
            }
            _ => {}
        }
    });
}

fn main() {
    let application_arguments = Args::parse();

    // Initialize the standardized logging system with tracing
    tracing_subscriber::fmt()
        .with_target(false)
        .with_env_filter(
            tracing_subscriber::EnvFilter::builder()
                .with_default_directive(tracing_subscriber::filter::LevelFilter::INFO.into())
                .from_env_lossy()
                .add_directive("wgpu_core=warn".parse().unwrap())
                .add_directive("wgpu_hal=warn".parse().unwrap())
                .add_directive("naga=warn".parse().unwrap()),
        )
        .init();

    pollster::block_on(run(application_arguments.script));
}
