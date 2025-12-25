//! Lua scripting engine for protein manipulation
//!
//! This module provides the `ScriptEngine` which embeds a Lua interpreter
//! and exposes an API for fetching, loading, and inspecting protein structures

mod protein;
mod selection;

use mlua::{Lua, Result};
use std::sync::{Arc, RwLock};
use tracing::{info, info_span};

use crate::protein::ProteinStore;

use crate::renderer::Camera;

pub use protein::LuaProtein;

/// Events sent from the script engine to the main application
pub enum ScriptEvent {
    /// Request to open a new window with a specific protein and camera
    NewWindow(Arc<RwLock<ProteinStore>>, Arc<RwLock<Camera>>),
}

/// The Lua script engine
pub struct ScriptEngine {
    lua: Lua,
}

impl ScriptEngine {
    /// Creates a new `ScriptEngine` and initializes the `pdb` global API table
    pub fn new(
        protein_data_store: Arc<RwLock<ProteinStore>>,
        camera: Arc<RwLock<Camera>>,
        event_sender: crossbeam_channel::Sender<ScriptEvent>,
    ) -> Result<Self> {
        let lua_runtime_instance = Lua::new();

        // Create pdb module
        let pdb_api_table = lua_runtime_instance.create_table()?;

        // pdb.fetch_protein_from_rcsb(code) fetches a protein from RCSB by its PDB identifier
        let cloned_event_sender = event_sender.clone();
        let cloned_global_store = protein_data_store.clone();
        pdb_api_table.set(
            "fetch_protein_from_rcsb",
            lua_runtime_instance.create_function(move |_lua, requested_pdb_code: String| {
                info!("Fetching PDB: {}...", requested_pdb_code);
                let new_store = Arc::new(RwLock::new(ProteinStore::new()));
                let mut locked_new_store = new_store.write().unwrap();
                match locked_new_store.fetch(&requested_pdb_code) {
                    Ok(shared_protein_handle) => {
                        // Also add to global store for listing/saving
                        cloned_global_store.write().unwrap().add(shared_protein_handle.clone());

                        let locked_protein_data = shared_protein_handle.read().unwrap();
                        info!(
                            "Loaded {} with {} atoms, {} chains",
                            locked_protein_data.name,
                            locked_protein_data.atom_count(),
                            locked_protein_data.chain_ids().len()
                        );
                        drop(locked_protein_data);

                        // Create a new camera for this window
                        let new_camera = Arc::new(RwLock::new(Camera::new(1.0))); // Aspect will be set on resize

                        // Signal to open a new window
                        let _ = cloned_event_sender.send(ScriptEvent::NewWindow(new_store.clone(), new_camera));

                        Ok(LuaProtein::new(shared_protein_handle.clone()))
                    }
                    Err(lua_runtime_error_message) => {
                        Err(mlua::Error::RuntimeError(lua_runtime_error_message))
                    }
                }
            })?,
        )?;

        // pdb.load_protein_from_local_file(path) loads a protein from a local file (PDB or mmCIF)
        let cloned_event_sender = event_sender.clone();
        let cloned_global_store = protein_data_store.clone();
        pdb_api_table.set(
            "load_protein_from_local_file",
            lua_runtime_instance.create_function(move |_lua, requested_file_path: String| {
                info!("Loading file: {}...", requested_file_path);
                let new_store = Arc::new(RwLock::new(ProteinStore::new()));
                let mut locked_new_store = new_store.write().unwrap();
                match locked_new_store.load(&requested_file_path) {
                    Ok(shared_protein_handle) => {
                        // Also add to global store for listing/saving
                        cloned_global_store.write().unwrap().add(shared_protein_handle.clone());

                        let locked_protein_data = shared_protein_handle.read().unwrap();
                        info!(
                            "Loaded {} with {} atoms, {} chains",
                            locked_protein_data.name,
                            locked_protein_data.atom_count(),
                            locked_protein_data.chain_ids().len()
                        );
                        drop(locked_protein_data);

                        // Create a new camera for this window
                        let new_camera = Arc::new(RwLock::new(Camera::new(1.0))); // Aspect will be set on resize

                        // Signal to open a new window
                        let _ = cloned_event_sender.send(ScriptEvent::NewWindow(new_store.clone(), new_camera));

                        Ok(LuaProtein::new(shared_protein_handle.clone()))
                    }
                    Err(lua_runtime_error_message) => {
                        Err(mlua::Error::RuntimeError(lua_runtime_error_message))
                    }
                }
            })?,
        )?;

        // pdb.get_loaded_protein_identifiers() returns a table of names for all currently loaded protein identifiers
        let cloned_protein_store_reference = protein_data_store.clone();
        pdb_api_table.set(
            "get_loaded_protein_identifiers",
            lua_runtime_instance.create_function(move |lua_context, ()| {
                let locked_protein_store = cloned_protein_store_reference.read().unwrap();
                let available_protein_names = locked_protein_store.list();
                lua_context.create_sequence_from(available_protein_names)
            })?,
        )?;

        lua_runtime_instance.globals().set("pdb", pdb_api_table)?;

        // Create camera module
        let camera_api_table = lua_runtime_instance.create_table()?;
        let cloned_camera = camera.clone();
        camera_api_table.set(
            "get_camera_world_position",
            lua_runtime_instance.create_function(move |_, ()| {
                let locked_camera_handle = cloned_camera.read().unwrap();
                let camera_world_position = locked_camera_handle.position();
                Ok((camera_world_position.x, camera_world_position.y, camera_world_position.z))
            })?,
        )?;

        let cloned_camera = camera.clone();
        camera_api_table.set(
            "set_camera_focus_target",
            lua_runtime_instance.create_function(move |_, (target_coordinate_x, target_coordinate_y, target_coordinate_z): (f32, f32, f32)| {
                let mut mutable_camera_handle = cloned_camera.write().unwrap();
                mutable_camera_handle.target = glam::Vec3::new(target_coordinate_x, target_coordinate_y, target_coordinate_z);
                Ok(())
            })?,
        )?;

        let cloned_camera = camera.clone();
        camera_api_table.set(
            "set_camera_spherical_parameters",
            lua_runtime_instance.create_function(
                move |_, (requested_camera_distance, requested_camera_yaw, requested_camera_pitch): (f32, f32, f32)| {
                    let mut mutable_camera_handle = cloned_camera.write().unwrap();
                    mutable_camera_handle.distance = requested_camera_distance;
                    mutable_camera_handle.yaw = requested_camera_yaw;
                    mutable_camera_handle.pitch = requested_camera_pitch;
                    Ok(())
                },
            )?,
        )?;

        lua_runtime_instance
            .globals()
            .set("camera", camera_api_table)?;

        // Create session module
        let session_api_table = lua_runtime_instance.create_table()?;
        let cloned_store = protein_data_store.clone();
        let cloned_camera = camera.clone();
        session_api_table.set(
            "save_application_session",
            lua_runtime_instance.create_function(move |_, output_file_path: String| {
                use std::io::Write;
                let mut session_output_file = std::fs::File::create(&output_file_path)
                    .map_err(|file_creation_error| mlua::Error::RuntimeError(file_creation_error.to_string()))?;

                writeln!(session_output_file, "-- Protein Viewer Session State").unwrap();

                // Save proteins
                let locked_protein_store = cloned_store.read().unwrap();
                for (protein_index, shared_protein_handle) in locked_protein_store.iter().enumerate() {
                    let locked_protein_data = shared_protein_handle.read().unwrap();
                    let lua_protein_variable_name = format!("protein_object_{}", protein_index);
                    match &locked_protein_data.source {
                        crate::protein::structure::ProteinSource::Rcsb(pdb_code) => {
                            writeln!(session_output_file, "local {} = pdb.fetch_protein_from_rcsb(\"{}\")", lua_protein_variable_name, pdb_code).unwrap();
                        }
                        crate::protein::structure::ProteinSource::File(absolute_file_path) => {
                            writeln!(session_output_file, "local {} = pdb.load_protein_from_local_file(\"{}\")", lua_protein_variable_name, absolute_file_path).unwrap();
                        }
                    }

                    let representation_mode_string = match locked_protein_data.representation {
                        crate::protein::structure::Representation::Spheres => "spheres",
                        crate::protein::structure::Representation::Backbone => "backbone",
                        crate::protein::structure::Representation::BackboneAndSpheres => "both",
                        crate::protein::structure::Representation::Sticks => "sticks",
                        crate::protein::structure::Representation::BallAndStick => "ball-and-stick",
                        crate::protein::structure::Representation::SpaceFilling => "space-filling",
                        crate::protein::structure::Representation::Lines => "lines",
                    };
                    writeln!(session_output_file, "{}:set_representation_mode(\"{}\")", lua_protein_variable_name, representation_mode_string).unwrap();

                    match locked_protein_data.color_scheme {
                        crate::protein::structure::ColorScheme::ByChain => {
                            writeln!(session_output_file, "{}:set_color_scheme_by_property(\"chain\")", lua_protein_variable_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::ByElement => {
                            writeln!(session_output_file, "{}:set_color_scheme_by_property(\"element\")", lua_protein_variable_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::ByBFactor => {
                            writeln!(session_output_file, "{}:set_color_scheme_by_property(\"bfactor\")", lua_protein_variable_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::BySecondary => {
                            writeln!(session_output_file, "{}:set_color_scheme_by_property(\"secondary\")", lua_protein_variable_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::Uniform(rgb_color_components) => {
                            writeln!(session_output_file, "{}:set_uniform_rgb_color({}, {}, {})", lua_protein_variable_name, rgb_color_components[0], rgb_color_components[1], rgb_color_components[2])
                                .unwrap()
                        }
                    }

                    if !locked_protein_data.visible {
                        writeln!(session_output_file, "{}:set_visibility_off()", lua_protein_variable_name).unwrap();
                    }
                }

                // Save camera state
                let locked_camera_handle = cloned_camera.read().unwrap();
                writeln!(
                    session_output_file,
                    "camera.set_camera_focus_target({}, {}, {})",
                    locked_camera_handle.target.x, locked_camera_handle.target.y, locked_camera_handle.target.z
                )
                .unwrap();
                writeln!(
                    session_output_file,
                    "camera.set_camera_spherical_parameters({}, {}, {})",
                    locked_camera_handle.distance, locked_camera_handle.yaw, locked_camera_handle.pitch
                )
                .unwrap();

                Ok(())
            })?,
        )?;

        session_api_table.set(
            "load_application_session",
            lua_runtime_instance.create_function(move |lua_runtime_context, session_file_path: String| {
                let session_script_content = std::fs::read_to_string(&session_file_path)
                    .map_err(|file_read_error| mlua::Error::RuntimeError(file_read_error.to_string()))?;
                lua_runtime_context.load(&session_script_content).exec()
            })?,
        )?;

        lua_runtime_instance
            .globals()
            .set("session", session_api_table)?;

        // Simple print override for cleaner output
        let _ = lua_runtime_instance.globals().set(
            "print",
            lua_runtime_instance.create_function(
                |_, variadic_print_arguments: mlua::Variadic<mlua::Value>| {
                    let formatted_output_strings: Vec<String> = variadic_print_arguments
                        .iter()
                        .map(|argument_value| match argument_value {
                            mlua::Value::Nil => "nil".to_string(),
                            mlua::Value::Boolean(boolean_value) => boolean_value.to_string(),
                            mlua::Value::Integer(integer_value) => integer_value.to_string(),
                            mlua::Value::Number(numeric_value) => format!("{:.4}", numeric_value),
                            mlua::Value::String(string_value) => {
                                string_value.to_str().unwrap_or("").to_string()
                            }
                            _ => format!("{:?}", argument_value),
                        })
                        .collect();
                    info!("{}", formatted_output_strings.join("\t"));
                    Ok(())
                },
            )?,
        );

        Ok(Self {
            lua: lua_runtime_instance,
        })
    }

    /// Executes a string as Lua code
    pub fn run_script(&self, lua_code_string: &str) -> Result<()> {
        let span = info_span!("lua");
        let _enter = span.enter();
        self.lua.load(lua_code_string).exec()
    }

    /// Reads a file and executes its content as Lua code
    pub fn run_file(&self, lua_file_path: &str) -> Result<()> {
        let span = info_span!("lua", script = lua_file_path);
        let _enter = span.enter();
        match std::fs::read_to_string(lua_file_path) {
            Ok(loaded_lua_code_string) => self.run_script(&loaded_lua_code_string),
            Err(file_read_error_message) => Err(mlua::Error::RuntimeError(format!(
                "Failed to read {}: {}",
                lua_file_path, file_read_error_message
            ))),
        }
    }
}
