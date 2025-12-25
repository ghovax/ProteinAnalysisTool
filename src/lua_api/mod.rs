//! Lua scripting engine for protein manipulation
//!
//! This module provides the `ScriptEngine` which embeds a Lua interpreter
//! and exposes an API for fetching, loading, and inspecting protein structures

mod protein;

use mlua::{Lua, Result};
use std::sync::{Arc, RwLock};

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

        // pdb.fetch(code) fetches a protein from RCSB by its PDB identifier
        let cloned_event_sender = event_sender.clone();
        let cloned_global_store = protein_data_store.clone();
        pdb_api_table.set(
            "fetch",
            lua_runtime_instance.create_function(move |_lua, requested_pdb_code: String| {
                println!("Fetching PDB: {}...", requested_pdb_code);
                let new_store = Arc::new(RwLock::new(ProteinStore::new()));
                let mut locked_new_store = new_store.write().unwrap();
                match locked_new_store.fetch(&requested_pdb_code) {
                    Ok(shared_protein_handle) => {
                        // Also add to global store for listing/saving
                        cloned_global_store.write().unwrap().add(shared_protein_handle.clone());

                        let locked_protein_data = shared_protein_handle.read().unwrap();
                        println!(
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

        // pdb.load(path) loads a protein from a local file (PDB or mmCIF)
        let cloned_event_sender = event_sender.clone();
        let cloned_global_store = protein_data_store.clone();
        pdb_api_table.set(
            "load",
            lua_runtime_instance.create_function(move |_lua, requested_file_path: String| {
                println!("Loading file: {}...", requested_file_path);
                let new_store = Arc::new(RwLock::new(ProteinStore::new()));
                let mut locked_new_store = new_store.write().unwrap();
                match locked_new_store.load(&requested_file_path) {
                    Ok(shared_protein_handle) => {
                        // Also add to global store for listing/saving
                        cloned_global_store.write().unwrap().add(shared_protein_handle.clone());

                        let locked_protein_data = shared_protein_handle.read().unwrap();
                        println!(
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

        // pdb.list() returns a table of names for all currently loaded protein identifiers
        let cloned_protein_store_reference = protein_data_store.clone();
        pdb_api_table.set(
            "list",
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
            "get_pos",
            lua_runtime_instance.create_function(move |_, ()| {
                let camera = cloned_camera.read().unwrap();
                let pos = camera.position();
                Ok((pos.x, pos.y, pos.z))
            })?,
        )?;

        let cloned_camera = camera.clone();
        camera_api_table.set(
            "set_target",
            lua_runtime_instance.create_function(move |_, (x, y, z): (f32, f32, f32)| {
                let mut camera = cloned_camera.write().unwrap();
                camera.target = glam::Vec3::new(x, y, z);
                Ok(())
            })?,
        )?;

        let cloned_camera = camera.clone();
        camera_api_table.set(
            "set_params",
            lua_runtime_instance.create_function(
                move |_, (dist, yaw, pitch): (f32, f32, f32)| {
                    let mut camera = cloned_camera.write().unwrap();
                    camera.distance = dist;
                    camera.yaw = yaw;
                    camera.pitch = pitch;
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
            "save",
            lua_runtime_instance.create_function(move |_, path: String| {
                use std::io::Write;
                let mut file = std::fs::File::create(&path)
                    .map_err(|e| mlua::Error::RuntimeError(e.to_string()))?;

                writeln!(file, "-- Protein Viewer Session State").unwrap();

                // Save proteins
                let store = cloned_store.read().unwrap();
                for (i, protein) in store.iter().enumerate() {
                    let p = protein.read().unwrap();
                    let var_name = format!("p{}", i);
                    match &p.source {
                        crate::protein::structure::ProteinSource::Rcsb(code) => {
                            writeln!(file, "local {} = pdb.fetch(\"{}\")", var_name, code).unwrap();
                        }
                        crate::protein::structure::ProteinSource::File(path) => {
                            writeln!(file, "local {} = pdb.load(\"{}\")", var_name, path).unwrap();
                        }
                    }

                    let repr = match p.representation {
                        crate::protein::structure::Representation::Spheres => "spheres",
                        crate::protein::structure::Representation::Backbone => "backbone",
                        crate::protein::structure::Representation::BackboneAndSpheres => "both",
                        crate::protein::structure::Representation::Sticks => "sticks",
                        crate::protein::structure::Representation::BallAndStick => "ball-and-stick",
                        crate::protein::structure::Representation::SpaceFilling => "space-filling",
                        crate::protein::structure::Representation::Lines => "lines",
                    };
                    writeln!(file, "{}:representation(\"{}\")", var_name, repr).unwrap();

                    match p.color_scheme {
                        crate::protein::structure::ColorScheme::ByChain => {
                            writeln!(file, "{}:color_by(\"chain\")", var_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::ByElement => {
                            writeln!(file, "{}:color_by(\"element\")", var_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::ByBFactor => {
                            writeln!(file, "{}:color_by(\"bfactor\")", var_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::BySecondary => {
                            writeln!(file, "{}:color_by(\"secondary\")", var_name).unwrap()
                        }
                        crate::protein::structure::ColorScheme::Uniform(c) => {
                            writeln!(file, "{}:color({}, {}, {})", var_name, c[0], c[1], c[2])
                                .unwrap()
                        }
                    }

                    if !p.visible {
                        writeln!(file, "{}:hide()", var_name).unwrap();
                    }
                }

                // Save camera
                let cam = cloned_camera.read().unwrap();
                writeln!(
                    file,
                    "camera.set_target({}, {}, {})",
                    cam.target.x, cam.target.y, cam.target.z
                )
                .unwrap();
                writeln!(
                    file,
                    "camera.set_params({}, {}, {})",
                    cam.distance, cam.yaw, cam.pitch
                )
                .unwrap();

                Ok(())
            })?,
        )?;

        session_api_table.set(
            "load",
            lua_runtime_instance.create_function(move |lua, path: String| {
                let code = std::fs::read_to_string(&path)
                    .map_err(|e| mlua::Error::RuntimeError(e.to_string()))?;
                lua.load(&code).exec()
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
                    println!("[Lua] {}", formatted_output_strings.join("\t"));
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
        self.lua.load(lua_code_string).exec()
    }

    /// Reads a file and executes its content as Lua code
    pub fn run_file(&self, lua_file_path: &str) -> Result<()> {
        match std::fs::read_to_string(lua_file_path) {
            Ok(loaded_lua_code_string) => self.run_script(&loaded_lua_code_string),
            Err(file_read_error_message) => Err(mlua::Error::RuntimeError(format!(
                "Failed to read {}: {}",
                lua_file_path, file_read_error_message
            ))),
        }
    }
}
