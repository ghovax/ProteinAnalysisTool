//! Lua scripting engine for protein manipulation
//!
//! This module provides the `ScriptEngine` which embeds a Lua interpreter
//! and exposes an API for fetching, loading, and inspecting protein structures

mod protein;

use mlua::{Lua, Result};
use std::sync::{Arc, RwLock};

use crate::protein::ProteinStore;

pub use protein::LuaProtein;

/// The Lua script engine
pub struct ScriptEngine {
    lua: Lua,
}

impl ScriptEngine {
    /// Creates a new `ScriptEngine` and initializes the `pdb` global API table
    pub fn new(protein_data_store: Arc<RwLock<ProteinStore>>) -> Result<Self> {
        let lua_runtime_instance = Lua::new();

        // Create pdb module
        let pdb_api_table = lua_runtime_instance.create_table()?;

        // pdb.fetch(code) fetches a protein from RCSB by its PDB identifier
        let cloned_protein_store_reference = protein_data_store.clone();
        pdb_api_table.set(
            "fetch",
            lua_runtime_instance.create_function(move |_lua, requested_pdb_code: String| {
                println!("Fetching PDB: {}...", requested_pdb_code);
                let mut locked_protein_store = cloned_protein_store_reference.write().unwrap();
                match locked_protein_store.fetch(&requested_pdb_code) {
                    Ok(shared_protein_handle) => {
                        let locked_protein_data = shared_protein_handle.read().unwrap();
                        println!(
                            "Loaded {} with {} atoms, {} chains",
                            locked_protein_data.name,
                            locked_protein_data.atom_count(),
                            locked_protein_data.chain_ids().len()
                        );
                        drop(locked_protein_data);
                        Ok(LuaProtein::new(shared_protein_handle.clone()))
                    }
                    Err(lua_runtime_error_message) => Err(mlua::Error::RuntimeError(lua_runtime_error_message)),
                }
            })?,
        )?;

        // pdb.load(path) loads a protein from a local file (PDB or mmCIF)
        let cloned_protein_store_reference = protein_data_store.clone();
        pdb_api_table.set(
            "load",
            lua_runtime_instance.create_function(move |_lua, requested_file_path: String| {
                println!("Loading file: {}...", requested_file_path);
                let mut locked_protein_store = cloned_protein_store_reference.write().unwrap();
                match locked_protein_store.load(&requested_file_path) {
                    Ok(shared_protein_handle) => {
                        let locked_protein_data = shared_protein_handle.read().unwrap();
                        println!(
                            "Loaded {} with {} atoms, {} chains",
                            locked_protein_data.name,
                            locked_protein_data.atom_count(),
                            locked_protein_data.chain_ids().len()
                        );
                        drop(locked_protein_data);
                        Ok(LuaProtein::new(shared_protein_handle.clone()))
                    }
                    Err(lua_runtime_error_message) => Err(mlua::Error::RuntimeError(lua_runtime_error_message)),
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

        // Simple print override for cleaner output
        let _ = lua_runtime_instance.globals().set(
            "print",
            lua_runtime_instance.create_function(|_, variadic_print_arguments: mlua::Variadic<mlua::Value>| {
                let formatted_output_strings: Vec<String> = variadic_print_arguments
                    .iter()
                    .map(|argument_value| match argument_value {
                        mlua::Value::Nil => "nil".to_string(),
                        mlua::Value::Boolean(boolean_value) => boolean_value.to_string(),
                        mlua::Value::Integer(integer_value) => integer_value.to_string(),
                        mlua::Value::Number(numeric_value) => format!("{:.4}", numeric_value),
                        mlua::Value::String(string_value) => string_value.to_str().unwrap_or("").to_string(),
                        _ => format!("{:?}", argument_value),
                    })
                    .collect();
                println!("[Lua] {}", formatted_output_strings.join("\t"));
                Ok(())
            })?,
        );

        Ok(Self { lua: lua_runtime_instance })
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
                lua_file_path,
                file_read_error_message
            ))),
        }
    }
}


