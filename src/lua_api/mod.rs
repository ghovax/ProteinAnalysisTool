mod protein;

use mlua::{Lua, Result};
use std::sync::{Arc, RwLock};

use crate::protein::ProteinStore;

pub use protein::LuaProtein;

pub struct ScriptEngine {
    lua: Lua,
}

impl ScriptEngine {
    pub fn new(store: Arc<RwLock<ProteinStore>>) -> Result<Self> {
        let lua = Lua::new();

        // Create pdb module
        let pdb = lua.create_table()?;

        // pdb.fetch(code) -> Protein
        let store_clone = store.clone();
        pdb.set(
            "fetch",
            lua.create_function(move |_lua, code: String| {
                println!("Fetching PDB: {}...", code);
                let mut store = store_clone.write().unwrap();
                match store.fetch(&code) {
                    Ok(protein_arc) => {
                        let protein = protein_arc.read().unwrap();
                        println!(
                            "Loaded {} with {} atoms, {} chains",
                            protein.name,
                            protein.atom_count(),
                            protein.chain_ids().len()
                        );
                        drop(protein);
                        Ok(LuaProtein::new(protein_arc.clone()))
                    }
                    Err(e) => Err(mlua::Error::RuntimeError(e)),
                }
            })?,
        )?;

        // pdb.load(path) -> Protein
        let store_clone = store.clone();
        pdb.set(
            "load",
            lua.create_function(move |_lua, path: String| {
                println!("Loading file: {}...", path);
                let mut store = store_clone.write().unwrap();
                match store.load(&path) {
                    Ok(protein_arc) => {
                        let protein = protein_arc.read().unwrap();
                        println!(
                            "Loaded {} with {} atoms, {} chains",
                            protein.name,
                            protein.atom_count(),
                            protein.chain_ids().len()
                        );
                        drop(protein);
                        Ok(LuaProtein::new(protein_arc.clone()))
                    }
                    Err(e) => Err(mlua::Error::RuntimeError(e)),
                }
            })?,
        )?;

        // pdb.list() -> table of names
        let store_clone = store.clone();
        pdb.set(
            "list",
            lua.create_function(move |lua, ()| {
                let store = store_clone.read().unwrap();
                let names = store.list();
                lua.create_sequence_from(names)
            })?,
        )?;

        lua.globals().set("pdb", pdb)?;

        // Simple print override for cleaner output
        lua.globals().set(
            "print",
            lua.create_function(|_, args: mlua::Variadic<mlua::Value>| {
                let output: Vec<String> = args
                    .iter()
                    .map(|v| match v {
                        mlua::Value::Nil => "nil".to_string(),
                        mlua::Value::Boolean(b) => b.to_string(),
                        mlua::Value::Integer(i) => i.to_string(),
                        mlua::Value::Number(n) => format!("{:.4}", n),
                        mlua::Value::String(s) => s.to_str().unwrap_or("").to_string(),
                        _ => format!("{:?}", v),
                    })
                    .collect();
                println!("[Lua] {}", output.join("\t"));
                Ok(())
            })?,
        )?;

        Ok(Self { lua })
    }

    pub fn run_script(&self, code: &str) -> Result<()> {
        self.lua.load(code).exec()
    }

    pub fn run_file(&self, path: &str) -> Result<()> {
        match std::fs::read_to_string(path) {
            Ok(code) => self.run_script(&code),
            Err(e) => Err(mlua::Error::RuntimeError(format!(
                "Failed to read {}: {}",
                path, e
            ))),
        }
    }
}
