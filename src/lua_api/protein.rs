use mlua::UserData;
use std::sync::{Arc, RwLock};

use crate::protein::{ProteinData, Representation, ColorScheme};

#[derive(Clone)]
pub struct LuaProtein {
    inner: Arc<RwLock<ProteinData>>,
}

impl LuaProtein {
    pub fn new(inner: Arc<RwLock<ProteinData>>) -> Self {
        Self { inner }
    }
}

impl UserData for LuaProtein {
    fn add_methods<'lua, M: mlua::UserDataMethods<'lua, Self>>(methods: &mut M) {
        // p:name() -> string
        methods.add_method("name", |_, this, ()| {
            let protein = this.inner.read().unwrap();
            Ok(protein.name.clone())
        });

        // p:atom_count() -> number
        methods.add_method("atom_count", |_, this, ()| {
            let protein = this.inner.read().unwrap();
            Ok(protein.atom_count())
        });

        // p:chains() -> table of chain IDs
        methods.add_method("chains", |lua, this, ()| {
            let protein = this.inner.read().unwrap();
            let chains = protein.chain_ids();
            lua.create_sequence_from(chains)
        });

        // p:center_of_mass() -> x, y, z
        methods.add_method("center_of_mass", |_, this, ()| {
            let protein = this.inner.read().unwrap();
            let com = protein.center_of_mass();
            Ok((com.x, com.y, com.z))
        });

        // p:bounding_box() -> min_x, min_y, min_z, max_x, max_y, max_z
        methods.add_method("bounding_box", |_, this, ()| {
            let protein = this.inner.read().unwrap();
            let (min, max) = protein.bounding_box();
            Ok((min.x, min.y, min.z, max.x, max.y, max.z))
        });

        // p:ca_count() -> number of alpha carbons
        methods.add_method("ca_count", |_, this, ()| {
            let protein = this.inner.read().unwrap();
            Ok(protein.ca_positions().len())
        });

        // p:residue_count() -> approximate residue count (via CA atoms)
        methods.add_method("residue_count", |_, this, ()| {
            let protein = this.inner.read().unwrap();
            Ok(protein.ca_positions().len())
        });

        // p:show() / p:hide() - visibility control
        methods.add_method_mut("show", |_, this, ()| {
            let mut protein = this.inner.write().unwrap();
            protein.visible = true;
            Ok(())
        });

        methods.add_method_mut("hide", |_, this, ()| {
            let mut protein = this.inner.write().unwrap();
            protein.visible = false;
            Ok(())
        });

        // p:info() -> formatted string
        methods.add_method("info", |_, this, ()| {
            let protein = this.inner.read().unwrap();
            let (min, max) = protein.bounding_box();
            let size = max - min;
            Ok(format!(
                "Protein: {}\n  Atoms: {}\n  Chains: {:?}\n  Residues: ~{}\n  Size: {:.1} x {:.1} x {:.1} A",
                protein.name,
                protein.atom_count(),
                protein.chain_ids(),
                protein.ca_positions().len(),
                size.x,
                size.y,
                size.z
            ))
        });

        // p:atoms() -> iterator (returns table of atom data for now)
        methods.add_method("atoms", |lua, this, ()| {
            let protein = this.inner.read().unwrap();
            let atoms_table = lua.create_table()?;

            let mut idx = 1;
            for atom in protein.pdb.atoms() {
                let atom_table = lua.create_table()?;
                let pos = atom.pos();
                atom_table.set("x", pos.0)?;
                atom_table.set("y", pos.1)?;
                atom_table.set("z", pos.2)?;
                atom_table.set("name", atom.name())?;
                atom_table.set("element", atom.element().map(|e| e.symbol()).unwrap_or("?"))?;
                atom_table.set("bfactor", atom.b_factor())?;
                atoms_table.set(idx, atom_table)?;
                idx += 1;
            }

            Ok(atoms_table)
        });

        // p:residues(chain_id) -> table of residue info
        methods.add_method("residues", |lua, this, chain_id: Option<String>| {
            let protein = this.inner.read().unwrap();
            let residues_table = lua.create_table()?;

            let mut idx = 1;
            for chain in protein.pdb.chains() {
                if let Some(ref filter) = chain_id {
                    if chain.id() != filter {
                        continue;
                    }
                }

                for residue in chain.residues() {
                    let res_table = lua.create_table()?;
                    res_table.set("chain", chain.id())?;
                    res_table.set("number", residue.serial_number())?;

                    // Get residue name from first conformer
                    if let Some(conformer) = residue.conformers().next() {
                        res_table.set("name", conformer.name())?;
                    }

                    residues_table.set(idx, res_table)?;
                    idx += 1;
                }
            }

            Ok(residues_table)
        });

        // p:representation(mode) - set representation mode
        // modes: "spheres", "backbone", "both"
        methods.add_method_mut("representation", |_, this, mode: String| {
            let mut protein = this.inner.write().unwrap();
            protein.representation = match mode.to_lowercase().as_str() {
                "spheres" | "sphere" => Representation::Spheres,
                "backbone" | "trace" | "line" | "lines" => Representation::Backbone,
                "both" | "all" => Representation::BackboneAndSpheres,
                _ => {
                    return Err(mlua::Error::RuntimeError(format!(
                        "Unknown representation: '{}'. Use 'spheres', 'backbone', or 'both'",
                        mode
                    )));
                }
            };
            Ok(())
        });

        // p:color_by(scheme) - set color scheme
        // schemes: "chain", "element", "bfactor", "secondary"
        methods.add_method_mut("color_by", |_, this, scheme: String| {
            let mut protein = this.inner.write().unwrap();
            protein.color_scheme = match scheme.to_lowercase().as_str() {
                "chain" | "chains" => ColorScheme::ByChain,
                "element" | "elements" | "atom" => ColorScheme::ByElement,
                "bfactor" | "b-factor" | "temperature" => ColorScheme::ByBFactor,
                "secondary" | "ss" | "structure" => ColorScheme::BySecondary,
                _ => {
                    return Err(mlua::Error::RuntimeError(format!(
                        "Unknown color scheme: '{}'. Use 'chain', 'element', 'bfactor', or 'secondary'",
                        scheme
                    )));
                }
            };
            Ok(())
        });

        // p:color(r, g, b) - set uniform color
        methods.add_method_mut("color", |_, this, (r, g, b): (f32, f32, f32)| {
            let mut protein = this.inner.write().unwrap();
            protein.color_scheme = ColorScheme::Uniform([r, g, b]);
            Ok(())
        });
    }
}
