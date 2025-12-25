//! Lua interface for protein data
//!
//! This module implements `mlua::UserData` for `LuaProtein`, allowing Lua scripts
//! to interact with and modify protein structures loaded in the application

use mlua::UserData;
use std::sync::{Arc, RwLock};

use crate::protein::{ColorScheme, ProteinData, Representation};

/// A Lua-exposed wrapper around a shared protein data object
#[derive(Clone)]
pub struct LuaProtein {
    inner: Arc<RwLock<ProteinData>>,
}

impl LuaProtein {
    /// Creates a new `LuaProtein` from a shared `ProteinData` handle
    pub fn new(inner: Arc<RwLock<ProteinData>>) -> Self {
        Self { inner }
    }
}

impl UserData for LuaProtein {
    fn add_methods<'lua, M: mlua::UserDataMethods<'lua, Self>>(methods: &mut M) {
        // p:name() returns the name/ID of the protein
        methods.add_method("name", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data.name.clone())
        });

        // p:atom_count() returns the total number of atoms in the protein
        methods.add_method("atom_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data.atom_count())
        });

        // p:chains() returns a list of all chain identifiers in the protein
        methods.add_method("chains", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let available_chain_identifiers = locked_protein_data.chain_ids();
            lua_context.create_sequence_from(available_chain_identifiers)
        });

        // p:center_of_mass() returns the center of mass of the protein
        methods.add_method("center_of_mass", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let center_of_mass_vector = locked_protein_data.center_of_mass();
            Ok((
                center_of_mass_vector.x,
                center_of_mass_vector.y,
                center_of_mass_vector.z,
            ))
        });

        // p:bounding_box() returns the axis-aligned bounding box of the protein
        methods.add_method("bounding_box", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let (minimum_coordinate_bound, maximum_coordinate_bound) =
                locked_protein_data.bounding_box();
            Ok((
                minimum_coordinate_bound.x,
                minimum_coordinate_bound.y,
                minimum_coordinate_bound.z,
                maximum_coordinate_bound.x,
                maximum_coordinate_bound.y,
                maximum_coordinate_bound.z,
            ))
        });

        // p:ca_count() returns the number of Alpha Carbon (CA) atoms
        methods.add_method("ca_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .len())
        });

        // p:residue_count() returns an approximate count of residues based on CA atoms
        methods.add_method("residue_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .len())
        });

        // p:show() and p:hide() control whether the protein is rendered
        methods.add_method_mut("show", |_, this, ()| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            mutable_protein_data.visible = true;
            Ok(())
        });

        methods.add_method_mut("hide", |_, this, ()| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            mutable_protein_data.visible = false;
            Ok(())
        });

        // p:info() returns a summary string with protein information
        methods.add_method("info", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let (minimum_coordinate_bound, maximum_coordinate_bound) = locked_protein_data.bounding_box();
            let bounding_box_dimensions = maximum_coordinate_bound - minimum_coordinate_bound;
            Ok(format!(
                "Protein: {}\n  Atoms: {}\n  Chains: {:?}\n  Residues: ~{}\n  Size: {:.1} x {:.1} x {:.1} A",
                locked_protein_data.name,
                locked_protein_data.atom_count(),
                locked_protein_data.chain_ids(),
                locked_protein_data.get_alpha_carbon_positions_and_chain_identifiers().len(),
                bounding_box_dimensions.x,
                bounding_box_dimensions.y,
                bounding_box_dimensions.z
            ))
        });

        // p:atoms() returns a table containing detailed information for every atom
        methods.add_method("atoms", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let lua_atoms_collection_table = lua_context.create_table()?;

            let mut atom_indexing_counter = 1;
            for current_atom_reference in locked_protein_data.pdb.atoms() {
                let lua_atom_data_table = lua_context.create_table()?;
                let atom_coordinates_tuple = current_atom_reference.pos();
                lua_atom_data_table.set("x", atom_coordinates_tuple.0)?;
                lua_atom_data_table.set("y", atom_coordinates_tuple.1)?;
                lua_atom_data_table.set("z", atom_coordinates_tuple.2)?;
                lua_atom_data_table.set("name", current_atom_reference.name())?;
                lua_atom_data_table.set(
                    "element",
                    current_atom_reference
                        .element()
                        .map(|element_reference| element_reference.symbol())
                        .unwrap_or("?"),
                )?;
                lua_atom_data_table.set("bfactor", current_atom_reference.b_factor())?;
                lua_atoms_collection_table.set(atom_indexing_counter, lua_atom_data_table)?;
                atom_indexing_counter += 1;
            }

            Ok(lua_atoms_collection_table)
        });

        // p:residues(chain_id) returns a table containing information for residues, optionally filtered by chain
        methods.add_method(
            "residues",
            |lua_context, this, chain_identifier_filter: Option<String>| {
                let locked_protein_data = this.inner.read().unwrap();
                let lua_residues_collection_table = lua_context.create_table()?;

                let mut residue_indexing_counter = 1;
                for current_chain_reference in locked_protein_data.pdb.chains() {
                    if let Some(ref filter_string) = chain_identifier_filter {
                        if current_chain_reference.id() != filter_string {
                            continue;
                        }
                    }

                    for current_residue_reference in current_chain_reference.residues() {
                        let lua_residue_data_table = lua_context.create_table()?;
                        lua_residue_data_table.set("chain", current_chain_reference.id())?;
                        lua_residue_data_table
                            .set("number", current_residue_reference.serial_number())?;

                        // Get residue name from first conformer
                        if let Some(current_conformer_reference) =
                            current_residue_reference.conformers().next()
                        {
                            lua_residue_data_table
                                .set("name", current_conformer_reference.name())?;
                        }

                        lua_residues_collection_table
                            .set(residue_indexing_counter, lua_residue_data_table)?;
                        residue_indexing_counter += 1;
                    }
                }

                Ok(lua_residues_collection_table)
            },
        );

        // p:representation(mode) sets the representation mode
        // Available modes are "spheres", "backbone", and "both"
        methods.add_method_mut(
            "representation",
            |_, this, requested_representation_mode: String| {
                let mut mutable_protein_data = this.inner.write().unwrap();
                mutable_protein_data.representation =
                    match requested_representation_mode.to_lowercase().as_str() {
                        "spheres" | "sphere" => Representation::Spheres,
                        "backbone" | "trace" | "line" | "lines" => Representation::Backbone,
                        "both" | "all" => Representation::BackboneAndSpheres,
                        _ => {
                            return Err(mlua::Error::RuntimeError(format!(
                        "Unknown representation: '{}'. Use 'spheres', 'backbone', or 'both'",
                        requested_representation_mode
                    )));
                        }
                    };
                Ok(())
            },
        );

        // p:color_by(scheme) sets the color scheme
        // Available schemes are "chain", "element", "bfactor", and "secondary"
        methods.add_method_mut("color_by", |_, this, requested_color_scheme_mode: String| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            mutable_protein_data.color_scheme = match requested_color_scheme_mode.to_lowercase().as_str() {
                "chain" | "chains" => ColorScheme::ByChain,
                "element" | "elements" | "atom" => ColorScheme::ByElement,
                "bfactor" | "b-factor" | "temperature" => ColorScheme::ByBFactor,
                "secondary" | "ss" | "structure" => ColorScheme::BySecondary,
                _ => {
                    return Err(mlua::Error::RuntimeError(format!(
                        "Unknown color scheme: '{}'. Use 'chain', 'element', 'bfactor', or 'secondary'",
                        requested_color_scheme_mode
                    )));
                }
            };
            Ok(())
        });

        // p:color(r, g, b) sets a uniform RGB color for the entire protein
        methods.add_method_mut(
            "color",
            |_,
             this,
             (red_color_component, green_color_component, blue_color_component): (
                f32,
                f32,
                f32,
            )| {
                let mut mutable_protein_data = this.inner.write().unwrap();
                mutable_protein_data.color_scheme = ColorScheme::Uniform([
                    red_color_component,
                    green_color_component,
                    blue_color_component,
                ]);
                Ok(())
            },
        );
    }
}
