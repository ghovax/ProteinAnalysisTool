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

impl<'lua> mlua::FromLua<'lua> for LuaProtein {
    fn from_lua(value: mlua::Value<'lua>, _lua: &'lua mlua::Lua) -> mlua::Result<Self> {
        match value {
            mlua::Value::UserData(ud) => Ok(ud.borrow::<Self>()?.clone()),
            _ => Err(mlua::Error::FromLuaConversionError {
                from: value.type_name(),
                to: "LuaProtein",
                message: Some("expected a LuaProtein object".to_string()),
            }),
        }
    }
}

impl UserData for LuaProtein {
    fn add_methods<'lua, M: mlua::UserDataMethods<'lua, Self>>(methods: &mut M) {
        // protein:name() returns the name/ID of the protein
        methods.add_method("name", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data.name.clone())
        });

        // protein:atom_count() returns the total number of atoms in the protein
        methods.add_method("atom_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data.atom_count())
        });

        // protein:chains() returns a list of all chain identifiers in the protein
        methods.add_method("chains", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let available_chain_identifiers = locked_protein_data.chain_ids();
            lua_context.create_sequence_from(available_chain_identifiers)
        });

        // protein:center_of_mass() returns the center of mass of the protein
        methods.add_method("center_of_mass", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let center_of_mass_vector = locked_protein_data.center_of_mass();
            Ok((
                center_of_mass_vector.x,
                center_of_mass_vector.y,
                center_of_mass_vector.z,
            ))
        });

        // protein:bounding_box() returns the axis-aligned bounding box of the protein
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

        // protein:ca_count() returns the number of Alpha Carbon (CA) atoms
        methods.add_method("ca_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .len())
        });

        // protein:residue_count() returns an approximate count of residues based on CA atoms
        methods.add_method("residue_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .len())
        });

        // protein:show() and protein:hide() control whether the protein is rendered
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

        // protein:info() returns a summary string with protein information
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

        // protein:atoms() returns a table containing detailed information for every atom
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

        // protein:residues(chain_id) returns a table containing information for residues, optionally filtered by chain
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

        // protein:select(query) returns a list of atom indices matching the PyMOL-style selection query
        methods.add_method("select", |lua_context, this, selection_query_string: String| {
            let locked_protein_data = this.inner.read().unwrap();
            let selection_set = locked_protein_data.select_atoms_from_string(&selection_query_string)
                .map_err(|error_message| mlua::Error::RuntimeError(error_message))?;
            
            let lua_selected_indices_table = lua_context.create_table()?;
            for (selection_index, atom_index) in selection_set.atom_indices.iter().enumerate() {
                lua_selected_indices_table.set(selection_index + 1, *atom_index)?;
            }
            Ok(lua_selected_indices_table)
        });

        // protein:representation(mode) sets the representation mode
        // Available modes are "spheres", "backbone", "both", "sticks", "ball-and-stick", "space-filling", and "lines"
        methods.add_method_mut(
            "representation",
            |_, this, requested_representation_mode: String| {
                let mut mutable_protein_data = this.inner.write().unwrap();
                mutable_protein_data.representation =
                    match requested_representation_mode.to_lowercase().as_str() {
                        "spheres" | "sphere" => Representation::Spheres,
                        "backbone" | "trace" | "line_trace" => Representation::Backbone,
                        "both" | "all" | "backbone_and_spheres" => Representation::BackboneAndSpheres,
                        "sticks" | "stick" | "cylinders" => Representation::Sticks,
                        "ball_and_stick" | "ball-and-stick" => Representation::BallAndStick,
                        "space_filling" | "space-filling" | "vdw" => Representation::SpaceFilling,
                        "lines" | "line" | "wireframe" => Representation::Lines,
                        _ => {
                            return Err(mlua::Error::RuntimeError(format!(
                        "Unknown representation: '{}'. Use 'spheres', 'backbone', 'sticks', 'ball-and-stick', etc.",
                        requested_representation_mode
                    )));
                        }
                    };
                Ok(())
            },
        );

        // protein:ramachandran_data() returns a list of Phi/Psi points for all residues
        methods.add_method("ramachandran_data", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let ramachandran_points_collection = crate::analysis::ramachandran::calculate_ramachandran_angles_for_protein(&locked_protein_data);
            
            let lua_ramachandran_points_table = lua_context.create_table()?;
            for (point_indexing_counter, current_point) in ramachandran_points_collection.into_iter().enumerate() {
                let lua_point_data_table = lua_context.create_table()?;
                lua_point_data_table.set("phi", current_point.phi_angle)?;
                lua_point_data_table.set("psi", current_point.psi_angle)?;
                lua_point_data_table.set("residue_name", current_point.residue_name)?;
                lua_point_data_table.set("residue_number", current_point.residue_number)?;
                lua_point_data_table.set("chain", current_point.chain_identifier)?;
                
                lua_ramachandran_points_table.set(point_indexing_counter + 1, lua_point_data_table)?;
            }
            Ok(lua_ramachandran_points_table)
        });

        // protein:hydrogen_bonds() returns a list of detected hydrogen bonds
        methods.add_method("hydrogen_bonds", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let identified_hydrogen_bonds = crate::analysis::hbonds::detect_hydrogen_bonds_in_protein(&locked_protein_data);
            
            let lua_hydrogen_bonds_table = lua_context.create_table()?;
            for (bond_indexing_counter, current_bond) in identified_hydrogen_bonds.into_iter().enumerate() {
                let lua_bond_data_table = lua_context.create_table()?;
                lua_bond_data_table.set("donor_index", current_bond.donor_atom_index)?;
                lua_bond_data_table.set("acceptor_index", current_bond.acceptor_atom_index)?;
                lua_bond_data_table.set("distance", current_bond.distance_in_angstroms)?;
                
                lua_hydrogen_bonds_table.set(bond_indexing_counter + 1, lua_bond_data_table)?;
            }
            Ok(lua_hydrogen_bonds_table)
        });

        // protein:rmsd(other_protein) calculates the RMSD between two proteins based on Alpha Carbons
        methods.add_method("rmsd", |_, this, other_lua_protein: LuaProtein| {
            let reference_protein_locked_data = this.inner.read().unwrap();
            let moving_protein_locked_data = other_lua_protein.inner.read().unwrap();
            
            let reference_alpha_carbon_positions: Vec<_> = reference_protein_locked_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .into_iter()
                .map(|(position, _)| position)
                .collect();
                
            let moving_alpha_carbon_positions: Vec<_> = moving_protein_locked_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .into_iter()
                .map(|(position, _)| position)
                .collect();
            
            crate::analysis::rmsd::calculate_rmsd_between_coordinate_sets(
                &reference_alpha_carbon_positions,
                &moving_alpha_carbon_positions
            ).map_err(|error_message| mlua::Error::RuntimeError(error_message))
        });

        // protein:superpose(reference_protein) superimposes this protein onto a reference using Kabsch algorithm
        methods.add_method_mut("superpose", |_, this, reference_lua_protein: LuaProtein| {
            let mut moving_protein_mutable_data = this.inner.write().unwrap();
            let reference_protein_locked_data = reference_lua_protein.inner.read().unwrap();
            
            let reference_alpha_carbon_positions: Vec<_> = reference_protein_locked_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .into_iter()
                .map(|(position, _)| position)
                .collect();
                
            let moving_alpha_carbon_positions: Vec<_> = moving_protein_mutable_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .into_iter()
                .map(|(position, _)| position)
                .collect();
                
            if reference_alpha_carbon_positions.len() != moving_alpha_carbon_positions.len() {
                return Err(mlua::Error::RuntimeError("Proteins must have the same number of Alpha Carbons for superposition".to_string()));
            }
            
            // Calculate centers of mass for centering
            let reference_center_of_mass = reference_protein_locked_data.center_of_mass();
            let moving_center_of_mass = moving_protein_mutable_data.center_of_mass();
            
            let reference_centered_positions: Vec<_> = reference_alpha_carbon_positions.iter()
                .map(|&pos| pos - reference_center_of_mass)
                .collect();
            let moving_centered_positions: Vec<_> = moving_alpha_carbon_positions.iter()
                .map(|&pos| pos - moving_center_of_mass)
                .collect();
            
            let optimal_rotation_matrix = crate::analysis::rmsd::compute_kabsch_optimal_rotation(
                &reference_centered_positions,
                &moving_centered_positions
            ).map_err(|error_message| mlua::Error::RuntimeError(error_message))?;
            
            // Apply transformation to all atoms in the moving protein
            for current_atom_reference in moving_protein_mutable_data.pdb.atoms_mut() {
                let current_position_tuple = current_atom_reference.pos();
                let current_position_vector = glam::Vec3::new(
                    current_position_tuple.0 as f32,
                    current_position_tuple.1 as f32,
                    current_position_tuple.2 as f32
                );
                
                // Centering, rotating, and moving to reference center
                let transformed_position_vector = optimal_rotation_matrix * (current_position_vector - moving_center_of_mass) + reference_center_of_mass;
                
                current_atom_reference.set_pos((
                    transformed_position_vector.x as f64,
                    transformed_position_vector.y as f64,
                    transformed_position_vector.z as f64
                ));
            }
            
            Ok(())
        });

        // protein:color_by(scheme) sets the color scheme
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

        // protein:color(r, g, b) sets a uniform RGB color for the entire protein
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
