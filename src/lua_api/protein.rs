//! Lua interface for protein data
//! 
//! This module implements `mlua::UserData` for `LuaProtein`, allowing Lua scripts
//! to interact with and modify protein structures loaded in the application

use mlua::UserData;
use std::sync::{Arc, RwLock};

use crate::protein::{ColorScheme, ProteinData, Representation};
use crate::lua_api::selection::LuaSelection;

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
        // protein:get_protein_name() returns the name/ID of the protein
        methods.add_method("get_protein_name", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data.name.clone())
        });

        // protein:get_total_atom_count() returns the total number of atoms in the protein
        methods.add_method("get_total_atom_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data.atom_count())
        });

        // protein:get_chain_identifiers() returns a list of all chain identifiers in the protein
        methods.add_method("get_chain_identifiers", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let available_chain_identifiers = locked_protein_data.chain_ids();
            lua_context.create_sequence_from(available_chain_identifiers)
        });

        // protein:calculate_center_of_mass() returns the center of mass of the protein
        methods.add_method("calculate_center_of_mass", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let center_of_mass_vector = locked_protein_data.center_of_mass();
            Ok((
                center_of_mass_vector.x,
                center_of_mass_vector.y,
                center_of_mass_vector.z,
            ))
        });

        // protein:calculate_bounding_box_dimensions() returns the axis-aligned bounding box of the protein
        methods.add_method("calculate_bounding_box_dimensions", |_, this, ()| {
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

        // protein:get_alpha_carbon_count() returns the number of Alpha Carbon (CA) atoms
        methods.add_method("get_alpha_carbon_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .len())
        });

        // protein:get_total_residue_count() returns an approximate count of residues based on CA atoms
        methods.add_method("get_total_residue_count", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(locked_protein_data
                .get_alpha_carbon_positions_and_chain_identifiers()
                .len())
        });

        // protein:set_visibility_on() and protein:set_visibility_off() control whether the protein is rendered
        methods.add_method_mut("set_visibility_on", |_, this, ()| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            mutable_protein_data.visible = true;
            Ok(())
        });

        methods.add_method_mut("set_visibility_off", |_, this, ()| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            mutable_protein_data.visible = false;
            Ok(())
        });

        // protein:get_summary_information() returns a summary string with protein information
        methods.add_method("get_summary_information", |_, this, ()| {
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

        // protein:get_all_atom_data() returns a table containing detailed information for every atom
        methods.add_method("get_all_atom_data", |lua_context, this, ()| {
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

        // protein:get_all_residue_data(chain_id) returns a table containing information for residues, optionally filtered by chain
        methods.add_method(
            "get_all_residue_data",
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

        // protein:select_atoms_by_query(query) returns a LuaSelection object matching the PyMOL-style selection query
        methods.add_method("select_atoms_by_query", |_, this, selection_query_string: String| {
            let locked_protein_data = this.inner.read().unwrap();
            let selection_set = locked_protein_data.select_atoms_from_string(&selection_query_string)
                .map_err(|error_message| mlua::Error::RuntimeError(error_message))?;
            
            Ok(LuaSelection::new(selection_set, this.inner.clone()))
        });

        // protein:calculate_solvent_accessible_surface_area() returns the total Solvent Accessible Surface Area in square Angstroms
        methods.add_method("calculate_solvent_accessible_surface_area", |_, this, (probe_radius, points): (Option<f32>, Option<usize>)| {
            let locked_protein_data = this.inner.read().unwrap();
            let probe_radius_value = probe_radius.unwrap_or(1.4);
            let number_of_points = points.unwrap_or(100);
            let sasa_value = crate::analysis::sasa::calculate_protein_solvent_accessible_surface_area(
                &locked_protein_data,
                probe_radius_value,
                number_of_points
            );
            Ok(sasa_value)
        });

        // protein:export_analysis_report_to_json(path) exports analysis data to a JSON file
        methods.add_method("export_analysis_report_to_json", |_, this, file_system_path: String| {
            let locked_protein_data = this.inner.read().unwrap();
            crate::analysis::export::export_protein_analysis_to_json(&locked_protein_data, &file_system_path)
                .map_err(|error_message| mlua::Error::RuntimeError(error_message))
        });

        // protein:export_residue_data_to_csv(path) exports analysis data to a CSV file
        methods.add_method("export_residue_data_to_csv", |_, this, file_system_path: String| {
            let locked_protein_data = this.inner.read().unwrap();
            crate::analysis::export::export_residue_data_to_csv(&locked_protein_data, &file_system_path)
                .map_err(|error_message| mlua::Error::RuntimeError(error_message))
        });

        // protein:get_sequence_for_chain(chain_id) returns the primary sequence of a chain
        methods.add_method("get_sequence_for_chain", |_, this, chain_identifier: String| {
            let locked_protein_data = this.inner.read().unwrap();
            crate::analysis::sequence::extract_sequence_from_chain(&locked_protein_data, &chain_identifier)
                .map_err(|error_message| mlua::Error::RuntimeError(error_message))
        });

        // protein:generate_fasta_formatted_sequence() returns the protein sequence in FASTA format
        methods.add_method("generate_fasta_formatted_sequence", |_, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            Ok(crate::analysis::sequence::generate_fasta_formatted_string(&locked_protein_data))
        });

        // protein:detect_salt_bridge_interactions() returns a list of detected salt bridges
        methods.add_method("detect_salt_bridge_interactions", |lua_context, this, distance_threshold: Option<f32>| {
            let locked_protein_data = this.inner.read().unwrap();
            let threshold_value = distance_threshold.unwrap_or(4.0);
            let detected_salt_bridges = crate::analysis::interactions::detect_salt_bridge_interactions(&locked_protein_data, threshold_value);
            
            let lua_salt_bridges_table = lua_context.create_table()?;
            for (bond_index, current_bridge) in detected_salt_bridges.into_iter().enumerate() {
                let bridge_data_table = lua_context.create_table()?;
                bridge_data_table.set("donor_index", current_bridge.donor_atom_index)?;
                bridge_data_table.set("acceptor_index", current_bridge.acceptor_atom_index)?;
                bridge_data_table.set("distance", current_bridge.inter_atomic_distance)?;
                lua_salt_bridges_table.set(bond_index + 1, bridge_data_table)?;
            }
            Ok(lua_salt_bridges_table)
        });

        // protein:perform_secondary_structure_assignment() re-calculates secondary structure assignment
        methods.add_method_mut("perform_secondary_structure_assignment", |_, this, ()| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            crate::analysis::dssp::assign_secondary_structure_to_protein(&mut mutable_protein_data);
            Ok(())
        });

        // protein:calculate_root_mean_square_fluctuation(selection) calculates Root-Mean-Square Fluctuation across all models for a selection
        methods.add_method("calculate_root_mean_square_fluctuation", |lua_context, this, selection: LuaSelection| {
            let locked_protein_data = this.inner.read().unwrap();
            
            // Collect coordinates for selected atoms across all models
            let mut model_coordinates_collection = Vec::new();
            for current_model in locked_protein_data.pdb.models() {
                let mut atom_positions_in_model = Vec::new();
                // This is slightly complex because selection indices are global or relative to a specific model?
                // Assuming selection indices are global indices into pdb.atoms()
                for &atom_index in &selection.inner_selection_set.atom_indices {
                    if let Some(atom_reference) = current_model.atoms().nth(atom_index) {
                        let pos = atom_reference.pos();
                        atom_positions_in_model.push(glam::Vec3::new(pos.0 as f32, pos.1 as f32, pos.2 as f32));
                    }
                }
                if !atom_positions_in_model.is_empty() {
                    model_coordinates_collection.push(atom_positions_in_model);
                }
            }

            if model_coordinates_collection.is_empty() {
                return Ok(lua_context.create_table()?);
            }

            let model_count = model_coordinates_collection.len();
            let atom_count = model_coordinates_collection[0].len();
            let rmsf_results_table = lua_context.create_table()?;

            for atom_index in 0..atom_count {
                let mut average_position_vector = glam::Vec3::ZERO;
                for model_index in 0..model_count {
                    average_position_vector += model_coordinates_collection[model_index][atom_index];
                }
                average_position_vector /= model_count as f32;

                let mut sum_of_squared_deviations = 0.0;
                for model_index in 0..model_count {
                    sum_of_squared_deviations += model_coordinates_collection[model_index][atom_index].distance_squared(average_position_vector);
                }
                let rmsf_value = (sum_of_squared_deviations / model_count as f32).sqrt();
                rmsf_results_table.set(atom_index + 1, rmsf_value)?;
            }

            Ok(rmsf_results_table)
        });

                // protein:set_representation_mode(mode) sets the representation mode
                // Available modes are "spheres", "backbone_trace", "backbone_and_spheres", "sticks", "ball_and_stick", "space_filling", and "lines"
                methods.add_method_mut(
                    "set_representation_mode",
                    |_, this, requested_representation_mode: String| {
                        let mut mutable_protein_data = this.inner.write().unwrap();
                        mutable_protein_data.representation =
                            match requested_representation_mode.to_lowercase().as_str() {
                                "spheres" => Representation::Spheres,
                                "backbone_trace" => Representation::Backbone,
                                "backbone_and_spheres" => Representation::BackboneAndSpheres,
                                "sticks" => Representation::Sticks,
                                "ball_and_stick" => Representation::BallAndStick,
                                "space_filling" => Representation::SpaceFilling,
                                "lines" => Representation::Lines,
                                _ => {
                                    return Err(mlua::Error::RuntimeError(format!(
                                "Unknown representation: '{}'. Use 'spheres', 'backbone_trace', 'sticks', 'ball_and_stick', etc.",
                                requested_representation_mode
                            )));
                                }
                            };
                        Ok(())
                    },
                );
                // protein:calculate_ramachandran_dihedral_angles() returns a list of Phi/Psi points for all residues
        methods.add_method("calculate_ramachandran_dihedral_angles", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let backbone_dihedral_results_collection = crate::analysis::dihedrals::calculate_all_backbone_dihedrals(&locked_protein_data);
            
            let lua_ramachandran_points_table = lua_context.create_table()?;
            let mut result_counting_index = 1;

            for (chain_identifier, residue_number, dihedral_angles) in backbone_dihedral_results_collection {
                if let (Some(phi_angle), Some(psi_angle)) = (dihedral_angles.phi_angle, dihedral_angles.psi_angle) {
                    let lua_point_data_table = lua_context.create_table()?;
                    lua_point_data_table.set("phi", phi_angle)?;
                    lua_point_data_table.set("psi", psi_angle)?;
                    lua_point_data_table.set("residue_number", residue_number)?;
                    lua_point_data_table.set("chain", chain_identifier)?;
                    
                    lua_ramachandran_points_table.set(result_counting_index, lua_point_data_table)?;
                    result_counting_index += 1;
                }
            }
            Ok(lua_ramachandran_points_table)
        });

        // protein:detect_geometric_hydrogen_bonds() returns a list of detected hydrogen bonds
        methods.add_method("detect_geometric_hydrogen_bonds", |lua_context, this, ()| {
            let locked_protein_data = this.inner.read().unwrap();
            let identified_hydrogen_bonds = crate::analysis::hydrogen_bonds::detect_hydrogen_bonds_in_protein(&locked_protein_data);
            
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

        // protein:calculate_molecular_surface(resolution, threshold) generates a molecular surface mesh
        methods.add_method_mut("calculate_molecular_surface", |_, this, (voxel_resolution, isosurface_threshold): (Option<f32>, Option<f32>)| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            let resolution_value = voxel_resolution.unwrap_or(0.5);
            let threshold_value = isosurface_threshold.unwrap_or(0.0);
            
            // These tables will be provided by the user as per request. 
            // For now, I'll use empty/dummy ones to keep the logic structure correct.
            let edge_intersection_lookup_table = [0u32; 256];
            let triangulation_lookup_table = [[-1i32; 16]; 256];

            mutable_protein_data.compute_molecular_surface_mesh(
                resolution_value, 
                threshold_value,
                &edge_intersection_lookup_table,
                &triangulation_lookup_table
            );
            Ok(())
        });

        // protein:set_surface_visibility(visible) controls whether the molecular surface is rendered
        methods.add_method_mut("set_surface_visibility", |_, this, should_be_visible: bool| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            mutable_protein_data.molecular_surface_mesh.is_surface_visible = should_be_visible;
            Ok(())
        });

        // protein:calculate_root_mean_square_deviation(other_protein) calculates the RMSD between two proteins based on Alpha Carbons from the first model
        // protein:calculate_root_mean_square_deviation(other_protein, selection) calculates the RMSD between two proteins, optionally restricted to a selection
        methods.add_method("calculate_root_mean_square_deviation", |_, this, (other_lua_protein, optional_selection): (LuaProtein, Option<LuaSelection>)| {
            let reference_protein_locked_data = this.inner.read().unwrap();
            let moving_protein_locked_data = other_lua_protein.inner.read().unwrap();
            
            let reference_alpha_carbon_positions = if let Some(selection) = optional_selection {
                let mut selected_positions = Vec::new();
                for &atom_index in &selection.inner_selection_set.atom_indices {
                    let atom_reference = reference_protein_locked_data.pdb.atoms().nth(atom_index).unwrap();
                    let atom_position_tuple = atom_reference.pos();
                    selected_positions.push(glam::Vec3::new(atom_position_tuple.0 as f32, atom_position_tuple.1 as f32, atom_position_tuple.2 as f32));
                }
                selected_positions
            } else {
                reference_protein_locked_data.get_alpha_carbon_positions_for_first_model()
            };
                
            let moving_alpha_carbon_positions = moving_protein_locked_data
                .get_alpha_carbon_positions_for_first_model();
            
            if reference_alpha_carbon_positions.len() != moving_alpha_carbon_positions.len() {
                return Err(mlua::Error::RuntimeError(format!(
                    "Coordinate sets must have the same length for RMSD calculation ({} vs {} atoms).",
                    reference_alpha_carbon_positions.len(),
                    moving_alpha_carbon_positions.len()
                )));
            }

            crate::analysis::rmsd::calculate_rmsd_between_coordinate_sets(
                &reference_alpha_carbon_positions,
                &moving_alpha_carbon_positions
            ).map_err(|error_message| mlua::Error::RuntimeError(error_message))
        });

        // protein:superimpose_onto_reference_structure(reference_protein, selection) superimposes this protein onto a reference using Kabsch algorithm
        methods.add_method_mut("superimpose_onto_reference_structure", |_, this, (reference_lua_protein, optional_selection): (LuaProtein, Option<LuaSelection>)| {
            let mut moving_protein_mutable_data = this.inner.write().unwrap();
            let reference_protein_locked_data = reference_lua_protein.inner.read().unwrap();
            
            let (reference_alpha_carbon_positions, moving_alpha_carbon_positions) = if let Some(selection) = optional_selection {
                let mut reference_positions = Vec::new();
                let mut moving_positions = Vec::new();
                for &atom_index in &selection.inner_selection_set.atom_indices {
                    let reference_atom_handle = reference_protein_locked_data.pdb.atoms().nth(atom_index).unwrap();
                    let moving_atom_handle = moving_protein_mutable_data.pdb.atoms().nth(atom_index).unwrap();
                    let reference_atom_position_tuple = reference_atom_handle.pos();
                    let moving_atom_position_tuple = moving_atom_handle.pos();
                    reference_positions.push(glam::Vec3::new(reference_atom_position_tuple.0 as f32, reference_atom_position_tuple.1 as f32, reference_atom_position_tuple.2 as f32));
                    moving_positions.push(glam::Vec3::new(moving_atom_position_tuple.0 as f32, moving_atom_position_tuple.1 as f32, moving_atom_position_tuple.2 as f32));
                }
                (reference_positions, moving_positions)
            } else {
                (
                    reference_protein_locked_data.get_alpha_carbon_positions_for_first_model(),
                    moving_protein_mutable_data.get_alpha_carbon_positions_for_first_model()
                )
            };
                
            if reference_alpha_carbon_positions.len() != moving_alpha_carbon_positions.len() {
                return Err(mlua::Error::RuntimeError(format!(
                    "Proteins must have the same number of Alpha Carbons for superposition ({} vs {} atoms).",
                    reference_alpha_carbon_positions.len(),
                    moving_alpha_carbon_positions.len()
                )));
            }
            
            // Calculate centers of mass for centering
            let reference_center_of_mass = reference_protein_locked_data.center_of_mass();
            let moving_center_of_mass = moving_protein_mutable_data.center_of_mass();
            
            let reference_centered_positions: Vec<_> = reference_alpha_carbon_positions.iter() 
                .map(|&position_vector| position_vector - reference_center_of_mass)
                .collect();
            let moving_centered_positions: Vec<_> = moving_alpha_carbon_positions.iter()
                .map(|&position_vector| position_vector - moving_center_of_mass)
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
                )).unwrap();
            }
            
            Ok(())
        });

        // protein:set_color_scheme_by_property(scheme) sets the color scheme
        // Available schemes are "chain_identifier", "chemical_element", "bfactor_value", and "secondary_structure"
        methods.add_method_mut("set_color_scheme_by_property", |_, this, requested_color_scheme_mode: String| {
            let mut mutable_protein_data = this.inner.write().unwrap();
            mutable_protein_data.color_scheme = match requested_color_scheme_mode.to_lowercase().as_str() {
                "chain_identifier" => ColorScheme::ByChain,
                "chemical_element" => ColorScheme::ByElement,
                "bfactor_value" => ColorScheme::ByBFactor,
                "secondary_structure" => ColorScheme::BySecondary,
                _ => {
                    return Err(mlua::Error::RuntimeError(format!(
                        "Unknown color scheme: '{}'. Use 'chain_identifier', 'chemical_element', 'bfactor_value', or 'secondary_structure'",
                        requested_color_scheme_mode
                    )));
                }
            };
            Ok(())
        });

        // protein:set_uniform_rgb_color(r, g, b) sets a uniform RGB color for the entire protein
        methods.add_method_mut(
            "set_uniform_rgb_color",
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