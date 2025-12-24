use glam::Vec3;
use pdbtbx::{PDB, ReadOptions, Format};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use super::fetch::{fetch_pdb, load_file, FileFormat};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Representation {
    #[default]
    Spheres,
    Backbone,
    BackboneAndSpheres,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum ColorScheme {
    #[default]
    ByChain,
    ByElement,
    ByBFactor,
    BySecondary,
    Uniform([f32; 3]),
}

pub struct ProteinData {
    pub pdb: PDB,
    pub name: String,
    pub visible: bool,
    pub representation: Representation,
    pub color_scheme: ColorScheme,
}

impl ProteinData {
    pub fn from_string(file_content_string: &str, name: &str, format: FileFormat) -> Result<Self, String> {
        let target_file_format = match format {
            FileFormat::Pdb => Format::Pdb,
            FileFormat::Cif => Format::Mmcif,
        };

        let (pdb, parsing_warning_messages) = ReadOptions::default()
            .set_level(pdbtbx::StrictnessLevel::Loose)
            .set_format(target_file_format)
            .read_raw(std::io::BufReader::new(std::io::Cursor::new(file_content_string.as_bytes())))
            .map_err(|pdb_parse_error| format!("Parse error: {:?}", pdb_parse_error))?;

        if !parsing_warning_messages.is_empty() {
            eprintln!("Warnings while parsing {}: {:?}", name, parsing_warning_messages);
        }

        Ok(Self {
            pdb,
            name: name.to_string(),
            visible: true,
            representation: Representation::default(),
            color_scheme: ColorScheme::default(),
        })
    }

    pub fn atom_count(&self) -> usize {
        self.pdb.atom_count()
    }

    pub fn chain_ids(&self) -> Vec<String> {
        self.pdb
            .chains()
            .map(|current_chain| current_chain.id().to_string())
            .collect()
    }

    pub fn center_of_mass(&self) -> Vec3 {
        let mut position_sum_vector = Vec3::ZERO;
        let mut atom_counting_index = 0;

        for current_atom_reference in self.pdb.atoms() {
            let atom_coordinates_tuple = current_atom_reference.pos();
            position_sum_vector += Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32);
            atom_counting_index += 1;
        }

        if atom_counting_index > 0 {
            position_sum_vector / atom_counting_index as f32
        } else {
            Vec3::ZERO
        }
    }

    pub fn bounding_box(&self) -> (Vec3, Vec3) {
        let mut minimum_coordinate_bound = Vec3::splat(f32::MAX);
        let mut maximum_coordinate_bound = Vec3::splat(f32::MIN);

        for current_atom_reference in self.pdb.atoms() {
            let atom_coordinates_tuple = current_atom_reference.pos();
            let atom_position_vector = Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32);
            minimum_coordinate_bound = minimum_coordinate_bound.min(atom_position_vector);
            maximum_coordinate_bound = maximum_coordinate_bound.max(atom_position_vector);
        }

        (minimum_coordinate_bound, maximum_coordinate_bound)
    }

    pub fn ca_positions(&self) -> Vec<(Vec3, String)> {
        let mut alpha_carbon_positions_collection = Vec::new();

        for current_chain_reference in self.pdb.chains() {
            let chain_id_string = current_chain_reference.id().to_string();
            for current_residue_reference in current_chain_reference.residues() {
                for current_conformer_reference in current_residue_reference.conformers() {
                    for current_atom_reference in current_conformer_reference.atoms() {
                        if current_atom_reference.name() == "CA" {
                            let atom_coordinates_tuple = current_atom_reference.pos();
                            alpha_carbon_positions_collection.push((
                                Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32),
                                chain_id_string.clone(),
                            ));
                        }
                    }
                }
            }
        }

        alpha_carbon_positions_collection
    }

    /// Returns backbone segments as pairs of (start, end, chain_id) for line rendering
    pub fn backbone_segments(&self) -> Vec<(Vec3, Vec3, String)> {
        let mut backbone_segment_collection = Vec::new();

        for current_chain_reference in self.pdb.chains() {
            let chain_id_string = current_chain_reference.id().to_string();
            let mut previous_alpha_carbon_position: Option<Vec3> = None;

            for current_residue_reference in current_chain_reference.residues() {
                for current_conformer_reference in current_residue_reference.conformers() {
                    for current_atom_reference in current_conformer_reference.atoms() {
                        if current_atom_reference.name() == "CA" {
                            let atom_coordinates_tuple = current_atom_reference.pos();
                            let current_alpha_carbon_position = Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32);

                            if let Some(previous_position_reference) = previous_alpha_carbon_position {
                                // Only connect if distance is reasonable (< 5 Angstroms)
                                if previous_position_reference.distance(current_alpha_carbon_position) < 5.0 {
                                    backbone_segment_collection.push((previous_position_reference, current_alpha_carbon_position, chain_id_string.clone()));
                                }
                            }
                            previous_alpha_carbon_position = Some(current_alpha_carbon_position);
                        }
                    }
                }
            }
        }

        backbone_segment_collection
    }

    /// Get B-factor range for coloring
    pub fn bfactor_range(&self) -> (f32, f32) {
        let mut minimum_bfactor_value = f32::MAX;
        let mut maximum_bfactor_value = f32::MIN;

        for current_atom_reference in self.pdb.atoms() {
            let atom_bfactor_value = current_atom_reference.b_factor() as f32;
            minimum_bfactor_value = minimum_bfactor_value.min(atom_bfactor_value);
            maximum_bfactor_value = maximum_bfactor_value.max(atom_bfactor_value);
        }

        (minimum_bfactor_value, maximum_bfactor_value)
    }

    /// Get CA atoms with their B-factors
    pub fn ca_with_bfactor(&self) -> Vec<(Vec3, String, f32)> {
        let mut alpha_carbon_bfactor_collection = Vec::new();

        for current_chain_reference in self.pdb.chains() {
            let chain_id_string = current_chain_reference.id().to_string();
            for current_residue_reference in current_chain_reference.residues() {
                for current_conformer_reference in current_residue_reference.conformers() {
                    for current_atom_reference in current_conformer_reference.atoms() {
                        if current_atom_reference.name() == "CA" {
                            let atom_coordinates_tuple = current_atom_reference.pos();
                            alpha_carbon_bfactor_collection.push((
                                Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32),
                                chain_id_string.clone(),
                                current_atom_reference.b_factor() as f32,
                            ));
                        }
                    }
                }
            }
        }

        alpha_carbon_bfactor_collection
    }
}

pub struct ProteinStore {
    proteins: HashMap<String, Arc<RwLock<ProteinData>>>,
}

impl ProteinStore {
    pub fn new() -> Self {
        Self {
            proteins: HashMap::new(),
        }
    }

    pub fn fetch(&mut self, pdb_identifier_code: &str) -> Result<Arc<RwLock<ProteinData>>, String> {
        let pdb_identifier_code = pdb_identifier_code.to_uppercase();

        if let Some(previously_fetched_protein) = self.proteins.get(&pdb_identifier_code) {
            return Ok(previously_fetched_protein.clone());
        }

        let fetch_operation_result = fetch_pdb(&pdb_identifier_code)?;
        let parsed_protein_data = ProteinData::from_string(&fetch_operation_result.content, &pdb_identifier_code, fetch_operation_result.format)?;
        let shared_protein_handle = Arc::new(RwLock::new(parsed_protein_data));
        self.proteins.insert(pdb_identifier_code, shared_protein_handle.clone());
        Ok(shared_protein_handle)
    }

    pub fn load(&mut self, file_system_path: &str) -> Result<Arc<RwLock<ProteinData>>, String> {
        let extracted_protein_name = std::path::Path::new(file_system_path)
            .file_stem()
            .and_then(|file_stem_string| file_stem_string.to_str())
            .unwrap_or("unknown")
            .to_string();

        if let Some(previously_loaded_protein) = self.proteins.get(&extracted_protein_name) {
            return Ok(previously_loaded_protein.clone());
        }

        let load_operation_result = load_file(file_system_path)?;
        let parsed_protein_data = ProteinData::from_string(&load_operation_result.content, &extracted_protein_name, load_operation_result.format)?;
        let shared_protein_handle = Arc::new(RwLock::new(parsed_protein_data));
        self.proteins.insert(extracted_protein_name, shared_protein_handle.clone());
        Ok(shared_protein_handle)
    }

    pub fn list(&self) -> Vec<String> {
        self.proteins.keys().cloned().collect()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Arc<RwLock<ProteinData>>> {
        self.proteins.values()
    }
}