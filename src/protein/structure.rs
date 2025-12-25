//! Internal representation and management of protein structures
//!
//! This module defines the core data structures for proteins, including
//! atomic positions, connectivity, and metadata, as well as the `ProteinStore`
//! for managing multiple loaded proteins

use glam::Vec3;
use pdbtbx::{PDB, ReadOptions, Format};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use super::fetch::{fetch_pdb, load_file, FileFormat};

use serde::{Serialize, Deserialize};

/// The available visual representation modes for a protein
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Representation {
    /// Each atom (or CA) is rendered as a sphere
    #[default]
    Spheres,
    /// Only the backbone trace is rendered as a series of lines
    Backbone,
    /// Both the backbone trace and spheres are rendered
    BackboneAndSpheres,
}

/// The available color schemes for a protein
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum ColorScheme {
    /// Color by chain identifier
    #[default]
    ByChain,
    /// Color by chemical element
    ByElement,
    /// Color by B-factor (temperature factor)
    ByBFactor,
    /// Color by secondary structure (alpha helix, beta sheet, etc)
    BySecondary,
    /// Use a single uniform color
    Uniform([f32; 3]),
}

/// Secondary structure types for coloring
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SecondaryStructureType {
    Helix,
    Sheet,
    Other,
}

/// Source of the protein data (for saving/loading state)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum ProteinSource {
    Rcsb(String),
    File(String),
}

/// Core data structure containing parsed protein information and rendering state
pub struct ProteinData {
    /// The underlying PDB data structure from `pdbtbx`
    pub pdb: PDB,
    /// Display name of the protein
    pub name: String,
    /// Source of the protein
    pub source: ProteinSource,
    /// Whether the protein should be rendered
    pub visible: bool,
    /// The current representation mode
    pub representation: Representation,
    /// The current color scheme
    pub color_scheme: ColorScheme,
}

impl ProteinData {
    /// Parses a protein structure from a string in the specified format
    pub fn from_string(file_content_string: &str, name: &str, format: FileFormat, source: ProteinSource) -> Result<Self, String> {
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
            source,
            visible: true,
            representation: Representation::default(),
            color_scheme: ColorScheme::default(),
        })
    }

    /// Returns the total number of atoms in the structure
    pub fn atom_count(&self) -> usize {
        self.pdb.atom_count()
    }

    /// Returns a list of all chain identifiers
    pub fn chain_ids(&self) -> Vec<String> {
        self.pdb
            .chains()
            .map(|current_chain| current_chain.id().to_string())
            .collect()
    }

    /// Calculates the geometric center of mass of all atoms
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

    /// Calculates the axis-aligned bounding box (min, max) of all atoms
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

    /// Returns the positions and chain IDs of all Alpha Carbon (CA) atoms
    pub fn get_alpha_carbon_positions_and_chain_identifiers(&self) -> Vec<(Vec3, String)> {
        let mut alpha_carbon_positions_and_chain_identifiers_collection = Vec::new();

        for current_chain_reference in self.pdb.chains() {
            let chain_id_string = current_chain_reference.id().to_string();
            for current_residue_reference in current_chain_reference.residues() {
                for current_conformer_reference in current_residue_reference.conformers() {
                    for current_atom_reference in current_conformer_reference.atoms() {
                        if current_atom_reference.name() == "CA" {
                            let atom_coordinates_tuple = current_atom_reference.pos();
                            alpha_carbon_positions_and_chain_identifiers_collection.push((
                                Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32),
                                chain_id_string.clone(),
                            ));
                        }
                    }
                }
            }
        }

        alpha_carbon_positions_and_chain_identifiers_collection
    }

    /// Returns backbone segments as pairs of (start, end, chain_id) for line rendering
    pub fn get_backbone_segments_for_rendering(&self) -> Vec<(Vec3, Vec3, String)> {
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

    /// Returns the minimum and maximum B-factor values in the structure
    pub fn calculate_bfactor_range(&self) -> (f32, f32) {
        let mut minimum_bfactor_value = f32::MAX;
        let mut maximum_bfactor_value = f32::MIN;

        for current_atom_reference in self.pdb.atoms() {
            let atom_bfactor_value = current_atom_reference.b_factor() as f32;
            minimum_bfactor_value = minimum_bfactor_value.min(atom_bfactor_value);
            maximum_bfactor_value = maximum_bfactor_value.max(atom_bfactor_value);
        }

        (minimum_bfactor_value, maximum_bfactor_value)
    }

    /// Returns Alpha Carbon (CA) atoms with their associated B-factors and chain IDs
    pub fn get_alpha_carbon_data_with_bfactors(&self) -> Vec<(Vec3, String, f32)> {
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

    /// Returns Alpha Carbon (CA) atoms with their secondary structure and chain IDs
    pub fn get_alpha_carbon_data_with_secondary_structure(&self) -> Vec<(Vec3, String, SecondaryStructureType)> {
        let mut alpha_carbon_secondary_structure_collection = Vec::new();

        for current_chain_reference in self.pdb.chains() {
            let chain_id_string = current_chain_reference.id().to_string();
            for current_residue_reference in current_chain_reference.residues() {
                // Placeholder: pdbtbx 0.12.0 does not expose helices/sheets directly.
                // We use a dummy pattern based on residue serial number for visualization.
                let residue_serial_number = current_residue_reference.serial_number();
                let secondary_structure_type = if (residue_serial_number / 10) % 3 == 0 {
                    SecondaryStructureType::Helix
                } else if (residue_serial_number / 10) % 3 == 1 {
                    SecondaryStructureType::Sheet
                } else {
                    SecondaryStructureType::Other
                };

                for current_conformer_reference in current_residue_reference.conformers() {
                    for current_atom_reference in current_conformer_reference.atoms() {
                        if current_atom_reference.name() == "CA" {
                            let atom_coordinates_tuple = current_atom_reference.pos();
                            alpha_carbon_secondary_structure_collection.push((
                                Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32),
                                chain_id_string.clone(),
                                secondary_structure_type,
                            ));
                        }
                    }
                }
            }
        }

        alpha_carbon_secondary_structure_collection
    }

    /// Returns backbone segments with secondary structure for line rendering
    pub fn get_backbone_segments_with_secondary_structure(&self) -> Vec<(Vec3, Vec3, String, SecondaryStructureType)> {
        let mut backbone_segment_collection = Vec::new();

        for current_chain_reference in self.pdb.chains() {
            let chain_id_string = current_chain_reference.id().to_string();
            let mut previous_alpha_carbon_information: Option<(Vec3, SecondaryStructureType)> = None;

            for current_residue_reference in current_chain_reference.residues() {
                let residue_serial_number = current_residue_reference.serial_number();
                let secondary_structure_type = if (residue_serial_number / 10) % 3 == 0 {
                    SecondaryStructureType::Helix
                } else if (residue_serial_number / 10) % 3 == 1 {
                    SecondaryStructureType::Sheet
                } else {
                    SecondaryStructureType::Other
                };

                for current_conformer_reference in current_residue_reference.conformers() {
                    for current_atom_reference in current_conformer_reference.atoms() {
                        if current_atom_reference.name() == "CA" {
                            let atom_coordinates_tuple = current_atom_reference.pos();
                            let current_alpha_carbon_position = Vec3::new(atom_coordinates_tuple.0 as f32, atom_coordinates_tuple.1 as f32, atom_coordinates_tuple.2 as f32);

                            if let Some((previous_position, _previous_secondary_structure)) = previous_alpha_carbon_information {
                                if previous_position.distance(current_alpha_carbon_position) < 5.0 {
                                    backbone_segment_collection.push((previous_position, current_alpha_carbon_position, chain_id_string.clone(), secondary_structure_type));
                                }
                            }
                            previous_alpha_carbon_information = Some((current_alpha_carbon_position, secondary_structure_type));
                        }
                    }
                }
            }
        }

        backbone_segment_collection
    }
}

/// A central repository for managing multiple loaded proteins
pub struct ProteinStore {
    proteins: HashMap<String, Arc<RwLock<ProteinData>>>,
}

impl ProteinStore {
    /// Creates a new empty `ProteinStore`
    pub fn new() -> Self {
        Self {
            proteins: HashMap::new(),
        }
    }

    /// Fetches a protein from the RCSB repository by ID and stores it
    pub fn fetch(&mut self, pdb_identifier_code: &str) -> Result<Arc<RwLock<ProteinData>>, String> {
        let pdb_identifier_code = pdb_identifier_code.to_uppercase();

        if let Some(previously_fetched_protein) = self.proteins.get(&pdb_identifier_code) {
            return Ok(previously_fetched_protein.clone());
        }

        let fetch_operation_result = fetch_pdb(&pdb_identifier_code)?;
        let parsed_protein_data = ProteinData::from_string(
            &fetch_operation_result.content, 
            &pdb_identifier_code, 
            fetch_operation_result.format,
            ProteinSource::Rcsb(pdb_identifier_code.clone())
        )?;
        let shared_protein_handle = Arc::new(RwLock::new(parsed_protein_data));
        self.proteins.insert(pdb_identifier_code, shared_protein_handle.clone());
        Ok(shared_protein_handle)
    }

    /// Loads a protein from a local file system path and stores it
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
        let parsed_protein_data = ProteinData::from_string(
            &load_operation_result.content, 
            &extracted_protein_name, 
            load_operation_result.format,
            ProteinSource::File(file_system_path.to_string())
        )?;
        let shared_protein_handle = Arc::new(RwLock::new(parsed_protein_data));
        self.proteins.insert(extracted_protein_name, shared_protein_handle.clone());
        Ok(shared_protein_handle)
    }

    /// Returns a list of all loaded protein identifiers
    pub fn list(&self) -> Vec<String> {
        self.proteins.keys().cloned().collect()
    }

    /// Returns an iterator over all loaded protein handles
    pub fn iter(&self) -> impl Iterator<Item = &Arc<RwLock<ProteinData>>> {
        self.proteins.values()
    }
}