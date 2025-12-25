//! Internal representation and management of protein structures
//!
//! This module defines the core data structures for proteins, including
//! atomic positions, connectivity, and metadata, as well as the `ProteinStore`
//! for managing multiple loaded proteins

use glam::Vec3;
use pdbtbx::{Format, ReadOptions, PDB, ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use tracing::warn;

use super::fetch_rcsb::{fetch_pdb, load_file, FileFormat};
use crate::selection::SelectionSet;

use serde::{Deserialize, Serialize};

/// The available visual representation modes for a protein
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Representation {
    /// Each atom (or CA) is rendered as a sphere
    Spheres,
    /// Only the backbone trace is rendered as a series of lines
    Backbone,
    /// Both the backbone trace and spheres are rendered
    BackboneAndSpheres,
    /// Cylinders for all covalent bonds
    Sticks,
    /// Spheres for atoms and cylinders for bonds
    #[default]
    BallAndStick,
    /// Van der Waals spheres for all atoms
    SpaceFilling,
    /// Simple lines for all covalent bonds
    Lines,
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
    /// The underlying PDB data structure from `pdbtbx`, if loaded
    pub underlying_pdb_data: Option<PDB>,
    /// Display name of the protein
    pub display_name: String,
    /// Source of the protein
    pub data_source: ProteinSource,
    /// Whether the protein should be rendered
    pub is_visible: bool,
    /// Whether the protein is currently being fetched or parsed
    pub is_currently_loading: bool,
    /// Error message if loading failed
    pub loading_error_message: Option<String>,
    /// The current representation mode
    pub visual_representation: Representation,
    /// The current color scheme
    pub active_color_scheme: ColorScheme,
    /// Cached set of covalent bonds identified in the structure
    pub identified_atom_bonds: std::collections::HashSet<super::bonds::AtomBond>,
    /// The generated molecular surface mesh for this protein
    pub molecular_surface_mesh: crate::surface::MolecularSurfaceMesh,
    /// A counter incremented whenever the protein data or its visual state changes
    pub structural_data_revision_number: usize,
}

impl ProteinData {
    /// Creates a new `ProteinData` in a pending loading state
    pub fn new_pending(display_name: &str, source: ProteinSource) -> Self {
        Self {
            underlying_pdb_data: None,
            display_name: display_name.to_string(),
            data_source: source,
            is_visible: true,
            is_currently_loading: true,
            loading_error_message: None,
            visual_representation: Representation::default(),
            active_color_scheme: ColorScheme::default(),
            identified_atom_bonds: std::collections::HashSet::new(),
            molecular_surface_mesh: crate::surface::MolecularSurfaceMesh::new(),
            structural_data_revision_number: 0,
        }
    }

    /// Parses a protein structure from a string in the specified format
    pub fn from_string(
        file_content_string: &str,
        display_name: &str,
        format: FileFormat,
        source: ProteinSource,
    ) -> Result<Self, String> {
        let target_file_format = match format {
            FileFormat::Pdb => Format::Pdb,
            FileFormat::Cif => Format::Mmcif,
        };

        let (pdb, parsing_warning_messages) = ReadOptions::default()
            .set_level(pdbtbx::StrictnessLevel::Loose)
            .set_format(target_file_format)
            .read_raw(std::io::BufReader::new(std::io::Cursor::new(
                file_content_string.as_bytes(),
            )))
            .map_err(|pdb_parse_error| format!("Parse error: {:?}", pdb_parse_error))?;

        if !parsing_warning_messages.is_empty() {
            warn!(
                "Warnings while parsing {}: {:?}",
                display_name, parsing_warning_messages
            );
        }

        let bond_calculator = super::bonds::BondCalculator::new(&pdb);
        let identified_atom_bonds = bond_calculator.calculate_all_bonds();

        Ok(Self {
            underlying_pdb_data: Some(pdb),
            display_name: display_name.to_string(),
            data_source: source,
            is_visible: true,
            is_currently_loading: false,
            loading_error_message: None,
            visual_representation: Representation::default(),
            active_color_scheme: ColorScheme::default(),
            identified_atom_bonds,
            molecular_surface_mesh: crate::surface::MolecularSurfaceMesh::new(),
            structural_data_revision_number: 1,
        })
    }

    /// Generates a molecular surface mesh for the protein structure
    pub fn compute_molecular_surface_mesh(
        &mut self,
        voxel_resolution_angstroms: f32,
        isosurface_threshold_value: f32,
        edge_intersection_lookup_table: &[u32; 256],
        triangulation_lookup_table: &[[i32; 16]; 256],
    ) {
        let distance_field_grid = crate::surface::distance_field::calculate_solvent_accessible_distance_field(
            self,
            voxel_resolution_angstroms,
            4.0, // Grid padding in angstroms
        );

        let (surface_vertices, surface_indices) = crate::surface::marching_cubes::extract_isosurface_from_distance_field(
            &distance_field_grid,
            isosurface_threshold_value,
            edge_intersection_lookup_table,
            triangulation_lookup_table,
        );

        self.molecular_surface_mesh.surface_vertices_collection = surface_vertices;
        self.molecular_surface_mesh.surface_indices_collection = surface_indices;
        self.molecular_surface_mesh.is_surface_visible = true;
        self.structural_data_revision_number += 1;
    }

    /// Returns the total number of atoms in the structure
    pub fn get_total_atom_count(&self) -> usize {
        self.underlying_pdb_data.as_ref().map(|pdb| pdb.atom_count()).unwrap_or(0)
    }

    /// Returns a list of all chain identifiers
    pub fn get_all_chain_identifiers(&self) -> Vec<String> {
        self.underlying_pdb_data.as_ref()
            .map(|pdb| pdb.chains().map(|current_chain| current_chain.id().to_string()).collect())
            .unwrap_or_default()
    }

    /// Selects atoms based on a PyMOL-style selection string
    pub fn select_atoms_from_string(&self, selection_query_string: &str) -> Result<SelectionSet, String> {
        let selection_expression = crate::selection::parser::parse_selection_expression(selection_query_string)?;
        let selection_evaluator = crate::selection::evaluator::Evaluator::new(self);
        Ok(selection_evaluator.evaluate_expression(&selection_expression))
    }

    /// Calculates the geometric center of mass of all atoms
    pub fn calculate_geometric_center_of_mass(&self) -> Vec3 {
        let mut position_sum_vector = Vec3::ZERO;
        let mut atom_counting_index = 0;

        if let Some(pdb) = &self.underlying_pdb_data {
            for current_atom_reference in pdb.atoms() {
                let atom_coordinates_tuple = current_atom_reference.pos();
                position_sum_vector += Vec3::new(
                    atom_coordinates_tuple.0 as f32,
                    atom_coordinates_tuple.1 as f32,
                    atom_coordinates_tuple.2 as f32,
                );
                atom_counting_index += 1;
            }
        }

        if atom_counting_index > 0 {
            position_sum_vector / atom_counting_index as f32
        } else {
            Vec3::ZERO
        }
    }

    /// Calculates the axis-aligned bounding box (min, max) of all atoms
    pub fn calculate_axis_aligned_bounding_box(&self) -> (Vec3, Vec3) {
        let mut minimum_coordinate_bound = Vec3::splat(f32::MAX);
        let mut maximum_coordinate_bound = Vec3::splat(f32::MIN);

        if let Some(pdb) = &self.underlying_pdb_data {
            for current_atom_reference in pdb.atoms() {
                let atom_coordinates_tuple = current_atom_reference.pos();
                let atom_position_vector = Vec3::new(
                    atom_coordinates_tuple.0 as f32,
                    atom_coordinates_tuple.1 as f32,
                    atom_coordinates_tuple.2 as f32,
                );
                minimum_coordinate_bound = minimum_coordinate_bound.min(atom_position_vector);
                maximum_coordinate_bound = maximum_coordinate_bound.max(atom_position_vector);
            }
        } else {
            return (Vec3::ZERO, Vec3::ZERO);
        }

        (minimum_coordinate_bound, maximum_coordinate_bound)
    }

    /// Returns the positions of all Alpha Carbon (CA) atoms in the first model only
    pub fn get_alpha_carbon_positions_for_first_model(&self) -> Vec<Vec3> {
        use rayon::prelude::*;

        if let Some(pdb) = &self.underlying_pdb_data {
            if let Some(first_model_reference) = pdb.models().next() {
                // Use par_bridge to parallelize the hierarchy traversal
                first_model_reference.chains()
                    .flat_map(|current_chain_reference| current_chain_reference.residues())
                    .flat_map(|current_residue_reference| current_residue_reference.conformers())
                    .flat_map(|current_conformer_reference| current_conformer_reference.atoms())
                    .par_bridge() // Parallelize the processing of collected atoms
                    .filter(|current_atom_reference| current_atom_reference.name().trim() == "CA")
                    .map(|current_atom_reference| {
                        let atom_coordinates_tuple = current_atom_reference.pos();
                        Vec3::new(
                            atom_coordinates_tuple.0 as f32,
                            atom_coordinates_tuple.1 as f32,
                            atom_coordinates_tuple.2 as f32,
                        )
                    })
                    .collect()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        }
    }

    /// Returns the positions and chain IDs of all Alpha Carbon (CA) atoms
    pub fn get_alpha_carbon_positions_and_chain_identifiers(&self) -> Vec<(Vec3, String)> {
        use rayon::prelude::*;

        if let Some(pdb) = &self.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .par_bridge()
                .filter(|current_hierarchy_reference| current_hierarchy_reference.atom().name().trim() == "CA")
                .map(|current_hierarchy_reference| {
                    let atom_coordinates_tuple = current_hierarchy_reference.atom().pos();
                    (
                        Vec3::new(
                            atom_coordinates_tuple.0 as f32,
                            atom_coordinates_tuple.1 as f32,
                            atom_coordinates_tuple.2 as f32,
                        ),
                        current_hierarchy_reference.chain().id().to_string(),
                    )
                })
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Returns backbone segments as pairs of (start, end, chain_id) for line rendering
    pub fn get_backbone_segments_for_rendering(&self) -> Vec<(Vec3, Vec3, String)> {
        let mut backbone_segment_collection = Vec::new();

        if let Some(pdb) = &self.underlying_pdb_data {
            for current_chain_reference in pdb.chains() {
                let chain_id_string = current_chain_reference.id().to_string();
                let mut previous_alpha_carbon_position: Option<Vec3> = None;

                for current_residue_reference in current_chain_reference.residues() {
                    for current_conformer_reference in current_residue_reference.conformers() {
                        for current_atom_reference in current_conformer_reference.atoms() {
                            if current_atom_reference.name().trim() == "CA" {
                                let atom_coordinates_tuple = current_atom_reference.pos();
                                let current_alpha_carbon_position = Vec3::new(
                                    atom_coordinates_tuple.0 as f32,
                                    atom_coordinates_tuple.1 as f32,
                                    atom_coordinates_tuple.2 as f32,
                                );

                                if let Some(previous_position_reference) =
                                    previous_alpha_carbon_position
                                {
                                    // Only connect if distance is reasonable (< 5 Angstroms)
                                    if previous_position_reference
                                        .distance(current_alpha_carbon_position)
                                        < 5.0
                                    {
                                        backbone_segment_collection.push((
                                            previous_position_reference,
                                            current_alpha_carbon_position,
                                            chain_id_string.clone(),
                                        ));
                                    }
                                }
                                previous_alpha_carbon_position = Some(current_alpha_carbon_position);
                            }
                        }
                    }
                }
            }
        }

        backbone_segment_collection
    }

    /// Returns the minimum and maximum B-factor values in the structure
    pub fn calculate_bfactor_range_in_structure(&self) -> (f32, f32) {
        use rayon::prelude::*;

        if let Some(pdb) = &self.underlying_pdb_data {
            let bfactor_values_collection: Vec<f32> = pdb.atoms()
                .par_bridge()
                .map(|current_atom_reference| current_atom_reference.b_factor() as f32)
                .collect();

            if bfactor_values_collection.is_empty() {
                return (0.0, 0.0);
            }

            let minimum_bfactor_value = *bfactor_values_collection.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            let maximum_bfactor_value = *bfactor_values_collection.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

            (minimum_bfactor_value, maximum_bfactor_value)
        } else {
            (0.0, 0.0)
        }
    }

    /// Returns Alpha Carbon (CA) atoms with their associated B-factors and chain IDs
    pub fn get_alpha_carbon_data_with_bfactors(&self) -> Vec<(Vec3, String, f32)> {
        use rayon::prelude::*;

        if let Some(pdb) = &self.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .par_bridge()
                .filter(|current_hierarchy_reference| current_hierarchy_reference.atom().name().trim() == "CA")
                .map(|current_hierarchy_reference| {
                    let atom_coordinates_tuple = current_hierarchy_reference.atom().pos();
                    (
                        Vec3::new(
                            atom_coordinates_tuple.0 as f32,
                            atom_coordinates_tuple.1 as f32,
                            atom_coordinates_tuple.2 as f32,
                        ),
                        current_hierarchy_reference.chain().id().to_string(),
                        current_hierarchy_reference.atom().b_factor() as f32,
                    )
                })
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Returns Alpha Carbon (CA) atoms with their secondary structure and chain IDs
    pub fn get_alpha_carbon_data_with_secondary_structure(
        &self,
    ) -> Vec<(Vec3, String, SecondaryStructureType)> {
        use rayon::prelude::*;

        if let Some(pdb) = &self.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .par_bridge()
                .filter(|current_hierarchy_reference| current_hierarchy_reference.atom().name().trim() == "CA")
                .map(|current_hierarchy_reference| {
                    let chain_id_string = current_hierarchy_reference.chain().id();
                    let residue_number = current_hierarchy_reference.residue().serial_number();
                    let atom_coordinates_tuple = current_hierarchy_reference.atom().pos();

                    let secondary_structure_type = if pdb.is_residue_in_helix(chain_id_string, residue_number) {
                        SecondaryStructureType::Helix
                    } else if pdb.is_residue_in_sheet(chain_id_string, residue_number) {
                        SecondaryStructureType::Sheet
                    } else {
                        SecondaryStructureType::Other
                    };

                    (
                        Vec3::new(
                            atom_coordinates_tuple.0 as f32,
                            atom_coordinates_tuple.1 as f32,
                            atom_coordinates_tuple.2 as f32,
                        ),
                        chain_id_string.to_string(),
                        secondary_structure_type,
                    )
                })
                .collect()
        } else {
            Vec::new()
        }
    }

    /// Returns backbone segments with secondary structure for line rendering
    pub fn get_backbone_segments_with_secondary_structure(
        &self,
    ) -> Vec<(Vec3, Vec3, String, SecondaryStructureType)> {
        let mut backbone_segment_collection = Vec::new();

        if let Some(pdb) = &self.underlying_pdb_data {
            for current_chain_reference in pdb.chains() {
                let chain_id_string = current_chain_reference.id().to_string();
                let mut previous_alpha_carbon_information: Option<(Vec3, SecondaryStructureType)> =
                    None;

                for current_residue_reference in current_chain_reference.residues() {
                    let residue_number = current_residue_reference.serial_number();
                    let secondary_structure_type = if pdb.is_residue_in_helix(&chain_id_string, residue_number) {
                        SecondaryStructureType::Helix
                    } else if pdb.is_residue_in_sheet(&chain_id_string, residue_number) {
                        SecondaryStructureType::Sheet
                    } else {
                        SecondaryStructureType::Other
                    };

                    for current_conformer_reference in current_residue_reference.conformers() {
                        for current_atom_reference in current_conformer_reference.atoms() {
                            if current_atom_reference.name().trim() == "CA" {
                                let atom_coordinates_tuple = current_atom_reference.pos();
                                let current_alpha_carbon_position = Vec3::new(
                                    atom_coordinates_tuple.0 as f32,
                                    atom_coordinates_tuple.1 as f32,
                                    atom_coordinates_tuple.2 as f32,
                                );

                                if let Some((previous_position, _previous_secondary_structure)) =
                                    previous_alpha_carbon_information
                                {
                                    if previous_position.distance(current_alpha_carbon_position) < 5.0 {
                                        // Use the secondary structure of the current residue for the segment
                                        // Or we could interpolate, but usually segments are assigned to the latter residue
                                        backbone_segment_collection.push((
                                            previous_position,
                                            current_alpha_carbon_position,
                                            chain_id_string.clone(),
                                            secondary_structure_type,
                                        ));
                                    }
                                }
                                previous_alpha_carbon_information =
                                    Some((current_alpha_carbon_position, secondary_structure_type));
                            }
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

    /// Fetches a protein from the RCSB repository by ID and stores it (Non-blocking)
    pub fn fetch(&mut self, pdb_identifier_code: &str) -> Result<Arc<RwLock<ProteinData>>, String> {
        let pdb_identifier_code_uppercase = pdb_identifier_code.to_uppercase();

        if let Some(previously_fetched_protein) = self.proteins.get(&pdb_identifier_code_uppercase) {
            return Ok(previously_fetched_protein.clone());
        }

        let pending_protein_data = ProteinData::new_pending(&pdb_identifier_code_uppercase, ProteinSource::Rcsb(pdb_identifier_code_uppercase.clone()));
        let shared_protein_handle = Arc::new(RwLock::new(pending_protein_data));
        self.proteins.insert(pdb_identifier_code_uppercase.clone(), shared_protein_handle.clone());

        // Spawn a background thread to handle network fetch and parsing
        let handle_clone = shared_protein_handle.clone();
        let code_clone = pdb_identifier_code_uppercase.clone();
        std::thread::spawn(move || {
            let fetch_result = fetch_pdb(&code_clone);
            let mut protein_data = handle_clone.write().unwrap();
            match fetch_result {
                Ok(result) => {
                    match ProteinData::from_string(&result.content, &code_clone.clone(), result.format, ProteinSource::Rcsb(code_clone)) {
                        Ok(new_data) => {
                            *protein_data = new_data;
                        }
                        Err(parse_error) => {
                            protein_data.is_currently_loading = false;
                            protein_data.loading_error_message = Some(parse_error);
                        }
                    }
                }
                Err(fetch_error) => {
                    protein_data.is_currently_loading = false;
                    protein_data.loading_error_message = Some(fetch_error);
                }
            }
            protein_data.structural_data_revision_number += 1;
        });

        Ok(shared_protein_handle)
    }

    /// Loads a protein from a local file system path and stores it (Non-blocking)
    pub fn load(&mut self, file_system_path: &str) -> Result<Arc<RwLock<ProteinData>>, String> {
        let extracted_protein_name = std::path::Path::new(file_system_path)
            .file_stem()
            .and_then(|file_stem_string| file_stem_string.to_str())
            .unwrap_or("unknown")
            .to_string();

        if let Some(previously_loaded_protein) = self.proteins.get(&extracted_protein_name) {
            return Ok(previously_loaded_protein.clone());
        }

        let pending_protein_data = ProteinData::new_pending(&extracted_protein_name, ProteinSource::File(file_system_path.to_string()));
        let shared_protein_handle = Arc::new(RwLock::new(pending_protein_data));
        self.proteins.insert(extracted_protein_name.clone(), shared_protein_handle.clone());

        let handle_clone = shared_protein_handle.clone();
        let path_clone = file_system_path.to_string();
        let name_clone = extracted_protein_name.clone();
        std::thread::spawn(move || {
            let load_result = load_file(&path_clone);
            let mut protein_data = handle_clone.write().unwrap();
            match load_result {
                Ok(result) => {
                    match ProteinData::from_string(&result.content, &name_clone, result.format, ProteinSource::File(path_clone)) {
                        Ok(new_data) => {
                            *protein_data = new_data;
                        }
                        Err(parse_error) => {
                            protein_data.is_currently_loading = false;
                            protein_data.loading_error_message = Some(parse_error);
                        }
                    }
                }
                Err(load_error) => {
                    protein_data.is_currently_loading = false;
                    protein_data.loading_error_message = Some(load_error);
                }
            }
            protein_data.structural_data_revision_number += 1;
        });

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

    /// Clears all proteins from the store
    pub fn clear(&mut self) {
        self.proteins.clear();
    }

    /// Adds an existing protein to the store
    pub fn add(&mut self, protein: Arc<RwLock<ProteinData>>) {
        let display_name = protein.read().unwrap().display_name.clone();
        self.proteins.insert(display_name, protein);
    }
}
