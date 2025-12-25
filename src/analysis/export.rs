//! Export of protein analysis data to JSON and CSV formats

use serde::{Serialize, Deserialize};
use crate::protein::structure::ProteinData;
use std::fs::File;
use std::io::Write;

/// A structured report of protein properties and analysis results
#[derive(Serialize, Deserialize, Debug)]
pub struct ProteinAnalysisReport {
    pub protein_name: String,
    pub total_atom_count: usize,
    pub chain_identifiers: Vec<String>,
    pub center_of_mass_coordinates: [f32; 3],
    pub bounding_box_dimensions: [f32; 3],
    pub residues: Vec<ResidueReportEntry>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct ResidueReportEntry {
    pub chain_id: String,
    pub residue_number: isize,
    pub residue_name: String,
    pub phi_angle: Option<f32>,
    pub psi_angle: Option<f32>,
    pub solvent_accessible_surface_area: Option<f32>,
}

/// Exports protein analysis data to a JSON file
pub fn export_protein_analysis_to_json(
    protein_data: &ProteinData,
    file_system_path: &str,
) -> Result<(), String> {
    let (min_bound, max_bound) = protein_data.bounding_box();
    let bbox_size = max_bound - min_bound;
    let com = protein_data.center_of_mass();

    let dihedrals = crate::analysis::dihedrals::calculate_all_backbone_dihedrals(protein_data);
    
    let mut residue_entries = Vec::new();
    for (chain_id, res_num, angles) in dihedrals {
        residue_entries.push(ResidueReportEntry {
            chain_id,
            residue_number: res_num,
            residue_name: "UNKNOWN".to_string(), // Would need lookup to be more precise
            phi_angle: angles.phi_angle,
            psi_angle: angles.psi_angle,
            solvent_accessible_surface_area: None,
        });
    }

    let report_data = ProteinAnalysisReport {
        protein_name: protein_data.name.clone(),
        total_atom_count: protein_data.atom_count(),
        chain_identifiers: protein_data.chain_ids(),
        center_of_mass_coordinates: [com.x, com.y, com.z],
        bounding_box_dimensions: [bbox_size.x, bbox_size.y, bbox_size.z],
        residues: residue_entries,
    };

    let serialized_json_string = serde_json::to_string_pretty(&report_data)
        .map_err(|json_error| format!("JSON serialization error: {}", json_error))?;

    let mut output_file_handle = File::create(file_system_path)
        .map_err(|io_error| format!("Failed to create file: {}", io_error))?;
    
    output_file_handle.write_all(serialized_json_string.as_bytes())
        .map_err(|io_error| format!("Failed to write to file: {}", io_error))?;

    Ok(())
}

/// Exports residue analysis data to a CSV file
pub fn export_residue_data_to_csv(
    protein_data: &ProteinData,
    file_system_path: &str,
) -> Result<(), String> {
    let dihedrals = crate::analysis::dihedrals::calculate_all_backbone_dihedrals(protein_data);
    
    let mut output_file_handle = File::create(file_system_path)
        .map_err(|io_error| format!("Failed to create file: {}", io_error))?;
    
    writeln!(output_file_handle, "Chain,ResidueNumber,Phi,Psi")
        .map_err(|io_error| format!("Failed to write CSV header: {}", io_error))?;

    for (chain_id, res_num, angles) in dihedrals {
        let phi_string = angles.phi_angle.map(|a| format!("{:.2}", a)).unwrap_or_else(|| "".to_string());
        let psi_string = angles.psi_angle.map(|a| format!("{:.2}", a)).unwrap_or_else(|| "".to_string());
        
        writeln!(output_file_handle, "{},{},{},{}", chain_id, res_num, phi_string, psi_string)
            .map_err(|io_error| format!("Failed to write CSV row: {}", io_error))?;
    }

    Ok(())
}
