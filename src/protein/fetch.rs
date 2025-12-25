//! Utilities for fetching protein data from remote and local sources
//!
//! This module handles downloading PDB/mmCIF files from RCSB and reading
//! them from the local file system

const RCSB_BASE_URL: &str = "https://files.rcsb.org/download";

/// Supported protein file formats
#[derive(Clone, Copy, Debug)]
pub enum FileFormat {
    /// Protein Data Bank format (.pdb)
    Pdb,
    /// Macromolecular Crystallographic Information File format (.cif)
    Cif,
}

/// The result of a fetch or load operation
pub struct FetchResult {
    /// The raw string content of the file
    pub content: String,
    /// The detected or specified format of the file
    pub format: FileFormat,
}

/// Fetches a protein structure from RCSB by its PDB ID
///
/// This function first attempts to download the mmCIF version, and falls back
/// to the PDB version if mmCIF is not available
pub fn fetch_pdb(pdb_identifier_code: &str) -> Result<FetchResult, String> {
    let pdb_identifier_code = pdb_identifier_code.to_uppercase();

    // Try mmCIF first (preferred format), fall back to PDB
    let url_fetch_attempts_collection = [
        (
            format!("{}/{}.cif", RCSB_BASE_URL, pdb_identifier_code),
            FileFormat::Cif,
        ),
        (
            format!("{}/{}.pdb", RCSB_BASE_URL, pdb_identifier_code),
            FileFormat::Pdb,
        ),
    ];

    for (target_download_url, expected_file_format) in &url_fetch_attempts_collection {
        match reqwest::blocking::get(target_download_url) {
            Ok(http_response_data) if http_response_data.status().is_success() => {
                let fetched_file_content =
                    http_response_data
                        .text()
                        .map_err(|network_response_error| {
                            format!("Failed to read response: {}", network_response_error)
                        })?;
                return Ok(FetchResult {
                    content: fetched_file_content,
                    format: *expected_file_format,
                });
            }
            _ => continue,
        }
    }

    Err(format!(
        "Failed to fetch PDB code '{}' from RCSB",
        pdb_identifier_code
    ))
}

/// Loads a protein structure from the local file system
///
/// The format is determined based on the file extension (.cif and .mmcif are
/// treated as CIF, everything else as PDB)
pub fn load_file(file_system_path: &str) -> Result<FetchResult, String> {
    let loaded_file_content =
        std::fs::read_to_string(file_system_path).map_err(|file_system_read_error| {
            format!(
                "Failed to read file '{}': {}",
                file_system_path, file_system_read_error
            )
        })?;

    // Determine format from extension
    let detected_file_format =
        if file_system_path.ends_with(".cif") || file_system_path.ends_with(".mmcif") {
            FileFormat::Cif
        } else {
            FileFormat::Pdb
        };

    Ok(FetchResult {
        content: loaded_file_content,
        format: detected_file_format,
    })
}
