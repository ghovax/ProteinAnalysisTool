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
/// This function first checks the local `.cache` folder. If not found, it attempts 
/// to download the mmCIF version, falling back to the PDB version, and saves the 
/// result to the cache.
pub fn fetch_pdb(pdb_identifier_code: &str) -> Result<FetchResult, String> {
    let pdb_identifier_code_uppercase = pdb_identifier_code.to_uppercase();
    let cache_directory_path = std::path::Path::new(".cache");

    // Ensure cache directory exists
    if !cache_directory_path.exists() {
        std::fs::create_dir_all(cache_directory_path).map_err(|io_error_encountered| {
            format!("Failed to create cache directory: {}", io_error_encountered)
        })?;
    }

    // 1. Check cache for mmCIF first, then PDB
    let cached_mmcif_file_path = cache_directory_path.join(format!("{}.cif", pdb_identifier_code_uppercase));
    if cached_mmcif_file_path.exists() {
        return load_file(cached_mmcif_file_path.to_str().unwrap());
    }

    let cached_pdb_file_path = cache_directory_path.join(format!("{}.pdb", pdb_identifier_code_uppercase));
    if cached_pdb_file_path.exists() {
        return load_file(cached_pdb_file_path.to_str().unwrap());
    }

    // 2. Not in cache, try downloading mmCIF first, then PDB
    let download_attempt_configuration_collection = [
        (
            format!("{}/{}.cif", RCSB_BASE_URL, pdb_identifier_code_uppercase),
            FileFormat::Cif,
            ".cif"
        ),
        (
            format!("{}/{}.pdb", RCSB_BASE_URL, pdb_identifier_code_uppercase),
            FileFormat::Pdb,
            ".pdb"
        ),
    ];

    for (target_remote_download_url, target_file_format, target_file_extension) in &download_attempt_configuration_collection {
        match reqwest::blocking::get(target_remote_download_url) {
            Ok(http_response_data) if http_response_data.status().is_success() => {
                let downloaded_file_content_string =
                    http_response_data
                        .text()
                        .map_err(|network_read_error| {
                            format!("Failed to read response: {}", network_read_error)
                        })?;
                
                // Save to cache
                let output_cache_file_full_path = cache_directory_path.join(format!("{}{}", pdb_identifier_code_uppercase, target_file_extension));
                std::fs::write(&output_cache_file_full_path, &downloaded_file_content_string).map_err(|cache_write_io_error| {
                    format!("Failed to write to cache: {}", cache_write_io_error)
                })?;

                return Ok(FetchResult {
                    content: downloaded_file_content_string,
                    format: *target_file_format,
                });
            }
            _ => continue,
        }
    }

    Err(format!(
        "Failed to fetch PDB code '{}' from RCSB",
        pdb_identifier_code_uppercase
    ))
}

/// Loads a protein structure from the local file system
///
/// The format is determined based on the file extension (.cif and .mmcif are
/// treated as CIF, everything else as PDB)
pub fn load_file(file_system_path_string: &str) -> Result<FetchResult, String> {
    let loaded_file_content_string =
        std::fs::read_to_string(file_system_path_string).map_err(|file_system_io_error| {
            format!(
                "Failed to read file '{}': {}",
                file_system_path_string, file_system_io_error
            )
        })?;

    // Determine format from extension
    let detected_protein_file_format =
        if file_system_path_string.ends_with(".cif") || file_system_path_string.ends_with(".mmcif") {
            FileFormat::Cif
        } else {
            FileFormat::Pdb
        };

    Ok(FetchResult {
        content: loaded_file_content_string,
        format: detected_protein_file_format,
    })
}
