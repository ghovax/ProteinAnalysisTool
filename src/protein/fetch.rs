const RCSB_BASE_URL: &str = "https://files.rcsb.org/download";

#[derive(Clone, Copy, Debug)]
pub enum FileFormat {
    Pdb,
    Cif,
}

pub struct FetchResult {
    pub content: String,
    pub format: FileFormat,
}

pub fn fetch_pdb(code: &str) -> Result<FetchResult, String> {
    let code = code.to_uppercase();

    // Try mmCIF first (preferred format), fall back to PDB
    let attempts = [
        (format!("{}/{}.cif", RCSB_BASE_URL, code), FileFormat::Cif),
        (format!("{}/{}.pdb", RCSB_BASE_URL, code), FileFormat::Pdb),
    ];

    for (url, format) in &attempts {
        match reqwest::blocking::get(url) {
            Ok(response) if response.status().is_success() => {
                let content = response.text().map_err(|e| format!("Failed to read response: {}", e))?;
                return Ok(FetchResult { content, format: *format });
            }
            _ => continue,
        }
    }

    Err(format!("Failed to fetch PDB code '{}' from RCSB", code))
}

pub fn load_file(path: &str) -> Result<FetchResult, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read file '{}': {}", path, e))?;

    // Determine format from extension
    let format = if path.ends_with(".cif") || path.ends_with(".mmcif") {
        FileFormat::Cif
    } else {
        FileFormat::Pdb
    };

    Ok(FetchResult { content, format })
}
