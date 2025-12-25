//! Protein sequence analysis and FASTA export

use crate::protein::structure::ProteinData;

/// Returns the primary sequence of a protein chain in one-letter code
pub fn extract_sequence_from_chain(
    protein_data: &ProteinData,
    chain_identifier: &str,
) -> Result<String, String> {
    let chain_reference = protein_data.pdb.chains()
        .find(|current_chain| current_chain.id() == chain_identifier)
        .ok_or_else(|| format!("Chain {} not found", chain_identifier))?;

    let mut sequence_string = String::new();
    for current_residue in chain_reference.residues() {
        let residue_name = current_residue.name().unwrap_or("UNK");
        let one_letter_code = match residue_name {
            "ALA" => 'A', "CYS" => 'C', "ASP" => 'D', "GLU" => 'E',
            "PHE" => 'F', "GLY" => 'G', "HIS" => 'H', "ILE" => 'I',
            "LYS" => 'K', "LEU" => 'L', "MET" => 'M', "ASN" => 'N',
            "PRO" => 'P', "GLN" => 'Q', "ARG" => 'R', "SER" => 'S',
            "THR" => 'T', "VAL" => 'V', "TRP" => 'W', "TYR" => 'Y',
            _ => 'X',
        };
        sequence_string.push(one_letter_code);
    }
    Ok(sequence_string)
}

/// Generates a FASTA-formatted string for the protein structure
pub fn generate_fasta_formatted_string(
    protein_data: &ProteinData,
) -> String {
    let mut fasta_output_string = String::new();
    for current_chain in protein_data.pdb.chains() {
        if let Ok(sequence) = extract_sequence_from_chain(protein_data, current_chain.id()) {
            fasta_output_string.push_str(&format!(">{}:{}
", protein_data.name, current_chain.id()));
            fasta_output_string.push_str(&sequence);
            fasta_output_string.push('\n');
        }
    }
    fasta_output_string
}
