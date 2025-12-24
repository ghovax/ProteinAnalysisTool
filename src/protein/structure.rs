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
    pub fn from_string(content: &str, name: &str, format: FileFormat) -> Result<Self, String> {
        let pdbtbx_format = match format {
            FileFormat::Pdb => Format::Pdb,
            FileFormat::Cif => Format::Mmcif,
        };

        let (pdb, errors) = ReadOptions::default()
            .set_level(pdbtbx::StrictnessLevel::Loose)
            .set_format(pdbtbx_format)
            .read_raw(std::io::BufReader::new(std::io::Cursor::new(content.as_bytes())))
            .map_err(|e| format!("Parse error: {:?}", e))?;

        if !errors.is_empty() {
            eprintln!("Warnings while parsing {}: {:?}", name, errors);
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
            .map(|c| c.id().to_string())
            .collect()
    }

    pub fn center_of_mass(&self) -> Vec3 {
        let mut sum = Vec3::ZERO;
        let mut count = 0;

        for atom in self.pdb.atoms() {
            let pos = atom.pos();
            sum += Vec3::new(pos.0 as f32, pos.1 as f32, pos.2 as f32);
            count += 1;
        }

        if count > 0 {
            sum / count as f32
        } else {
            Vec3::ZERO
        }
    }

    pub fn bounding_box(&self) -> (Vec3, Vec3) {
        let mut min = Vec3::splat(f32::MAX);
        let mut max = Vec3::splat(f32::MIN);

        for atom in self.pdb.atoms() {
            let pos = atom.pos();
            let p = Vec3::new(pos.0 as f32, pos.1 as f32, pos.2 as f32);
            min = min.min(p);
            max = max.max(p);
        }

        (min, max)
    }

    pub fn ca_positions(&self) -> Vec<(Vec3, String)> {
        let mut positions = Vec::new();

        for chain in self.pdb.chains() {
            let chain_id = chain.id().to_string();
            for residue in chain.residues() {
                for conformer in residue.conformers() {
                    for atom in conformer.atoms() {
                        if atom.name() == "CA" {
                            let pos = atom.pos();
                            positions.push((
                                Vec3::new(pos.0 as f32, pos.1 as f32, pos.2 as f32),
                                chain_id.clone(),
                            ));
                        }
                    }
                }
            }
        }

        positions
    }

    /// Returns backbone segments as pairs of (start, end, chain_id) for line rendering
    pub fn backbone_segments(&self) -> Vec<(Vec3, Vec3, String)> {
        let mut segments = Vec::new();

        for chain in self.pdb.chains() {
            let chain_id = chain.id().to_string();
            let mut prev_ca: Option<Vec3> = None;

            for residue in chain.residues() {
                for conformer in residue.conformers() {
                    for atom in conformer.atoms() {
                        if atom.name() == "CA" {
                            let pos = atom.pos();
                            let current = Vec3::new(pos.0 as f32, pos.1 as f32, pos.2 as f32);

                            if let Some(prev) = prev_ca {
                                // Only connect if distance is reasonable (< 5 Angstroms)
                                if prev.distance(current) < 5.0 {
                                    segments.push((prev, current, chain_id.clone()));
                                }
                            }
                            prev_ca = Some(current);
                        }
                    }
                }
            }
        }

        segments
    }

    /// Get B-factor range for coloring
    pub fn bfactor_range(&self) -> (f32, f32) {
        let mut min = f32::MAX;
        let mut max = f32::MIN;

        for atom in self.pdb.atoms() {
            let bf = atom.b_factor() as f32;
            min = min.min(bf);
            max = max.max(bf);
        }

        (min, max)
    }

    /// Get CA atoms with their B-factors
    pub fn ca_with_bfactor(&self) -> Vec<(Vec3, String, f32)> {
        let mut positions = Vec::new();

        for chain in self.pdb.chains() {
            let chain_id = chain.id().to_string();
            for residue in chain.residues() {
                for conformer in residue.conformers() {
                    for atom in conformer.atoms() {
                        if atom.name() == "CA" {
                            let pos = atom.pos();
                            positions.push((
                                Vec3::new(pos.0 as f32, pos.1 as f32, pos.2 as f32),
                                chain_id.clone(),
                                atom.b_factor() as f32,
                            ));
                        }
                    }
                }
            }
        }

        positions
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

    pub fn fetch(&mut self, code: &str) -> Result<Arc<RwLock<ProteinData>>, String> {
        let code = code.to_uppercase();

        if let Some(existing) = self.proteins.get(&code) {
            return Ok(existing.clone());
        }

        let result = fetch_pdb(&code)?;
        let protein = ProteinData::from_string(&result.content, &code, result.format)?;
        let arc = Arc::new(RwLock::new(protein));
        self.proteins.insert(code, arc.clone());
        Ok(arc)
    }

    pub fn load(&mut self, path: &str) -> Result<Arc<RwLock<ProteinData>>, String> {
        let name = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        if let Some(existing) = self.proteins.get(&name) {
            return Ok(existing.clone());
        }

        let result = load_file(path)?;
        let protein = ProteinData::from_string(&result.content, &name, result.format)?;
        let arc = Arc::new(RwLock::new(protein));
        self.proteins.insert(name, arc.clone());
        Ok(arc)
    }

    pub fn list(&self) -> Vec<String> {
        self.proteins.keys().cloned().collect()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Arc<RwLock<ProteinData>>> {
        self.proteins.values()
    }
}
