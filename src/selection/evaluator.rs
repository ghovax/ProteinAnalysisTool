//! Evaluation of selection expressions

use std::collections::HashSet;
use glam::Vec3;
use pdbtbx::{ContainsAtomConformerResidue, ContainsAtomConformerResidueChain, ContainsAtomConformerResidueChainModel};

use super::{SelectionExpression, SelectionSet};
use crate::protein::structure::ProteinData;

pub struct Evaluator<'a> {
    protein: &'a ProteinData,
}

impl<'a> Evaluator<'a> {
    pub fn new(protein: &'a ProteinData) -> Self {
        Self { protein }
    }

    pub fn evaluate_expression(&self, selection_expression: &SelectionExpression) -> SelectionSet {
        let selected_atom_indices = match selection_expression {
            SelectionExpression::All => (0..self.protein.get_total_atom_count()).collect(),
            SelectionExpression::None => HashSet::new(),
            SelectionExpression::Chain(chain_identifier) => self.select_atoms_by_chain_identifier(chain_identifier),
            SelectionExpression::ResidueName(residue_name) => self.select_atoms_by_residue_name(residue_name),
            SelectionExpression::ResidueNumber(residue_number) => self.select_atoms_by_residue_number(*residue_number),
            SelectionExpression::ResidueRange(start_number, end_number) => self.select_atoms_by_residue_number_range(*start_number, *end_number),
            SelectionExpression::AtomName(atom_name) => self.select_atoms_by_atom_name(atom_name),
            SelectionExpression::Element(element_symbol) => self.select_atoms_by_element_symbol(element_symbol),
            SelectionExpression::Model(model_index) => self.select_atoms_by_model_index(*model_index),
            SelectionExpression::Backbone => self.select_backbone_atoms(),
            SelectionExpression::Sidechain => self.select_sidechain_atoms(),
            SelectionExpression::Helix => self.select_atoms_in_helices(),
            SelectionExpression::Sheet => self.select_atoms_in_sheets(),
            SelectionExpression::And(left_expression, right_expression) => {
                let left_set = self.evaluate_expression(left_expression);
                let right_set = self.evaluate_expression(right_expression);
                left_set.intersection(&right_set).atom_indices
            }
            SelectionExpression::Or(left_expression, right_expression) => {
                let left_set = self.evaluate_expression(left_expression);
                let right_set = self.evaluate_expression(right_expression);
                left_set.union(&right_set).atom_indices
            }
            SelectionExpression::Not(inner_expression) => {
                let all_atom_indices: HashSet<usize> = (0..self.protein.get_total_atom_count()).collect();
                let inner_set = self.evaluate_expression(inner_expression);
                all_atom_indices.iter().filter(|index| !inner_set.contains(**index)).cloned().collect()
            }
            SelectionExpression::Within(distance_threshold, inner_expression) => self.select_atoms_within_distance(*distance_threshold, inner_expression),
        };
        SelectionSet::new(selected_atom_indices)
    }

    fn select_atoms_by_chain_identifier(&self, target_chain_identifier: &str) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .enumerate()
                .filter(|(_, atom_hierarchy)| atom_hierarchy.chain().id() == target_chain_identifier)
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_by_residue_name(&self, target_residue_name: &str) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .enumerate()
                .filter(|(_, atom_hierarchy)| atom_hierarchy.residue().name().map(|n| n.trim()) == Some(target_residue_name))
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_by_residue_number(&self, target_residue_number: isize) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .enumerate()
                .filter(|(_, atom_hierarchy)| atom_hierarchy.residue().serial_number() == target_residue_number)
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_by_residue_number_range(&self, start_residue_number: isize, end_residue_number: isize) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .enumerate()
                .filter(|(_, atom_hierarchy)| {
                    let current_residue_number = atom_hierarchy.residue().serial_number();
                    current_residue_number >= start_residue_number && current_residue_number <= end_residue_number
                })
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_by_atom_name(&self, target_atom_name: &str) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms()
                .enumerate()
                .filter(|(_, atom_reference)| atom_reference.name().trim() == target_atom_name)
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_by_element_symbol(&self, target_element_symbol: &str) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms()
                .enumerate()
                .filter(|(_, atom_reference)| atom_reference.element().map_or(false, |element| element.symbol() == target_element_symbol))
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_by_model_index(&self, target_model_index: usize) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .enumerate()
                .filter(|(_, atom_hierarchy)| {
                    // Models in pdbtbx might be 0-based or use their serial number. 
                    // We'll compare with model serial number if available, or just use iteration order if needed.
                    // Standard PDB MODEL records use 1-based indices.
                    atom_hierarchy.model().serial_number() == target_model_index
                })
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_backbone_atoms(&self) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            let backbone_atom_names = ["N", "CA", "C", "O", "OXT"];
            pdb.atoms()
                .enumerate()
                .filter(|(_, atom_reference)| backbone_atom_names.contains(&atom_reference.name().trim()))
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_sidechain_atoms(&self) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            let backbone_atom_names = ["N", "CA", "C", "O", "OXT"];
            pdb.atoms()
                .enumerate()
                .filter(|(_, atom_reference)| !backbone_atom_names.contains(&atom_reference.name().trim()))
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_in_helices(&self) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .enumerate()
                .filter(|(_, atom_hierarchy)| pdb.is_residue_in_helix(atom_hierarchy.chain().id(), atom_hierarchy.residue().serial_number()))
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_in_sheets(&self) -> HashSet<usize> {
        if let Some(pdb) = &self.protein.underlying_pdb_data {
            pdb.atoms_with_hierarchy()
                .enumerate()
                .filter(|(_, atom_hierarchy)| pdb.is_residue_in_sheet(atom_hierarchy.chain().id(), atom_hierarchy.residue().serial_number()))
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }

    fn select_atoms_within_distance(&self, distance_threshold: f32, source_selection_expression: &SelectionExpression) -> HashSet<usize> {
        let source_selection_set = self.evaluate_expression(source_selection_expression);
        if source_selection_set.count() == 0 {
            return HashSet::new();
        }

        if let Some(pdb) = &self.protein.underlying_pdb_data {
            // Get coordinates of all atoms in the source selection set
            let source_atom_coordinates: Vec<Vec3> = source_selection_set.atom_indices.iter()
                .map(|atom_index| {
                    let atom_reference = pdb.atoms().nth(*atom_index).unwrap();
                    let position_tuple = atom_reference.pos();
                    Vec3::new(position_tuple.0 as f32, position_tuple.1 as f32, position_tuple.2 as f32)
                })
                .collect();

            let distance_threshold_squared = distance_threshold * distance_threshold;
            pdb.atoms()
                .enumerate()
                .filter(|(_, target_atom_reference)| {
                    let position_tuple = target_atom_reference.pos();
                    let target_atom_position = Vec3::new(position_tuple.0 as f32, position_tuple.1 as f32, position_tuple.2 as f32);
                    source_atom_coordinates.iter().any(|&source_coordinate| source_coordinate.distance_squared(target_atom_position) <= distance_threshold_squared)
                })
                .map(|(index, _)| index)
                .collect()
        } else {
            HashSet::new()
        }
    }
}
