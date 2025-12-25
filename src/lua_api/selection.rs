//! Lua interface for atom selections

use mlua::UserData;
use std::sync::{Arc, RwLock};
use glam::Vec3;
use crate::selection::SelectionSet;
use crate::protein::structure::ProteinData;

/// A Lua-exposed wrapper around a SelectionSet with a reference to its protein
#[derive(Clone)]
pub struct LuaSelection {
    pub inner_selection_set: SelectionSet,
    pub protein_reference: Arc<RwLock<ProteinData>>,
}

impl LuaSelection {
    pub fn new(selection_set: SelectionSet, protein_reference: Arc<RwLock<ProteinData>>) -> Self {
        Self { 
            inner_selection_set: selection_set,
            protein_reference 
        }
    }

    /// Calculates the geometric centroid of the selected atoms
    pub fn calculate_centroid(&self) -> Vec3 {
        let locked_protein_data = self.protein_reference.read().unwrap();
        let mut position_sum_vector = Vec3::ZERO;
        let mut selected_atoms_count = 0;

        if let Some(pdb) = &locked_protein_data.underlying_pdb_data {
            for &atom_index in &self.inner_selection_set.atom_indices {
                if let Some(atom_reference) = pdb.atoms().nth(atom_index) {
                    let position_tuple = atom_reference.pos();
                    position_sum_vector += Vec3::new(
                        position_tuple.0 as f32,
                        position_tuple.1 as f32,
                        position_tuple.2 as f32,
                    );
                    selected_atoms_count += 1;
                }
            }
        }

        if selected_atoms_count > 0 {
            position_sum_vector / selected_atoms_count as f32
        } else {
            Vec3::ZERO
        }
    }
}

impl UserData for LuaSelection {
    fn add_methods<'lua, M: mlua::UserDataMethods<'lua, Self>>(methods: &mut M) {
        // selection:get_selected_atom_count() returns the number of selected atoms
        methods.add_method("get_selected_atom_count", |_, this, ()| {
            Ok(this.inner_selection_set.count())
        });

        // selection:calculate_geometric_centroid() returns the (x, y, z) coordinates of the selection's geometric center
        methods.add_method("calculate_geometric_centroid", |_, this, ()| {
            let centroid_vector = this.calculate_centroid();
            Ok((centroid_vector.x, centroid_vector.y, centroid_vector.z))
        });

        // selection:calculate_distance_to_other_selection(other) returns the distance between centroids of two selections
        methods.add_method("calculate_distance_to_other_selection", |_, this, other_selection: LuaSelection| {
            let centroid_a = this.calculate_centroid();
            let centroid_b = other_selection.calculate_centroid();
            Ok(centroid_a.distance(centroid_b))
        });

        // selection:calculate_angle_between_three_selections(sel2, sel3) returns the angle between three selection centroids in degrees
        methods.add_method("calculate_angle_between_three_selections", |_, this, (selection_two, selection_three): (LuaSelection, LuaSelection)| {
            let vertex_a = this.calculate_centroid();
            let vertex_b = selection_two.calculate_centroid();
            let vertex_c = selection_three.calculate_centroid();
            
            let vector_ba = (vertex_a - vertex_b).normalize();
            let vector_bc = (vertex_c - vertex_b).normalize();
            Ok(vector_ba.dot(vector_bc).acos().to_degrees())
        });

        // selection:calculate_dihedral_between_four_selections(sel2, sel3, sel4) returns the dihedral angle between four selection centroids in degrees
        methods.add_method("calculate_dihedral_between_four_selections", |_, this, (selection_two, selection_three, selection_four): (LuaSelection, LuaSelection, LuaSelection)| {
            let vertex_a = this.calculate_centroid();
            let vertex_b = selection_two.calculate_centroid();
            let vertex_c = selection_three.calculate_centroid();
            let vertex_d = selection_four.calculate_centroid();
            
            Ok(crate::analysis::dihedrals::calculate_dihedral_angle_between_points(vertex_a, vertex_b, vertex_c, vertex_d))
        });

        // selection:translate_selected_atoms(x, y, z) moves all selected atoms by the given offset
        methods.add_method_mut("translate_selected_atoms", |_, this, (offset_x, offset_y, offset_z): (f32, f32, f32)| {
            let mut mutable_protein_data = this.protein_reference.write().unwrap();
            let translation_vector = Vec3::new(offset_x, offset_y, offset_z);
            
            if let Some(pdb) = &mut mutable_protein_data.underlying_pdb_data {
                for &atom_index in &this.inner_selection_set.atom_indices {
                    if let Some(atom_reference) = pdb.atoms_mut().nth(atom_index) {
                        let current_pos = atom_reference.pos();
                        let new_pos = Vec3::new(current_pos.0 as f32, current_pos.1 as f32, current_pos.2 as f32) + translation_vector;
                        atom_reference.set_pos((new_pos.x as f64, new_pos.y as f64, new_pos.z as f64)).unwrap();
                    }
                }
            }
            mutable_protein_data.structural_data_revision_number += 1;
            Ok(())
        });

        // selection:rotate_selected_atoms(angle_degrees, axis_x, axis_y, axis_z, origin_x, origin_y, origin_z)
        // Rotates the selection around an axis passing through an origin point.
        methods.add_method_mut("rotate_selected_atoms", |_, this, (angle_degrees, axis_x, axis_y, axis_z, origin_x, origin_y, origin_z): (f32, f32, f32, f32, f32, f32, f32)| {
            let mut mutable_protein_data = this.protein_reference.write().unwrap();
            let rotation_axis = Vec3::new(axis_x, axis_y, axis_z).normalize();
            let rotation_origin = Vec3::new(origin_x, origin_y, origin_z);
            let rotation_quaternion = glam::Quat::from_axis_angle(rotation_axis, angle_degrees.to_radians());
            
            if let Some(pdb) = &mut mutable_protein_data.underlying_pdb_data {
                for &atom_index in &this.inner_selection_set.atom_indices {
                    if let Some(atom_reference) = pdb.atoms_mut().nth(atom_index) {
                        let current_pos_tuple = atom_reference.pos();
                        let current_pos_vector = Vec3::new(current_pos_tuple.0 as f32, current_pos_tuple.1 as f32, current_pos_tuple.2 as f32);
                        
                        let relative_pos = current_pos_vector - rotation_origin;
                        let rotated_relative_pos = rotation_quaternion * relative_pos;
                        let new_world_pos = rotated_relative_pos + rotation_origin;
                        
                        atom_reference.set_pos((new_world_pos.x as f64, new_world_pos.y as f64, new_world_pos.z as f64)).unwrap();
                    }
                }
            }
            mutable_protein_data.structural_data_revision_number += 1;
            Ok(())
        });

        // selection:create_union_with_other_selection(other) returns a new selection combining both
        methods.add_method("create_union_with_other_selection", |_, this, other_selection: LuaSelection| {
            Ok(LuaSelection::new(
                this.inner_selection_set.union(&other_selection.inner_selection_set),
                this.protein_reference.clone()
            ))
        });

        // selection:create_intersection_with_other_selection(other) returns a new selection of atoms in both
        methods.add_method("create_intersection_with_other_selection", |_, this, other_selection: LuaSelection| {
            Ok(LuaSelection::new(
                this.inner_selection_set.intersection(&other_selection.inner_selection_set),
                this.protein_reference.clone()
            ))
        });
    }
}

impl<'lua> mlua::FromLua<'lua> for LuaSelection {
    fn from_lua(value: mlua::Value<'lua>, _lua: &'lua mlua::Lua) -> mlua::Result<Self> {
        match value {
            mlua::Value::UserData(user_data_handle) => Ok(user_data_handle.borrow::<Self>()?.clone()),
            _ => Err(mlua::Error::FromLuaConversionError {
                from: value.type_name(),
                to: "LuaSelection",
                message: Some("Expected a LuaSelection object".to_string()),
            }),
        }
    }
}