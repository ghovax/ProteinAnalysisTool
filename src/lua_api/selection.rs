//! Lua interface for atom selections

use mlua::UserData;
use crate::selection::SelectionSet;

/// A Lua-exposed wrapper around a SelectionSet
#[derive(Clone)]
pub struct LuaSelection {
    pub inner_selection_set: SelectionSet,
}

impl LuaSelection {
    pub fn new(selection_set: SelectionSet) -> Self {
        Self { inner_selection_set: selection_set }
    }
}

impl UserData for LuaSelection {
    fn add_methods<'lua, M: mlua::UserDataMethods<'lua, Self>>(methods: &mut M) {
        // selection:count() returns the number of selected atoms
        methods.add_method("count", |_, this, ()| {
            Ok(this.inner_selection_set.count())
        });

        // selection:union(other) returns a new selection combining both
        methods.add_method("union", |_, this, other_selection: LuaSelection| {
            Ok(LuaSelection::new(this.inner_selection_set.union(&other_selection.inner_selection_set)))
        });

        // selection:intersection(other) returns a new selection of atoms in both
        methods.add_method("intersection", |_, this, other_selection: LuaSelection| {
            Ok(LuaSelection::new(this.inner_selection_set.intersection(&other_selection.inner_selection_set)))
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
