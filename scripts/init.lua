-- Protein viewer - example script.
-- This script runs on startup and when modified (hot-reload).

print("Protein analysis script.")

-- Fetch a protein from the PDB.
local fetched_protein_data = pdb.fetch("1AKE")

-- Print basic information.
print(fetched_protein_data:info())

-- Set representation mode.
-- Options: "spheres", "backbone", "both".
fetched_protein_data:representation("both")

-- Set color scheme.
-- Options: "chain", "element", "bfactor", "secondary".
fetched_protein_data:color_by("bfactor")

-- You can also set a uniform color:
-- fetched_protein_data:color(1.0, 0.5, 0.2)  -- Orange.

print("Controls.")
print("  Click: Select/deselect atoms (highlighted yellow).")
print("  Mouse drag: Rotate.")
print("  Scroll: Zoom.")
print("  D: Measure distances between selected atoms.")
print("  X: Clear selection.")
print("  1: Spheres only.")
print("  2: Backbone only.")
print("  3: Both.")
print("  C: Color by chain.")
print("  B: Color by B-factor.")
print("  E: Color by element.")
print("  R: Reset camera.")
print("  Esc: Quit.")
print("")
print("Edit scripts/init.lua to reload!")