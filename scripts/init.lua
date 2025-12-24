-- Protein Viewer - Example Script
-- This script runs on startup and when modified (hot-reload)

print("=== Protein Analysis Script ===")

-- Fetch a protein from the PDB
local protein = pdb.fetch("1AKE")

-- Print basic information
print(protein:info())

-- Set representation mode
-- Options: "spheres", "backbone", "both"
protein:representation("both")

-- Set color scheme
-- Options: "chain", "element", "bfactor", "secondary"
protein:color_by("bfactor")

-- You can also set a uniform color:
-- protein:color(1.0, 0.5, 0.2)  -- orange

print("=== Controls ===")
print("  Mouse drag: Rotate")
print("  Scroll: Zoom")
print("  1: Spheres only")
print("  2: Backbone only")
print("  3: Both")
print("  C: Color by chain")
print("  B: Color by B-factor")
print("  R: Reset camera")
print("  Esc: Quit")
print("")
print("Edit scripts/init.lua to reload!")
