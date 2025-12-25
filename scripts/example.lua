-- Example script
-- This script runs on startup and when modified (hot-reload)

-- Fetch a protein from the PDB
local fetched_protein_data = pdb.fetch("8W1K")

-- Set representation mode
-- Options: "spheres", "backbone", "both"
fetched_protein_data:representation("both")

-- Set color scheme
-- Options: "chain", "element", "bfactor", "secondary"
fetched_protein_data:color_by("bfactor")

-- You can also set a uniform color
-- fetched_protein_data:color(1.0, 0.5, 0.2)  -- Orange
