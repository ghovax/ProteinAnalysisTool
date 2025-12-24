-- Protein Viewer - Example Script
-- This script runs on startup and when modified (hot-reload)

print("=== Protein Analysis Script ===")

-- Fetch a protein from the PDB
-- Try some interesting examples:
--   1AKE - Adenylate Kinase (small, good for testing)
--   1HHO - Hemoglobin
--   4HHB - Hemoglobin (larger)
--   6LU7 - SARS-CoV-2 Main Protease

local protein = pdb.fetch("1AKE")

-- Print basic information
print(protein:info())

-- Get specific data
local chains = protein:chains()
print("Chain IDs:", table.concat(chains, ", "))

-- Get center of mass
local cx, cy, cz = protein:center_of_mass()
print("Center of mass:", cx, cy, cz)

-- Get bounding box
local minx, miny, minz, maxx, maxy, maxz = protein:bounding_box()
print("Bounding box:")
print("  Min:", minx, miny, minz)
print("  Max:", maxx, maxy, maxz)

-- Count residues per chain
for _, chain_id in ipairs(chains) do
    local residues = protein:residues(chain_id)
    local count = 0
    for _ in pairs(residues) do count = count + 1 end
    print("Chain " .. chain_id .. ": " .. count .. " residues")
end

print("=== Ready ===")
