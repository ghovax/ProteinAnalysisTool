-- Advanced Protein Viewer API Test Script
-- This script demonstrates and tests the newly implemented biotechnology toolkit features.

-- 1. Fetching proteins for superposition comparison
-- Using two similar structures of Ubiquitin (1UBQ and 1YYF)
local reference_protein = pdb.fetch("1UBQ")
local moving_protein = pdb.fetch("1YYF")

-- 2. Setting up initial visual representations
reference_protein:representation("ball-and-stick")
reference_protein:color_by("bfactor") -- Color by Helix/Sheet/Loop

moving_protein:representation("sticks")
-- moving_protein:color(0.4, 0.7, 1.0) -- Set a custom light blue color

-- 3. Testing RMSD and Kabsch Superposition
local initial_rmsd_value = moving_protein:rmsd(reference_protein)
print("Initial RMSD before superposition: " .. string.format("%.3f", initial_rmsd_value) .. " Angstroms")

print("Performing Kabsch superposition...")
moving_protein:superpose(reference_protein)

local final_rmsd_value = moving_protein:rmsd(reference_protein)
print("Final RMSD after superposition: " .. string.format("%.3f", final_rmsd_value) .. " Angstroms")

-- 4. Testing the Advanced Selection Language
print("Evaluating selection queries...")
local helix_atom_indices = reference_protein:select("helix")
print("Number of atoms identified in helices: " .. #helix_atom_indices)

local sheet_atom_indices = reference_protein:select("sheet")
print("Number of atoms identified in beta sheets: " .. #sheet_atom_indices)

-- Complex spatial query: atoms within 5A of any Valine residue
local proximity_selection = reference_protein:select("within 5.0 of resn VAL")
print("Atoms within 5.0A of VAL residues: " .. #proximity_selection)

-- 5. Testing Parallelized Surface Rendering (Marching Cubes)
print("Computing molecular surface mesh (Voxel size: 0.5A)...")
reference_protein:calculate_molecular_surface(0.5, 0.0)
reference_protein:set_surface_visibility(true)

-- 6. Testing Hydrogen Bond Detection
local detected_hydrogen_bonds = reference_protein:hydrogen_bonds()
print("Detected " .. #detected_hydrogen_bonds .. " geometric hydrogen bonds in reference protein")

-- 7. Testing Structural Analysis (Ramachandran Data)
local ramachandran_analysis_data = reference_protein:ramachandran_data()
if #ramachandran_analysis_data > 0 then
    local first_residue_entry = ramachandran_analysis_data[1]
    print("Analysis for residue: " .. first_residue_entry.residue_name .. " " .. first_residue_entry.residue_number)
    print("       Phi: " .. string.format("%.2f", first_residue_entry.phi) .. "°, Psi: " .. string.format("%.2f", first_residue_entry.psi) .. "°")
end

-- 8. Camera and Session Control
camera.set_target(0.0, 0.0, 0.0)
camera.set_params(75.0, 0.8, 0.4)

-- Optional: Save the current session state
-- session.save("scripts/test_session.lua")