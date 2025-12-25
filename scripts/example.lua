-- This script demonstrates and tests the newly implemented biotechnology toolkit features.

-- 1. Fetching proteins for superposition comparison
-- Using two similar structures of Ubiquitin (1UBQ and 1YYF)
-- 1YYF is an NMR structure with multiple models, 1UBQ is a single model.
local reference_protein = pdb.fetch("1YYF")
local moving_protein = pdb.fetch("1YYF")

-- 2. Setting up initial visual representations
reference_protein:representation("ball-and-stick")
reference_protein:color_by("bfactor") -- Color by Helix/Sheet/Loop

moving_protein:representation("sticks")
moving_protein:color_by("secondary") -- Color by Helix/Sheet/Loop

-- 3. Testing RMSD and Kabsch Superposition with Selections
-- We select only Model 1 and Alpha Carbons to ensure identical length (76 atoms)
local common_selection_query = "model 1 and name CA"
local reference_selection = reference_protein:select(common_selection_query)
local moving_selection = moving_protein:select(common_selection_query)

print("Reference atoms selected: " .. reference_selection:count())
print("Moving atoms selected: " .. moving_selection:count())

-- Calculate RMSD using the selection
local initial_rmsd_value = moving_protein:rmsd(reference_protein, moving_selection)
print("Initial RMSD (Selection) before superposition: " .. string.format("%.3f", initial_rmsd_value) .. " Angstroms")

print("Performing Kabsch superposition using selection...")
moving_protein:superpose(reference_protein, moving_selection)

local final_rmsd_value = moving_protein:rmsd(reference_protein, moving_selection)
print("Final RMSD (Selection) after superposition: " .. string.format("%.3f", final_rmsd_value) .. " Angstroms")

-- 4. Testing the Advanced Selection Language
print("Evaluating selection queries...")
local helix_selection = reference_protein:select("helix")
print("Number of atoms identified in helices: " .. helix_selection:count())

local sheet_selection = reference_protein:select("sheet")
print("Number of atoms identified in beta sheets: " .. sheet_selection:count())

-- Complex spatial query: atoms within 5A of any Valine residue
local proximity_selection = reference_protein:select("within 5.0 of resn VAL")
print("Atoms within 5.0A of VAL residues: " .. proximity_selection:count())

-- 5. Testing Parallelized Surface Rendering (Marching Cubes)
-- print("Computing molecular surface mesh (Voxel size: 0.5A)...")
-- reference_protein:calculate_molecular_surface(0.5, 0.0)
-- reference_protein:set_surface_visibility(true)

-- 6. Testing Hydrogen Bond Detection
local detected_hydrogen_bonds = reference_protein:hydrogen_bonds()
print("Detected " .. #detected_hydrogen_bonds .. " geometric hydrogen bonds in reference protein")

-- 7. Testing Structural Analysis (Ramachandran Data)
local ramachandran_analysis_data = reference_protein:ramachandran_data()
if #ramachandran_analysis_data > 0 then
    local first_residue_entry = ramachandran_analysis_data[1]
    print("Analysis for residue: index " .. first_residue_entry.residue_number)
    print("       Phi: " .. string.format("%.2f", first_residue_entry.phi) .. "°, Psi: " .. string.format("%.2f", first_residue_entry.psi) .. "°")
end

-- 8. Camera and Session Control
camera.set_target(0.0, 0.0, 0.0)
camera.set_params(75.0, 0.8, 0.4)
