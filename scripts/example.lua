-- This script demonstrates and tests the newly implemented biotechnology toolkit features.

-- 1. Fetching proteins for superposition comparison
-- Using two similar structures of Ubiquitin (1UBQ and 1YYF)
-- 1YYF is an NMR structure with multiple models, 1UBQ is a single model.
local reference_protein = pdb.fetch_protein_from_rcsb("1YYF")
local moving_protein = pdb.fetch_protein_from_rcsb("1YYF")
pdb.fetch_protein_from_rcsb("1UBQ")
pdb.fetch_protein_from_rcsb("9P9K")

-- 2. Setting up initial visual representations
reference_protein:set_representation_mode("ball_and_stick")
reference_protein:set_color_scheme_by_property("bfactor_value") -- Color by B-factor

moving_protein:set_representation_mode("sticks")
moving_protein:set_color_scheme_by_property("secondary_structure") -- Color by Helix/Sheet/Loop

-- 3. Testing RMSD and Kabsch Superposition with Selections
-- We select only Model 1 and Alpha Carbons to ensure identical length (76 atoms)
local common_selection_query = "model 1 and name CA"
local reference_selection = reference_protein:select_atoms_by_query(common_selection_query)
local moving_selection = moving_protein:select_atoms_by_query(common_selection_query)

print("Reference atoms selected: " .. reference_selection:get_selected_atom_count())
print("Moving atoms selected: " .. moving_selection:get_selected_atom_count())

-- Calculate RMSD using the selection
local initial_rmsd_value = moving_protein:calculate_root_mean_square_deviation(reference_protein, moving_selection)
print("Initial RMSD (Selection) before superposition: " .. string.format("%.3f", initial_rmsd_value) .. " Angstroms")

print("Performing Kabsch superposition using selection...")
moving_protein:superimpose_onto_reference_structure(reference_protein, moving_selection)

local final_rmsd_value = moving_protein:calculate_root_mean_square_deviation(reference_protein, moving_selection)
print("Final RMSD (Selection) after superposition: " .. string.format("%.3f", final_rmsd_value) .. " Angstroms")

-- 4. Testing the Advanced Selection Language
print("Evaluating selection queries...")
local helix_selection = reference_protein:select_atoms_by_query("helix")
print("Number of atoms identified in helices: " .. helix_selection:get_selected_atom_count())

local sheet_selection = reference_protein:select_atoms_by_query("sheet")
print("Number of atoms identified in beta sheets: " .. sheet_selection:get_selected_atom_count())

-- Complex spatial query: atoms within 5A of any Valine residue
local proximity_selection = reference_protein:select_atoms_by_query("within 5.0 of resn VAL")
print("Atoms within 5.0A of VAL residues: " .. proximity_selection:get_selected_atom_count())

-- 5. Testing Parallelized Surface Rendering (Marching Cubes)
-- print("Computing molecular surface mesh (Voxel size: 0.5A)...")
-- reference_protein:calculate_molecular_surface(0.5, 0.0)
-- reference_protein:set_surface_visibility(true)

-- 6. Testing Hydrogen Bond Detection
local detected_hydrogen_bonds = reference_protein:detect_geometric_hydrogen_bonds()
print("Detected " .. #detected_hydrogen_bonds .. " geometric hydrogen bonds in reference protein")

-- 7. Testing Structural Analysis (Ramachandran Data)
local ramachandran_analysis_data = reference_protein:calculate_ramachandran_dihedral_angles()
if #ramachandran_analysis_data > 0 then
    local first_residue_entry = ramachandran_analysis_data[1]
    print("Analysis for residue: index " .. first_residue_entry.residue_number)
    print("Resulting residue's angles: Phi: " .. string.format("%.2f", first_residue_entry.phi) .. "°, Psi: " .. string.format("%.2f", first_residue_entry.psi) .. "°")
end

-- 8. Testing New Advanced Toolkit Features

print("\n--- Advanced Toolkit Features ---")

-- SASA Calculation
local sasa_value = reference_protein:calculate_solvent_accessible_surface_area()
print("Total SASA of reference protein: " .. string.format("%.2f", sasa_value) .. " A^2")

-- Sequence Analysis
local seq_a = reference_protein:get_sequence_for_chain("A")
print("Sequence of Chain A (first 10 residues): " .. seq_a:sub(1, 10) .. "...")

-- Selection Centroid Geometry
local sel1 = reference_protein:select_atoms_by_query("resi 1-10")
local sel2 = reference_protein:select_atoms_by_query("resi 20-30")
local dist = sel1:calculate_distance_to_other_selection(sel2)
print("Distance between centroids of resi 1-10 and 20-30: " .. string.format("%.2f", dist) .. " A")

-- Salt Bridge Detection
local salt_bridges = reference_protein:detect_salt_bridge_interactions()
print("Detected " .. #salt_bridges .. " potential salt bridges")

-- Atomic Manipulation (Testing translation)
print("Translating selection resi 1-5 by 10A in X...")
local manip_sel = reference_protein:select_atoms_by_query("resi 1-5")
manip_sel:translate_selected_atoms(10.0, 0.0, 0.0)

-- Data Export
print("Exporting analysis to JSON and CSV...")
reference_protein:export_analysis_report_to_json("analysis_report.json")
reference_protein:export_residue_data_to_csv("residues_data.csv")

-- FASTA Export
local fasta_data = reference_protein:generate_fasta_formatted_sequence()
print("Generated FASTA header: " .. fasta_data:sub(1, 20) .. "...")

-- 9. Camera and Session Control
camera.set_camera_focus_target(0.0, 0.0, 0.0)
camera.set_camera_spherical_parameters(75.0, 0.8, 0.4)
