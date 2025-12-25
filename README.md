# Protein Viewer

A high-performance protein visualization tool built with Rust, WGPU, and Lua.

## Features

### 1. Advanced Selection Language
A PyMOL-style selection engine for precise atom filtering.
*   **Keywords**: `chain`, `resn`, `resi`, `name`, `elem`, `backbone`, `sidechain`, `helix`, `sheet`, `all`, `none`.
*   **Operators**: `and`, `or`, `not`, and parentheses for grouping.
*   **Spatial Queries**: `within 5.0 of <selection>` for proximity-based selection.

### 2. Full Secondary Structure Support
Real-time parsing of `HELIX` and `SHEET` records from PDB files for accurate structural annotation and coloring.

### 3. Comprehensive Rendering Modes
High-quality visualization of molecular structures:
*   **Spheres**: Alpha carbon (CA) or all-atom spheres.
*   **Sticks**: Cylindrical bonds between atoms.
*   **Ball & Stick**: Combination of atom spheres and bond cylinders.
*   **Space-Filling**: Van der Waals (VdW) radius-based spheres.
*   **Backbone Trace**: Smooth line segments through the protein backbone.

### 4. Parallelized Surface Rendering
Fast molecular surface generation using a parallelized Marching Cubes algorithm and Signed Distance Fields (SDF).
*   **Solvent-Accessible Surface (SAS)** based on VdW radii.
*   Adjustable voxel resolution and isosurface threshold.

### 5. RMSD & Protein Superposition
Compare and align structures using the **Kabsch Algorithm**.
*   Calculate Root-Mean-Square Deviation (RMSD) between protein pairs.
*   Optimal rigid-body superposition (rotation + translation).

### 6. Geometric Hydrogen Bond Detection
Automated identification of hydrogen bonds based on donor-acceptor distances and geometric criteria.

### 7. Physical & Chemical Interactions
Detection of complex biochemical interactions beyond simple bonding.
*   **Salt Bridges**: Identification of charged residue interactions (Arg/Lys and Asp/Glu).

### 8. Structural Analysis
*   **Dihedral Angles**: Precise calculation of Phi, Psi, and Omega angles for every residue.
*   **Ramachandran Data**: Export backbone torsion data for validation and plotting.
*   **SASA**: Solvent Accessible Surface Area calculation using the Shrake-Rupley algorithm.
*   **RMSF**: Root-Mean-Square Fluctuation analysis across structural ensembles (NMR models or MD snapshots).

### 9. Data Manipulation & Export
Tools for programmatically modifying and exporting structural data.
*   **Atomic Manipulation**: Direct translation and rotation of selections.
*   **Structured Export**: Built-in support for JSON, CSV, and FASTA formats.

## Controls

### Mouse

*   **Left Click & Drag**: Rotate camera (Orbit)
*   **Right Click & Drag**: Translate camera (Pan)
*   **Left Click (on Atom)**: Toggle atom selection (Max 2 atoms)
*   **Scroll**: Zoom in/out

### Keyboard Shortcuts

#### Representation Modes
*   `1`: Spheres only (Alpha Carbons)
*   `2`: Backbone trace only
*   `3`: Both Spheres and Backbone
*   `4`: Sticks (Cylindrical bonds)
*   `5`: Ball and Stick (Default)
*   `6`: Space-Filling (Van der Waals spheres)
*   `7`: Lines (Wireframe bonds)

#### Coloring Schemes
*   `C`: Color by Chain identifier
*   `B`: Color by B-factor (Temperature factor)
*   `E`: Color by chemical Element
*   `S`: Color by Secondary structure (Helix/Sheet/Loop)

#### Analysis and Surface
*   `M`: Calculate distance between selected atoms (Logged to terminal console)
*   `F`: Toggle Molecular Surface visibility (if computed)
*   `R`: Reset Camera focus, clear selections and measurements
*   `Esc`: Quit application

## Scripting API (Lua)

The application embeds a Lua engine for automation and analysis. Scripts support hot-reloading.

### Global `pdb` Module
*   `pdb.fetch_protein_from_rcsb(code)`: Fetches a protein from RCSB by its PDB identifier.
*   `pdb.load_protein_from_local_file(path)`: Loads a protein from a local file (PDB or mmCIF).
*   `pdb.get_loaded_protein_identifiers()`: Returns a list of all currently loaded protein identifiers.

### Global `camera` Module
*   `camera.get_camera_world_position()`: Returns the current camera world coordinates (x, y, z).
*   `camera.set_camera_focus_target(x, y, z)`: Sets the focus point of the camera.
*   `camera.set_camera_spherical_parameters(dist, yaw, pitch)`: Sets the spherical coordinates of the camera.

### Global `session` Module
*   `session.save_application_session(path)`: Saves current state to a Lua script.
*   `session.load_application_session(path)`: Restores a session by executing a previously saved script.

### Protein Object Methods (`protein:method`)
*   `protein:get_protein_name()`: Returns the name/ID of the protein.
*   `protein:get_total_atom_count()`: Returns the total number of atoms.
*   `protein:select_atoms_by_query(query)`: Returns a Selection object matching the PyMOL-style query.
*   `protein:set_representation_mode(mode)`: Sets visual mode (`spheres`, `backbone`, `sticks`, etc.).
*   `protein:set_color_scheme_by_property(scheme)`: Sets coloring scheme (`chain`, `element`, `bfactor`, `secondary`).
*   `protein:set_uniform_rgb_color(r, g, b)`: Sets a uniform RGB color.
*   `protein:calculate_molecular_surface(res, threshold)`: Generates a molecular surface mesh.
*   `protein:calculate_solvent_accessible_surface_area()`: Returns total SASA in square Angstroms.
*   `protein:calculate_root_mean_square_deviation(other, selection)`: Calculates RMSD.
*   `protein:calculate_root_mean_square_fluctuation(selection)`: Calculates RMSF across all models.
*   `protein:calculate_ramachandran_dihedral_angles()`: Returns Phi/Psi angles for all residues.
*   `protein:superimpose_onto_reference_structure(ref, selection)`: Aligns structure using Kabsch algorithm.
*   `protein:detect_geometric_hydrogen_bonds()`: Detects and returns hydrogen bonds.
*   `protein:detect_salt_bridge_interactions()`: Identifies charged residue interactions.
*   `protein:export_analysis_report_to_json(path)`: Exports detailed data to JSON.
*   `protein:export_residue_data_to_csv(path)`: Exports residue angles to CSV.
*   `protein:get_sequence_for_chain(id)` / `protein:generate_fasta_formatted_sequence()`: Sequence analysis and export.
*   `protein:set_visibility_on()` / `protein:set_visibility_off()`: Controls rendering.
*   `protein:get_all_atom_data()` / `protein:get_all_residue_data()`: Returns detailed structural tables.

### Selection Object Methods (`selection:method`)
*   `selection:get_selected_atom_count()`: Returns number of atoms.
*   `selection:calculate_geometric_centroid()`: Returns (x, y, z) coordinates.
*   `selection:calculate_distance_to_other_selection(other)`: Distance between selection centroids.
*   `selection:calculate_angle_between_three_selections(s2, s3)`: Angle between three centroids.
*   `selection:calculate_dihedral_between_four_selections(s2, s3, s4)`: Dihedral between four centroids.
*   `selection:translate_selected_atoms(x, y, z)`: Moves atoms by given offset.
*   `selection:rotate_selected_atoms(angle, axis_x...origin_x...)`: Rotates atoms around an axis.
*   `selection:create_union_with_other_selection(other)`: Returns new combined selection.

## Terminal Logging
The application uses a standardized logging system. All analysis results, file loading status, and Lua `print()` output are directed to the terminal console. Run with `RUST_LOG=info` to see all messages.

## Hot-Reloading
The application watches for changes to `.lua` files in the `scripts/` and root directories, re-executing them automatically upon save.
