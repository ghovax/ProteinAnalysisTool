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

### 7. Structural Analysis
*   **Dihedral Angles**: Precise calculation of Phi, Psi, and Omega angles for every residue.
*   **Ramachandran Data**: Export backbone torsion data for validation and plotting.

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
*   `pdb.fetch(code)`: Fetches a protein from RCSB by its PDB identifier.
*   `pdb.load(path)`: Loads a protein from a local file (PDB or mmCIF).
*   `pdb.list()`: Returns a table of names for all currently loaded protein identifiers.

### Global `camera` Module
*   `camera.get_pos()`: Returns the current camera world coordinates (x, y, z).
*   `camera.set_target(x, y, z)`: Sets the focus point of the camera.
*   `camera.set_params(dist, yaw, pitch)`: Sets the spherical coordinates of the camera.

### Global `session` Module
*   `session.save(path)`: Saves the current application state (proteins and camera) to a Lua script.
*   `session.load(path)`: Restores a session by executing a previously saved script.

### Protein Object Methods (`protein:method`)
*   `protein:select(query)`: Selects atoms using the PyMOL-style query language.
*   `protein:representation(mode)`: Sets visual mode (`spheres`, `backbone`, `both`, `sticks`, `ball-and-stick`, `space-filling`, `lines`).
*   `protein:calculate_molecular_surface(resolution, threshold)`: Generates a molecular surface mesh.
*   `protein:set_surface_visibility(visible)`: Toggles surface rendering.
*   `protein:rmsd(other_protein)`: Calculates RMSD relative to another protein.
*   `protein:superpose(reference_protein)`: Aligns this protein to a reference using Kabsch algorithm.
*   `protein:ramachandran_data()`: Returns Phi/Psi angles for all residues.
*   `protein:hydrogen_bonds()`: Detects and returns hydrogen bonds.
*   `protein:color_by(scheme)`: Sets coloring scheme (`chain`, `element`, `bfactor`, `secondary`).
*   `protein:color(r, g, b)`: Sets a uniform RGB color.
*   `protein:atoms()` / `protein:residues()`: Returns detailed atomic and residue information.

## Terminal Logging
The application uses a standardized logging system. All analysis results (like distance measurements), file loading status, and Lua `print()` output are directed to the terminal console. Run with `RUST_LOG=info` to see all messages.

## Hot-Reloading

The application watches the `scripts/` directory and the current directory for changes to `.lua` files. Saving a script will automatically re-execute it in the running application.