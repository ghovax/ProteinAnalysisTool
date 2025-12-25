# Protein Viewer

A high-performance protein visualization tool built with Rust, WGPU, and Lua.

## Features

### 1. Advanced Selection Language

A PyMOL-style selection engine for precise atom filtering.

- **Keywords**: `chain`, `resn`, `resi`, `name`, `elem`, `backbone`, `sidechain`, `helix`, `sheet`, `all`, `none`.
- **Operators**: `and`, `or`, `not`, and parentheses for grouping.
- **Spatial Queries**: `within 5.0 of <selection>` for proximity-based selection.

### 2. Full Secondary Structure Support

Real-time parsing of `HELIX` and `SHEET` records from PDB files for accurate structural annotation and coloring.

### 3. Comprehensive Rendering Modes

High-quality visualization of molecular structures:

- **Spheres**: Alpha carbon (CA) or all-atom spheres.
- **Sticks**: Cylindrical bonds between atoms.
- **Ball & Stick**: Combination of atom spheres and bond cylinders.
- **Space-Filling**: Van der Waals (VdW) radius-based spheres.
- **Backbone Trace**: Smooth line segments through the protein backbone.

### 4. Parallelized Surface Rendering

Fast molecular surface generation using a parallelized Marching Cubes algorithm and Signed Distance Fields (SDF).

- **Solvent-Accessible Surface (SAS)** based on VdW radii.
- Adjustable voxel resolution and isosurface threshold.

### 5. RMSD & Protein Superposition

Compare and align structures using the **Kabsch Algorithm**.

- Calculate Root-Mean-Square Deviation (RMSD) between protein pairs.
- Optimal rigid-body superposition (rotation + translation).

### 6. Geometric Hydrogen Bond Detection

Automated identification of hydrogen bonds based on donor-acceptor distances and geometric criteria.

### 7. Structural Analysis

- **Dihedral Angles**: Precise calculation of Phi, Psi, and Omega angles for every residue.
- **Ramachandran Data**: Export backbone torsion data for validation and plotting.

## Controls

### Mouse

- **Left Click & Drag**: Rotate camera
- **Scroll**: Zoom in/out

### Keyboard Shortcuts

- `1`: Set representation to Spheres only
- `2`: Set representation to Backbone only
- `3`: Set representation to Both (Backbone and Spheres)
- `C`: Color by Chain
- `B`: Color by B-factor
- `E`: Color by Element
- `S`: Color by Secondary structure
- `R`: Reset Camera focus and selection
- `M`: Calculate distance between selected atoms
- `Esc`: Quit application

## Scripting API (Lua)

The application embeds a Lua engine for automation and analysis. Scripts support hot-reloading.

### Global `pdb` Module

- `pdb.fetch(code)`: Fetches a protein from RCSB by its PDB identifier
- `pdb.load(path)`: Loads a protein from a local file (PDB or mmCIF)
- `pdb.list()`: Returns a table of names for all currently loaded protein identifiers

### Protein Object Methods (`protein:method`)

- `protein:select(query)`: Selects atoms using the PyMOL-style query language.
- `protein:representation(mode)`: Sets visual mode (`spheres`, `backbone`, `both`, `sticks`, `ball-and-stick`, `space-filling`, `lines`).
- `protein:calculate_molecular_surface(resolution, threshold)`: Generates a molecular surface mesh.
- `protein:set_surface_visibility(visible)`: Toggles surface rendering.
- `protein:rmsd(other_protein)`: Calculates RMSD relative to another protein.
- `protein:superpose(reference_protein)`: Aligns this protein to a reference.
- `protein:ramachandran_data()`: Returns Phi/Psi angles for all residues.
- `protein:hydrogen_bonds()`: Detects and returns hydrogen bonds.
- `protein:color_by(scheme)`: Sets coloring scheme (`chain`, `element`, `bfactor`, `secondary`).
- `protein:color(r, g, b)`: Sets a uniform RGB color.
- `protein:atoms()` / `protein:residues()`: Returns detailed atomic and residue information.

## Hot-Reloading

The application watches the `scripts/` directory and the current directory for changes to `.lua` files. Saving a script will automatically re-execute it in the running application.
