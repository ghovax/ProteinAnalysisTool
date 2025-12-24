# Protein Viewer

A high-performance protein visualization tool built with Rust, WGPU, and Lua.

## Controls

### Mouse
* **Left Click & Drag**: Rotate camera
* **Scroll**: Zoom in/out

### Keyboard Shortcuts
* **1**: Set representation to Spheres only
* **2**: Set representation to Backbone only
* **3**: Set representation to Both (Backbone and Spheres)
* **C**: Color by Chain
* **B**: Color by B-factor
* **E**: Color by Element
* **R**: Reset Camera focus
* **Esc**: Quit application

## Scripting API (Lua)

The application embeds a Lua engine for automation and analysis. Scripts can be placed in `scripts/init.lua` and support hot-reloading.

### Global `pdb` Module
* `pdb.fetch(code)`: Fetches a protein from RCSB by its PDB identifier
* `pdb.load(path)`: Loads a protein from a local file (PDB or mmCIF)
* `pdb.list()`: Returns a table of names for all currently loaded protein identifiers

### Protein Object Methods (`p:method`)
* `p:name()`: Returns the name/ID of the protein
* `p:info()`: Returns a summary string with protein information
* `p:atom_count()`: Returns the total number of atoms
* `p:residue_count()`: Returns an approximate count of residues
* `p:chains()`: Returns a table of chain identifiers
* `p:show()` / `p:hide()`: Controls visibility
* `p:representation(mode)`: Sets visual mode ("spheres", "backbone", "both")
* `p:color_by(scheme)`: Sets coloring scheme ("chain", "element", "bfactor", "secondary")
* `p:color(r, g, b)`: Sets a uniform RGB color
* `p:atoms()`: Returns a detailed table of all atoms
* `p:residues(chain_id)`: Returns a table of residues, optionally filtered by chain

## Hot-Reloading
The application watches the `scripts/` directory and the current directory for changes to `.lua` files. Saving a script will automatically re-execute it in the running application.
