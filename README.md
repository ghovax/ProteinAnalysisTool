# Protein Viewer

A high-performance protein visualization tool built with Rust, WGPU, and Lua.

## Controls

### Mouse

* **Left Click & Drag**: Rotate camera
* **Scroll**: Zoom in/out

### Keyboard Shortcuts

* `1`: Set representation to Spheres only
* `2`: Set representation to Backbone only
* `3`: Set representation to Both (Backbone and Spheres)
* `C`: Color by Chain
* `B`: Color by B-factor
* `E`: Color by Element
* `S`: Color by Secondary structure
* `R`: Reset Camera focus and selection
* `M`: Calculate distance between selected atoms
* `Esc`: Quit application

## Scripting API (Lua)

The application embeds a Lua engine for automation and analysis. Scripts can be placed in `scripts/init.lua` and support hot-reloading.

### Global `pdb` Module

* `pdb.fetch(code)`: Fetches a protein from RCSB by its PDB identifier
* `pdb.load(path)`: Loads a protein from a local file (PDB or mmCIF)
* `pdb.list()`: Returns a table of names for all currently loaded protein identifiers

### Protein Object Methods (`protein:method`)
* `protein:name()`: Returns the name/ID of the protein
* `protein:info()`: Returns a summary string with protein information
* `protein:atom_count()`: Returns the total number of atoms
* `protein:residue_count()`: Returns an approximate count of residues
* `protein:chains()`: Returns a table of chain identifiers
* `protein:show()` / `protein:hide()`: Controls visibility
* `protein:representation(mode)`: Sets visual mode ("spheres", "backbone", "both")
* `protein:color_by(scheme)`: Sets coloring scheme ("chain", "element", "bfactor", "secondary")
* `protein:color(r, g, b)`: Sets a uniform RGB color
* `protein:atoms()`: Returns a detailed table of all atoms
* `protein:residues(chain_id)`: Returns a table of residues, optionally filtered by chain

## Hot-Reloading

The application watches the `scripts/` directory and the current directory for changes to `.lua` files. Saving a script will automatically re-execute it in the running application.
