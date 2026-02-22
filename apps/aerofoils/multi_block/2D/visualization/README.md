# Airfoil 2D Visualization

Generate an XDMF wrapper for ParaView from OpenSBLI multi-block 2D output.

## Usage

From `apps/aerofoils/multi_block/2D/visualization`:

```bash
python3 generate_xdmf.py --input ../opensbli_output.h5 --grid ../data.h5 --output read_airfoil_2d.xdmf
```

This writes two files:
- `airfoil_2d_viz.h5` (cropped to physical points; halo layers removed)
- `read_airfoil_2d.xdmf` (references `airfoil_2d_viz.h5`)

Then open `read_airfoil_2d.xdmf` in ParaView.

## Notes

- The script expects block groups `opensbliblock00`, `opensbliblock01`, etc.
- Coordinates are read from `data.h5`.
- Solution fields are read from `opensbli_output.h5`.
- The script crops arrays to `block*np*` sizes from the solution file to remove halo/ghost layers.
