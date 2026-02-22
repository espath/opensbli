#!/usr/bin/env python3
"""Generate ParaView-ready XDMF + cropped HDF5 for OpenSBLI 2D multi-block output."""

from __future__ import annotations

import argparse
from pathlib import Path
import h5py
import numpy as np


def _crop_center(arr, ny_target, nx_target):
    ny, nx = arr.shape
    if ny_target > ny or nx_target > nx:
        raise RuntimeError(f"Target shape {(ny_target, nx_target)} larger than source {(ny, nx)}")
    oy = (ny - ny_target) // 2
    ox = (nx - nx_target) // 2
    return arr[oy : oy + ny_target, ox : ox + nx_target]


def build_files(sol_h5: Path, grid_h5: Path, viz_h5: Path, xdmf_path: Path) -> None:
    with h5py.File(sol_h5, "r") as sf, h5py.File(grid_h5, "r") as gf, h5py.File(viz_h5, "w") as vf:
        blocks = sorted([k for k in sf.keys() if k.startswith("opensbliblock")])
        if not blocks:
            raise RuntimeError(f"No opensbliblock* groups in {sol_h5}")

        lines = [
            '<?xml version="1.0" ?>',
            '<Xdmf Version="3.0">',
            "  <Domain>",
            '    <Grid Name="OpenSBLI_2D_MultiBlock" GridType="Collection" CollectionType="Spatial">',
        ]

        for block_name in blocks:
            block_id = int(block_name[-2:])
            x_name = f"x0_B{block_id}"
            y_name = f"x1_B{block_id}"

            # Use physical sizes from constants in solution file.
            nx_target = int(np.asarray(sf[f"block{block_id}np0"][()]).reshape(-1)[0])
            ny_target = int(np.asarray(sf[f"block{block_id}np1"][()]).reshape(-1)[0])

            gx = gf[block_name][x_name][:]
            gy = gf[block_name][y_name][:]
            x = _crop_center(gx, ny_target, nx_target)
            y = _crop_center(gy, ny_target, nx_target)

            out_block = vf.create_group(block_name)
            out_block.create_dataset(x_name, data=x)
            out_block.create_dataset(y_name, data=y)

            field_names = sorted(
                [k for k in sf[block_name].keys() if k.endswith(f"_B{block_id}") and k not in (x_name, y_name)]
            )
            common_names = {}
            for fld in field_names:
                out_block.create_dataset(fld, data=_crop_center(sf[block_name][fld][:], ny_target, nx_target))
                # Map block-specific names (e.g. rho_B0) to common names (rho) for ParaView Merge Blocks.
                suffix = f"_B{block_id}"
                common_names[fld] = fld[:-len(suffix)] if fld.endswith(suffix) else fld

            lines += [
                f'      <Grid Name="{block_name}" GridType="Uniform">',
                f'        <Topology TopologyType="2DSMesh" Dimensions="{ny_target} {nx_target}"/>',
                '        <Geometry GeometryType="X_Y">',
                f'          <DataItem Dimensions="{ny_target} {nx_target}" NumberType="Float" Precision="8" Format="HDF">{viz_h5.name}:/{block_name}/{x_name}</DataItem>',
                f'          <DataItem Dimensions="{ny_target} {nx_target}" NumberType="Float" Precision="8" Format="HDF">{viz_h5.name}:/{block_name}/{y_name}</DataItem>',
                "        </Geometry>",
            ]

            for fld in field_names:
                attr_name = common_names[fld]
                lines += [
                    f'        <Attribute Name="{attr_name}" AttributeType="Scalar" Center="Node">',
                    f'          <DataItem Dimensions="{ny_target} {nx_target}" NumberType="Float" Precision="8" Format="HDF">{viz_h5.name}:/{block_name}/{fld}</DataItem>',
                    "        </Attribute>",
                ]
            lines += ["      </Grid>"]

        lines += ["    </Grid>", "  </Domain>", "</Xdmf>"]

    xdmf_path.write_text("\n".join(lines) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="../opensbli_output.h5")
    p.add_argument("--grid", default="../data.h5")
    p.add_argument("--output", default="read_airfoil_2d.xdmf")
    p.add_argument("--viz-h5", default="airfoil_2d_viz.h5")
    args = p.parse_args()

    sol_h5 = Path(args.input).resolve()
    grid_h5 = Path(args.grid).resolve()
    xdmf_path = Path(args.output).resolve()
    viz_h5 = Path(args.viz_h5).resolve()

    if not sol_h5.exists():
        raise FileNotFoundError(sol_h5)
    if not grid_h5.exists():
        raise FileNotFoundError(grid_h5)

    build_files(sol_h5, grid_h5, viz_h5, xdmf_path)
    print(f"Wrote {viz_h5}")
    print(f"Wrote {xdmf_path}")


if __name__ == "__main__":
    main()
