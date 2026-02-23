#!/usr/bin/env python3
"""Generate ParaView-ready mesh-only XDMF for OpenSBLI structured 2D grids."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py


def _crop_by_halo(arr, d_m, d_p):
    # OpenSBLI metadata: d_m negative left halo, d_p positive right halo.
    left_i = int(abs(d_m[0])) if d_m is not None else 0
    left_j = int(abs(d_m[1])) if d_m is not None and len(d_m) > 1 else 0
    right_i = int(d_p[0]) if d_p is not None else 0
    right_j = int(d_p[1]) if d_p is not None and len(d_p) > 1 else 0
    return arr[left_j : arr.shape[0] - right_j, left_i : arr.shape[1] - right_i]


def build_files(mesh_h5: Path, viz_h5: Path, xdmf_path: Path, block_name: str = "opensbliblock00", crop_halo: bool = True) -> None:
    with h5py.File(mesh_h5, "r") as f, h5py.File(viz_h5, "w") as vf:
        if block_name not in f:
            raise RuntimeError(f"Block group '{block_name}' not found in {mesh_h5}")
        g = f[block_name]

        x_names = sorted([k for k in g.keys() if k.startswith("x0_")])
        y_names = sorted([k for k in g.keys() if k.startswith("x1_")])
        if not x_names or not y_names:
            raise RuntimeError(f"Could not find x0_*/x1_* datasets in {mesh_h5}:{block_name}")

        x_name = x_names[0]
        y_name = y_names[0]
        x = g[x_name][:]
        y = g[y_name][:]

        if crop_halo:
            d_m = g[x_name].attrs.get("d_m")
            d_p = g[x_name].attrs.get("d_p")
            x = _crop_by_halo(x, d_m, d_p)
            y = _crop_by_halo(y, d_m, d_p)

        out_block = vf.create_group(block_name)
        out_block.create_dataset(x_name, data=x)
        out_block.create_dataset(y_name, data=y)
        ny, nx = x.shape

    lines = [
        '<?xml version="1.0" ?>',
        '<Xdmf Version="3.0">',
        "  <Domain>",
        f'    <Grid Name="{block_name}" GridType="Uniform">',
        f'      <Topology TopologyType="2DSMesh" Dimensions="{ny} {nx}"/>',
        '      <Geometry GeometryType="X_Y">',
        f'        <DataItem Dimensions="{ny} {nx}" NumberType="Float" Precision="8" Format="HDF">{viz_h5.name}:/{block_name}/{x_name}</DataItem>',
        f'        <DataItem Dimensions="{ny} {nx}" NumberType="Float" Precision="8" Format="HDF">{viz_h5.name}:/{block_name}/{y_name}</DataItem>',
        "      </Geometry>",
        "    </Grid>",
        "  </Domain>",
        "</Xdmf>",
        "",
    ]
    xdmf_path.write_text("\n".join(lines))


def main() -> None:
    p = argparse.ArgumentParser(description="Generate XDMF wrapper for OpenSBLI mesh HDF5")
    p.add_argument("--input", required=True, help="Input mesh HDF5 path")
    p.add_argument("--output", help="Output XDMF path (default: input stem + .xdmf)")
    p.add_argument("--viz-h5", help="Output visualization HDF5 (default: input stem + _viz.h5)")
    p.add_argument("--block", default="opensbliblock00", help="Block group name")
    p.add_argument("--keep-halo", action="store_true", help="Do not crop halo layers")
    args = p.parse_args()

    mesh_h5 = Path(args.input).resolve()
    if not mesh_h5.exists():
        raise FileNotFoundError(mesh_h5)

    xdmf = Path(args.output).resolve() if args.output else mesh_h5.with_suffix(".xdmf")
    viz_h5 = Path(args.viz_h5).resolve() if args.viz_h5 else mesh_h5.with_name(mesh_h5.stem + "_viz.h5")

    build_files(mesh_h5, viz_h5, xdmf, block_name=args.block, crop_halo=not args.keep_halo)
    print(f"Wrote {viz_h5}")
    print(f"Wrote {xdmf}")


if __name__ == "__main__":
    main()
