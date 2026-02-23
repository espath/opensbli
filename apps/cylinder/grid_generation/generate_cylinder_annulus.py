#!/usr/bin/env python3
"""Generate structured annulus grids for cylinder cases and write OpenSBLI data.h5."""

import argparse
from pathlib import Path

import h5py
import numpy as np
from opensbli import SimulationBlock, output_hdf5


def stretch_inner_only(n, beta):
    """Monotone map in [0,1] with clustering only near 0 for beta>1."""
    s = np.linspace(0.0, 1.0, n)
    if beta <= 1.0:
        return s
    # Power-law map: eta=s**beta
    # beta>1 => small spacing near inner wall (s=0), monotone growth outward.
    return s**beta


def fill_halos_linear(arr, nhalo):
    """Linear extrapolation halos in both directions."""
    out = np.zeros((arr.shape[0] + 2 * nhalo, arr.shape[1] + 2 * nhalo))
    out[nhalo:-nhalo, nhalo:-nhalo] = arr

    # i-direction halos
    for h in range(1, nhalo + 1):
        out[nhalo - h, nhalo:-nhalo] = out[nhalo, nhalo:-nhalo] - h * (
            out[nhalo + 1, nhalo:-nhalo] - out[nhalo, nhalo:-nhalo]
        )
        out[-nhalo - 1 + h, nhalo:-nhalo] = out[-nhalo - 1, nhalo:-nhalo] + h * (
            out[-nhalo - 1, nhalo:-nhalo] - out[-nhalo - 2, nhalo:-nhalo]
        )

    # j-direction halos
    for h in range(1, nhalo + 1):
        out[:, nhalo - h] = out[:, nhalo] - h * (out[:, nhalo + 1] - out[:, nhalo])
        out[:, -nhalo - 1 + h] = out[:, -nhalo - 1] + h * (
            out[:, -nhalo - 1] - out[:, -nhalo - 2]
        )

    return out


def build_annulus(nr, ntheta, r_inner, r_outer, theta0, theta1, beta_r):
    eta = stretch_inner_only(nr, beta_r)  # clustered near wall (r_inner) only
    r = r_inner + (r_outer - r_inner) * eta
    th = np.linspace(theta0, theta1, ntheta)

    rr, tt = np.meshgrid(r, th, indexing="ij")
    x = rr * np.cos(tt)
    y = rr * np.sin(tt)
    return x, y


def write_opensbli_data(output_path, x, y, nhalo):
    xh = fill_halos_linear(x, nhalo).T
    yh = fill_halos_linear(y, nhalo).T

    tmp = output_path.with_suffix(".tmp.h5")
    with h5py.File(tmp, "w") as f:
        f.create_dataset("x0", data=xh)
        f.create_dataset("x1", data=yh)

    b = SimulationBlock(2, block_number=0)
    with h5py.File(tmp, "r") as f:
        arrays = [f["x0"], f["x1"]]
        array_names = ["x0", "x1"]
        npoints = [x.shape[0], x.shape[1]]
        halos = [(-nhalo, nhalo), (-nhalo, nhalo)]
        output_hdf5(arrays, array_names, halos, npoints, b, **{"filename": str(output_path)})

    tmp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description="Generate structured cylinder-annulus data.h5")
    parser.add_argument("--nr", type=int, default=361, help="Radial points (interior)")
    parser.add_argument("--ntheta", type=int, default=361, help="Angular points (interior)")
    parser.add_argument("--r-inner", type=float, default=0.5, help="Cylinder radius")
    parser.add_argument("--r-outer", type=float, default=20.0, help="Far-field radius")
    parser.add_argument("--theta0", type=float, default=np.pi / 2.0, help="Start angle [rad]")
    parser.add_argument("--theta1", type=float, default=np.pi, help="End angle [rad]")
    parser.add_argument(
        "--beta-r",
        type=float,
        default=2.0,
        help="Radial clustering exponent (>1 clusters near inner wall only)",
    )
    parser.add_argument("--nhalo", type=int, default=5, help="Halo width")
    parser.add_argument("--output", default="data.h5", help="Output HDF5 path")
    args = parser.parse_args()

    if args.r_outer <= args.r_inner:
        raise ValueError("--r-outer must be > --r-inner")
    if args.theta1 <= args.theta0:
        raise ValueError("--theta1 must be > --theta0")

    x, y = build_annulus(
        args.nr, args.ntheta, args.r_inner, args.r_outer, args.theta0, args.theta1, args.beta_r
    )
    out = Path(args.output).resolve()
    write_opensbli_data(out, x, y, args.nhalo)
    print(f"Wrote {out}")
    print(f"Interior size: ({args.nr}, {args.ntheta}), halo={args.nhalo}")
    print(f"Theta range: [{args.theta0}, {args.theta1}] rad")


if __name__ == "__main__":
    main()
