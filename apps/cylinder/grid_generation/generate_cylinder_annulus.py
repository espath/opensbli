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


def extend_1d_linear(values, nhalo):
    """Linear halo extension in index space for 1D coordinates."""
    out = np.zeros(values.size + 2 * nhalo)
    out[nhalo:-nhalo] = values
    d0 = values[1] - values[0]
    d1 = values[-1] - values[-2]
    for h in range(1, nhalo + 1):
        out[nhalo - h] = values[0] - h * d0
        out[-nhalo - 1 + h] = values[-1] + h * d1
    return out


def build_annulus(nr, ntheta, r_inner, r_outer, theta0, theta1, beta_r):
    eta = stretch_inner_only(nr, beta_r)  # clustered near wall (r_inner) only
    r = r_inner + (r_outer - r_inner) * eta
    span = theta1 - theta0
    endpoint = not np.isclose(abs(span), 2.0 * np.pi, rtol=0.0, atol=1.0e-12)
    th = np.linspace(theta0, theta1, ntheta, endpoint=endpoint)

    rr, tt = np.meshgrid(r, th, indexing="ij")
    x = rr * np.cos(tt)
    y = rr * np.sin(tt)
    return x, y


def build_annulus_with_halos(nr, ntheta, r_inner, r_outer, theta0, theta1, beta_r, nhalo):
    """Build annulus coordinates including halos in parametric (r,theta) space."""
    eta = stretch_inner_only(nr, beta_r)
    r = r_inner + (r_outer - r_inner) * eta
    span = theta1 - theta0
    full_circle = np.isclose(abs(span), 2.0 * np.pi, rtol=0.0, atol=1.0e-12)
    endpoint = not full_circle
    th = np.linspace(theta0, theta1, ntheta, endpoint=endpoint)

    # Radial halos: linear extension in r.
    r_h = extend_1d_linear(r, nhalo)

    # Angular halos: wrap for full circle, linear extension in theta otherwise.
    if full_circle:
        dth = span / ntheta
        th_h = theta0 + dth * np.arange(-nhalo, ntheta + nhalo)
    else:
        th_h = extend_1d_linear(th, nhalo)

    rr, tt = np.meshgrid(r_h, th_h, indexing="ij")
    xh = rr * np.cos(tt)
    yh = rr * np.sin(tt)
    return xh, yh


def write_opensbli_data(output_path, xh, yh, npoints, nhalo, transpose=False):
    if transpose:
        xh = xh.T
        yh = yh.T

    tmp = output_path.with_suffix(".tmp.h5")
    with h5py.File(tmp, "w") as f:
        f.create_dataset("x0", data=xh)
        f.create_dataset("x1", data=yh)

    b = SimulationBlock(2, block_number=0)
    with h5py.File(tmp, "r") as f:
        arrays = [f["x0"], f["x1"]]
        array_names = ["x0", "x1"]
        if transpose:
            npoints = [npoints[1], npoints[0]]
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
        default=1.2,
        help="Radial clustering exponent (>1 clusters near inner wall only)",
    )
    parser.add_argument("--nhalo", type=int, default=5, help="Halo width")
    parser.add_argument(
        "--transpose",
        action="store_true",
        help="Transpose mesh arrays before writing OpenSBLI datasets.",
    )
    parser.add_argument("--output", default="data.h5", help="Output HDF5 path")
    args = parser.parse_args()

    if args.r_outer <= args.r_inner:
        raise ValueError("--r-outer must be > --r-inner")
    if args.theta1 <= args.theta0:
        raise ValueError("--theta1 must be > --theta0")

    xh, yh = build_annulus_with_halos(
        args.nr, args.ntheta, args.r_inner, args.r_outer, args.theta0, args.theta1, args.beta_r, args.nhalo
    )
    out = Path(args.output).resolve()
    write_opensbli_data(out, xh, yh, [args.nr, args.ntheta], args.nhalo, transpose=args.transpose)
    print(f"Wrote {out}")
    print(f"Interior size: ({args.nr}, {args.ntheta}), halo={args.nhalo}")
    print(f"Transposed: {args.transpose}")
    print(f"Theta range: [{args.theta0}, {args.theta1}] rad")


if __name__ == "__main__":
    main()
