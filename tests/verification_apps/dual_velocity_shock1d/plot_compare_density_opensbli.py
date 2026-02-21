#!/usr/bin/env python3
import argparse
from pathlib import Path

import h5py
import numpy as np


def _format_ma_for_file(ma_val: float) -> str:
    return f"{ma_val:.2f}".replace(".", "_")


def _read_dataset(group, name):
    ds = group[name]
    d_m = ds.attrs.get("d_m", None)
    if d_m is None:
        return ds[()]
    shape = ds.shape
    read_start = [abs(int(d)) for d in d_m]
    read_end = [int(s - abs(int(d))) for d, s in zip(d_m, shape)]
    if len(read_end) == 1:
        return ds[read_start[0] : read_end[0]]
    if len(read_end) == 2:
        return ds[read_start[0] : read_end[0], read_start[1] : read_end[1]]
    if len(read_end) == 3:
        return ds[
            read_start[0] : read_end[0],
            read_start[1] : read_end[1],
            read_start[2] : read_end[2],
        ]
    raise RuntimeError(f"Unsupported dataset rank for {name}: {len(read_end)}")


def _load_opensbli_density(h5_path: Path):
    with h5py.File(h5_path, "r") as f:
        g = f["opensbliblock00"]
        x = np.asarray(_read_dataset(g, "x0_B0")).reshape(-1)
        rho = np.asarray(_read_dataset(g, "rho_B0")).reshape(-1)
    order = np.argsort(x)
    return x[order], rho[order]


def _normalize_center(xs, rhos):
    rho1 = rhos[0]
    rho2 = rhos[-1]
    rho_n = (rhos - rho1) / (rho2 - rho1)
    mid_idx = np.argmin(np.abs(rho_n - 0.5))
    x_center = xs[mid_idx]
    return xs - x_center, rho_n


def _load_experiment(path: Path, merge_eps: float, bin_width: float):
    data = np.loadtxt(path)
    xs = data[:, 0]
    rhos = data[:, 1]
    order = np.argsort(xs)
    xs = xs[order]
    rhos = rhos[order]

    if xs.size <= 1:
        return xs, rhos

    if bin_width > 0:
        bins = np.round(xs / bin_width) * bin_width
        unique_bins = {}
        for x, rho, b in zip(xs, rhos, bins):
            unique_bins.setdefault(b, []).append((x, rho))
        merged_x = []
        merged_rho = []
        for b in sorted(unique_bins.keys()):
            vals = unique_bins[b]
            merged_x.append(np.mean([v[0] for v in vals]))
            merged_rho.append(np.mean([v[1] for v in vals]))
        return np.array(merged_x), np.array(merged_rho)

    if merge_eps <= 0:
        diffs = np.diff(xs)
        med = np.median(diffs)
        merge_eps = max(1e-6, 1.25 * med)

    merged_x = []
    merged_rho = []
    i = 0
    n = xs.size
    while i < n:
        j = i + 1
        while j < n and abs(xs[j] - xs[i]) <= merge_eps:
            j += 1
        merged_x.append(xs[i:j].mean())
        merged_rho.append(rhos[i:j].mean())
        i = j
    return np.array(merged_x), np.array(merged_rho)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ma", type=float, default=1.55, help="Mach number for selecting experimental file")
    parser.add_argument("--Ml", default=None, help="Ml label for legend/title")
    parser.add_argument("--opensbli", default="build/opensbli_output.h5", help="OpenSBLI output h5 path")
    parser.add_argument("--opensbli2", default=None, help="Second OpenSBLI output h5 path (optional)")
    parser.add_argument("--label2", default="OpenSBLI NSF ($\\mathrm{M}_\\ell=0$)", help="Legend label for --opensbli2")
    parser.add_argument("--exp", default=None, help="Experimental data txt path")
    parser.add_argument("--exp-merge-eps", type=float, default=-1.0)
    parser.add_argument("--exp-bin-width", type=float, default=0.0)
    parser.add_argument("--out", default=None, help="Output PDF path")
    args = parser.parse_args()

    ma_file = _format_ma_for_file(args.Ma)
    exp_path = Path(args.exp) if args.exp else Path(f"experimental_shock_{ma_file}.txt")
    if not exp_path.exists():
        raise SystemExit(f"Experimental file not found: {exp_path}")

    h5_path = Path(args.opensbli)
    if not h5_path.exists():
        raise SystemExit(f"OpenSBLI output not found: {h5_path}")

    sim_x, sim_rho = _load_opensbli_density(h5_path)
    sim_x, sim_rho_n = _normalize_center(sim_x, sim_rho)

    sim2_x = None
    sim2_rho_n = None
    if args.opensbli2:
        h5_path2 = Path(args.opensbli2)
        if not h5_path2.exists():
            raise SystemExit(f"Second OpenSBLI output not found: {h5_path2}")
        sim2_x, sim2_rho = _load_opensbli_density(h5_path2)
        sim2_x, sim2_rho_n = _normalize_center(sim2_x, sim2_rho)

    exp_x, exp_rho = _load_experiment(exp_path, args.exp_merge_eps, args.exp_bin_width)
    exp_x, exp_rho_n = _normalize_center(exp_x, exp_rho)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
        }
    )

    fig, ax = plt.subplots()
    if args.Ml is None:
        label = "OpenSBLI Dual-velocity"
    else:
        label = rf"OpenSBLI Dual-velocity ($\mathrm{{M}}_\ell={args.Ml}$)"
    ax.plot(sim_x, sim_rho_n, linewidth=1.5, label=label)
    if sim2_x is not None:
        ax.plot(sim2_x, sim2_rho_n, linewidth=1.5, linestyle="--", label=args.label2)
    ax.scatter(exp_x, exp_rho_n, s=12, marker="o", label="Experiment", zorder=3)
    ax.set_xlabel(r"$x/\lambda$")
    ax.set_ylabel(r"$\varrho$")
    title = rf"Shock density profile ($\mathrm{{Ma}}={args.Ma:g}$)"
    ax.set_title(title)
    ax.set_xlim(-10.0, 10.0)
    ax.set_xticks(np.arange(-10.0, 10.1, 2.5))
    ax.legend()
    fig.tight_layout()

    if args.out:
        out_path = Path(args.out)
    else:
        if args.Ml is None:
            out_path = Path(f"shock_compare_M{ma_file}_opensbli_vs_exp.pdf")
        else:
            out_path = Path(f"shock_compare_M{ma_file}_Ml{args.Ml}_opensbli_vs_exp.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
