#!/usr/bin/env python3
"""Create a metagene profile from a tab-separated gene-by-bin signal matrix."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_matrix(path):
    matrix = np.loadtxt(path, delimiter="\t", comments="#", ndmin=2)
    if matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[1] < 2:
        raise ValueError("Signal matrix must have at least one gene row and two bin columns.")
    if not np.isfinite(matrix).all():
        raise ValueError("Signal matrix contains non-finite values.")
    return matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="TSV numeric matrix: rows=genes, columns=normalized bins")
    parser.add_argument("--out", default="plots/metagene_plot.png", help="Output PNG path")
    parser.add_argument("--demo", action="store_true", help="Use synthetic example data (not biological data)")
    args = parser.parse_args()
    if not args.input and not args.demo:
        parser.error("provide --input MATRIX.tsv or explicitly use --demo")
    if args.input and args.demo:
        parser.error("--input and --demo cannot be used together")
    if args.demo:
        rng = np.random.default_rng(42)
        matrix = rng.random((50, 100)) * np.linspace(0.5, 1.5, 100)
    else:
        matrix = load_matrix(args.input)
    profile = matrix.mean(axis=0)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.plot(np.linspace(0, 100, profile.size), profile, linewidth=2)
    ax.set(title="Metagene profile" + (" (synthetic demo)" if args.demo else ""),
           xlabel="Normalized gene body (%)", ylabel="Mean signal")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Saved metagene plot: {out}")


if __name__ == "__main__":
    main()
