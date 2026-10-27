#!/usr/bin/env python3
"""Plot signal profiles for ORFs from a TSV with orf, position, and signal columns."""
import argparse
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_profiles(path):
    profiles = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"orf", "position", "signal"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Input TSV must have columns: orf, position, signal")
        for row in reader:
            name = row["orf"].strip()
            try:
                position, signal = float(row["position"]), float(row["signal"])
            except (TypeError, ValueError):
                raise ValueError("position and signal must be numeric") from None
            if not name or not np.isfinite([position, signal]).all():
                raise ValueError("ORF names must be non-empty and values finite")
            profiles.setdefault(name, []).append((position, signal))
    if not profiles:
        raise ValueError("No ORF observations found")
    return {name: sorted(values) for name, values in profiles.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="TSV with orf, position, signal columns")
    parser.add_argument("--out", default="plots/metaorf_plot.png")
    parser.add_argument("--demo", action="store_true", help="Use synthetic example data")
    args = parser.parse_args()
    if not args.input and not args.demo:
        parser.error("provide --input ORF.tsv or explicitly use --demo")
    if args.input and args.demo:
        parser.error("--input and --demo cannot be used together")
    if args.demo:
        rng = np.random.default_rng(1)
        profiles = {f"ORF_{i+1}": list(zip(range(length), rng.random(length) * np.linspace(1, 2, length)))
                    for i, length in enumerate((150, 300, 120))}
    else:
        profiles = read_profiles(args.input)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for name, values in profiles.items():
        positions, signal = zip(*values)
        ax.plot(positions, signal, label=name, linewidth=1.5)
    ax.set(title="ORF signal profiles" + (" (synthetic demo)" if args.demo else ""),
           xlabel="Position in ORF", ylabel="Signal")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Saved metaORF plot: {out}")


if __name__ == "__main__":
    main()
