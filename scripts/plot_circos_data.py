#!/usr/bin/env python3
"""Create a circular chromosome-level bar plot from chromosome and value columns."""
import argparse
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_values(path):
    labels, values = [], []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chromosome", "value"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Input TSV must have columns: chromosome, value")
        for row in reader:
            label = row["chromosome"].strip()
            try:
                value = float(row["value"])
            except (TypeError, ValueError):
                raise ValueError("value must be numeric") from None
            if not label or not np.isfinite(value) or value < 0:
                raise ValueError("chromosome labels must be non-empty and values finite/non-negative")
            labels.append(label)
            values.append(value)
    if not labels:
        raise ValueError("No chromosome values found")
    if len(set(labels)) != len(labels):
        raise ValueError("Chromosome labels must be unique")
    return labels, np.asarray(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="TSV with chromosome and value columns")
    parser.add_argument("--out", default="plots/circos_mutation_plot.png")
    parser.add_argument("--demo", action="store_true", help="Use synthetic example data")
    args = parser.parse_args()
    if not args.input and not args.demo:
        parser.error("provide --input chromosome_values.tsv or explicitly use --demo")
    if args.input and args.demo:
        parser.error("--input and --demo cannot be used together")
    if args.demo:
        labels = [f"chr{i}" for i in range(1, 13)]
        values = np.random.default_rng(12).integers(50, 300, size=len(labels))
    else:
        labels, values = read_values(args.input)
    theta = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    width = 2 * np.pi / max(len(labels), 1) * 0.8
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"polar": True})
    ax.bar(theta, values, width=width, alpha=0.85)
    for angle, value, label in zip(theta, values, labels):
        ax.text(angle, value + max(values.max() * 0.04, 0.1), label,
                ha="center", va="center", fontsize=8)
    ax.set_title("Chromosome-level values" + (" (synthetic demo)" if args.demo else ""), pad=25)
    ax.set_yticklabels([])
    ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Saved circular plot: {out}")


if __name__ == "__main__":
    main()
