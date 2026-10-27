#!/usr/bin/env python3
"""Plot strand-aware average signal across gene bodies from a BigWig and BED6 genes.

Input BED columns: chrom, start, end, name, score, strand. Coordinates follow
BED's 0-based, half-open convention. The output is a gene-body metaprofile,
not a TSS-centered profile.
"""
import argparse
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_bed6(path):
    genes = []
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"BED line {line_number} must have at least 6 tab-separated columns")
            chrom, start, end, name, _score, strand = fields[:6]
            try:
                start, end = int(start), int(end)
            except ValueError:
                raise ValueError(f"BED line {line_number}: start/end must be integers") from None
            if not chrom or not name or start < 0 or end <= start:
                raise ValueError(f"BED line {line_number}: invalid chromosome, name, or interval")
            if strand not in ("+", "-"):
                raise ValueError(f"BED line {line_number}: strand must be '+' or '-'")
            genes.append((chrom, start, end, name, strand))
    if not genes:
        raise ValueError("No valid BED6 gene intervals were found")
    return genes


def summarize_bigwig(bigwig_path, genes_bed, bins=100, min_coverage=0.8):
    """Return per-gene binned signal matrix and labels; missing signal is NaN."""
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("pyBigWig is required; install it with python -m pip install pyBigWig") from exc
    genes = read_bed6(genes_bed)
    matrix, labels = [], []
    bw = pyBigWig.open(str(bigwig_path))
    if bw is None:
        raise RuntimeError(f"Could not open BigWig file: {bigwig_path}")
    try:
        chrom_sizes = bw.chroms()
        for chrom, start, end, name, strand in genes:
            if chrom not in chrom_sizes:
                continue
            chrom_length = chrom_sizes[chrom]
            start, end = max(0, start), min(end, chrom_length)
            if end <= start:
                continue
            edges = np.linspace(start, end, bins + 1, dtype=int)
            row = []
            for left, right in zip(edges[:-1], edges[1:]):
                if right <= left:
                    row.append(np.nan)
                    continue
                values = np.asarray(bw.values(chrom, int(left), int(right), numpy=True), dtype=float)
                valid = np.isfinite(values)
                if valid.size == 0 or valid.mean() < min_coverage:
                    row.append(np.nan)
                else:
                    row.append(float(np.mean(values[valid])))
            if strand == "-":
                row.reverse()
            if np.isfinite(row).any():
                matrix.append(row)
                labels.append(name)
    finally:
        bw.close()
    if not matrix:
        raise ValueError("No genes had usable signal. Check chromosome naming, BED coordinates, and BigWig genome build.")
    return np.asarray(matrix, dtype=float), labels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bigwig", required=True, help="BigWig signal track (.bw/.bigWig)")
    parser.add_argument("--genes-bed", required=True, help="BED6 gene intervals on the same genome build")
    parser.add_argument("--out", default="plots/bigwig_metagene.png", help="Output PNG path")
    parser.add_argument("--bins", type=int, default=100, help="Number of bins per gene body (default: 100)")
    parser.add_argument("--min-coverage", type=float, default=0.8, help="Minimum fraction of non-missing bases per bin (0-1)")
    parser.add_argument("--matrix-out", help="Optional TSV path for per-gene binned values")
    args = parser.parse_args()
    if args.bins < 2:
        parser.error("--bins must be at least 2")
    if not 0 <= args.min_coverage <= 1:
        parser.error("--min-coverage must be between 0 and 1")
    matrix, labels = summarize_bigwig(args.bigwig, args.genes_bed, args.bins, args.min_coverage)
    mean = np.nanmean(matrix, axis=0)
    n = np.sum(np.isfinite(matrix), axis=0)
    sem = np.divide(np.nanstd(matrix, axis=0, ddof=1), np.sqrt(n),
                    out=np.zeros_like(mean), where=n > 1)
    x = np.linspace(0, 100, args.bins)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(x, mean, linewidth=2, label=f"Mean signal (n={len(labels)} genes)")
    ax.fill_between(x, mean - sem, mean + sem, alpha=0.2, label="±1 SEM")
    ax.set(title="BigWig metagene profile", xlabel="Gene body (%)", ylabel="Mean BigWig signal")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    if args.matrix_out:
        matrix_path = Path(args.matrix_out)
        matrix_path.parent.mkdir(parents=True, exist_ok=True)
        with matrix_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene"] + [f"bin_{i+1}" for i in range(args.bins)])
            for name, row in zip(labels, matrix):
                writer.writerow([name] + ["" if not np.isfinite(v) else f"{v:.8g}" for v in row])
    print(f"Used {len(labels)} genes; saved plot: {out}")


if __name__ == "__main__":
    main()
