#!/usr/bin/env python3
"""Plot mutation impact scores from a TSV with domain and impact_score columns."""
import argparse
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_scores(path):
    domains, scores = [], []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"domain", "impact_score"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("Input TSV must have columns: domain, impact_score")
        for row in reader:
            domain = row["domain"].strip()
            try:
                score = float(row["impact_score"])
            except (TypeError, ValueError):
                raise ValueError("impact_score must be numeric") from None
            if not domain or not np.isfinite(score):
                raise ValueError("domain must be non-empty and scores finite")
            domains.append(domain)
            scores.append(score)
    if not domains:
        raise ValueError("No mutation impact scores found")
    return domains, np.asarray(scores)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="TSV with domain and impact_score columns")
    parser.add_argument("--out", default="plots/mutation_impact_plot.png")
    parser.add_argument("--demo", action="store_true", help="Use synthetic example data")
    args = parser.parse_args()
    if not args.input and not args.demo:
        parser.error("provide --input impact.tsv or explicitly use --demo")
    if args.input and args.demo:
        parser.error("--input and --demo cannot be used together")
    if args.demo:
        domains = ["N-term", "DBD", "Linker", "CTD", "Tail"]
        scores = np.random.default_rng(7).uniform(0, 100, len(domains))
    else:
        domains, scores = read_scores(args.input)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(domains, scores)
    ax.set(title="Mutation impact by domain" + (" (synthetic demo)" if args.demo else ""),
           ylabel="Impact score", xlabel="Protein domain")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Saved mutation impact plot: {out}")


if __name__ == "__main__":
    main()
