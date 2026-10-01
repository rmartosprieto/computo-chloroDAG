#!/usr/bin/env python3
"""
Regenerate the manuscript's scoring tables on a COMMON variable set.

The published `Tables_Scoring/BIC_scoring_summary_*.csv` /
`LogL_scoring_summary_*.csv` are produced with `DiscreteBayesianNetwork(edges)`,
which silently drops isolated events; models are then compared on different
variable sets. This script rewrites the four CSVs using the scores computed by
`run_scores_imputed.py` on the full, common variable set (isolated events
included), keeping the exact column schema used by the manuscript figures:

    ;Reference model;Causal model;Scores;Difference;Best BIC
    ;Reference model;Causal model;Scores;Difference;Best Log-Likelihood

Run:  python comparisons/ablation/make_scoring_tables.py
"""
from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = HERE / "results" / "ablation_scores_imputed.csv"
DEST = ROOT / "Tables_Scoring"
METHODS = ["HC", "PC", "LiNGAM", "NOTEARS"]


def main():
    with open(SRC) as fh:
        rows = list(csv.DictReader(fh))
    for gene in ("ndhB", "ndhD"):
        ref_bic = float([r for r in rows if r["dataset"] == gene
                         and r["method"] == "Reference"][0]["BIC"])
        ref_logl = float([r for r in rows if r["dataset"] == gene
                          and r["method"] == "Reference"][0]["logL"])
        for metric, ref, fname, best_col in (
            ("BIC", ref_bic, f"BIC_scoring_summary_{gene}.csv", "Best BIC"),
            ("logL", ref_logl, f"LogL_scoring_summary_{gene}.csv", "Best Log-Likelihood"),
        ):
            lines = [["", "Reference model", "Causal model", "Scores", "Difference", best_col]]
            for i, m in enumerate(METHODS):
                score = float([r for r in rows if r["dataset"] == gene
                               and r["method"] == m and r["variant"] == "step2_chron"][0][metric])
                best = "Reference Model" if ref > score else m
                lines.append([i, ref, m, score, abs(score - ref), best])
            with open(DEST / fname, "w", newline="") as fh:
                csv.writer(fh, delimiter=";").writerows(lines)
            print(f"wrote {DEST / fname}")


if __name__ == "__main__":
    main()
