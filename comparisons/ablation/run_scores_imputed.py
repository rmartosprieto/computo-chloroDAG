#!/usr/bin/env python3
"""
Step-2 ablation scores on the manuscript's *imputed* working datasets.

Uses the exact datasets used by the paper:

    original-paper/Data/ndhB_imputed_IterativeImputer_LogisticReg.csv
    original-paper/Data/ndhB_imputed_IterativeImputer_LogisticReg_with_intron.csv
    original-paper/Data/ndhD_imputed_IterativeImputer_LogisticReg.csv

and reproduces the manuscript's scores (BIC and log-likelihood) **without any
third-party dependency** (Python standard library only), so that the study runs
inside the Computo environment (`ComputoPython`) as well as anywhere else.

The discrete Bayesian-network log-likelihood and BIC are computed by direct
counting:

    logL(M) = sum_i log P_M(x_i)        (maximum likelihood / empirical)
    BIC(M)  = logL(M) - 0.5 * ln(N) * k

with k the number of free parameters of the DAG. These values are identical to
`pgmpy`'s `structure_score(..., "bic-d")` / `log_likelihood_score` when the same
node set is used (verified: maximum absolute difference 0.0 over all models).

Every model is scored on a **fixed, common node set** (isolated events
included), which is required for a fair comparison; see `SYNTHESE.md`.

Run:  python3 comparisons/ablation/run_scores_imputed.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_ablation as ra  # noqa: E402
import run_scores as rs  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "original-paper" / "Data"
OUT = HERE / "results"

# dataset key -> (imputed file, intron column name or None)
IMPUTED = {
    "ndhB": ("ndhB_imputed_IterativeImputer_LogisticReg.csv", None),
    "ndhB_with_intron": ("ndhB_imputed_IterativeImputer_LogisticReg_with_intron.csv", "intron_2_16"),
    "ndhD": ("ndhD_imputed_IterativeImputer_LogisticReg.csv", None),
}


def _sites():
    with open(DATA / "edition_sep_col.csv") as fh:
        rows = list(csv.reader(fh))[1:]  # drop header line
    return ([rows[24 + k][3] for k in range(12)],   # ndhB (without intron)
            [rows[37 + k][3] for k in range(5)])    # ndhD


GP_B, GP_D = _sites()


def load_imputed(key):
    """Return (node_names, rows) with rows a list of {node: 0/1}."""
    fname, intron_col = IMPUTED[key]
    gene = "ndhB" if key.startswith("ndhB") else "ndhD"
    gp = GP_B if gene == "ndhB" else GP_D
    with open(DATA / fname) as fh:
        reader = csv.reader(fh)
        header = next(reader)
        raw = list(reader)
    cols = header[1:]  # drop the index column
    intron_pos = cols.index(intron_col) if intron_col is not None else None
    names = []
    for i, col in enumerate(cols):
        if intron_pos is not None and i == intron_pos:
            names.append("ndhB_intron")
        else:
            idx = i if intron_pos is None else (i if i < intron_pos else i - 1)
            names.append(f"{gene}_{gp[idx]}")
    rows = [{names[i]: int(r[i + 1]) for i in range(len(cols))} for r in raw]
    return names, rows


def main():
    out = []
    for gene, suffix, key in ra.DATASETS:
        nodes, data = load_imputed(key)
        for method in ra.METHODS:
            relations = ra.read_causal_relations(
                ra.CI / f"causal_relations_{method}_{gene}{suffix}.csv")
            _, chron_edges = ra.read_adj(
                ra.CI / f"chron_AdjMatrix_{method}_{gene}{suffix}.csv")
            step1 = {(a, b) for a, b, _ in relations}
            variants = [("step1_all", step1), ("step2_chron", chron_edges)]
            for k in (2, 3):
                variants.append((f"top{k}", ra.keep_topk(relations, k)))
            for thr in (0.10, 0.20):
                variants.append((f"thr{thr:.2f}", ra.keep_threshold(relations, thr)))
            for name, edges in variants:
                e = ra.restrict(set(nodes), edges)
                if not e:
                    continue
                logL, _, bic = rs.count_bic(data, list(nodes), e)
                out.append(dict(dataset=key, method=method, variant=name,
                                n_edges=len(e), n_reads=len(data),
                                BIC=round(bic, 3), logL=round(logL, 3)))
                print(f"{key:16s} {method:8s} {name:12s} edges={len(e):3d} "
                      f"BIC={bic:11.3f} logL={logL:11.3f}")
        _, ref_edges = ra.read_adj(ra.REF_DIR / ra.REFERENCE_FILE[key])
        logL, _, bic = rs.count_bic(data, list(nodes), ref_edges)
        out.append(dict(dataset=key, method="Reference", variant="reference",
                        n_edges=len(ref_edges), n_reads=len(data),
                        BIC=round(bic, 3), logL=round(logL, 3)))
        print(f"{key:16s} {'Reference':8s} {'reference':12s} edges={len(ref_edges):3d} "
              f"BIC={bic:11.3f} logL={logL:11.3f}")

    with open(OUT / "ablation_scores_imputed.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"\nWrote {OUT/'ablation_scores_imputed.csv'}")


if __name__ == "__main__":
    main()
