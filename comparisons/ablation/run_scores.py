#!/usr/bin/env python3
"""
BIC / log-likelihood arm of the Step-2 ablation.

Reuses the graph variants built by `run_ablation.py` and scores them on the
same complete-case data used by that script, with a direct, dependency-free
implementation of the discrete Bayesian-network log-likelihood and BIC:

    logL(M) = sum_i log P_M(x_i)        (maximum-likelihood / empirical)
    BIC(M)  = logL(M) - 0.5 * ln(N) * k

where k is the number of free parameters of the DAG (higher BIC is better,
matching the manuscript's convention).

Why not pgmpy
-------------
The manuscript scores models with `pgmpy`. With the current pgmpy (1.1.2) the
`BIC(df).score(model)` helper returns values that can exceed the true MLE
log-likelihood (e.g. a model with no edges scores 0.0), so it is not used here.
The count-based formulas above are the textbook definitions and are validated
by the nesting property logL(supergraph) >= logL(subgraph).

Proxy caveat
------------
The manuscript scores on its imputed working dataset, which is not stored in
this repository; our complete-case scoring cannot reproduce the absolute values
of `Tables_Scoring/*` (different N and different preprocessing). The ablation is
therefore an *internal* step1-vs-step2 comparison, not a replacement for the
paper's Table.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_ablation as ra  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = ROOT / "Tables_Raw_Data"
OUT = HERE / "results"

SITES = {
    "ndhB": dict(zip([f"ed{i}" for i in range(1, 13)],
                     [94622, 94999, 95225, 95608, 95644, 95650,
                      96419, 96439, 96457, 96579, 96698, 97016])),
    "ndhD": dict(zip([f"ed{i}" for i in range(1, 6)],
                     [116281, 116290, 116494, 116785, 117166])),
}
SITES["ndhB"]["intron"] = "intron"
DATA_FILE = {"ndhB": "ndhB", "ndhB_with_intron": "ndhB", "ndhD": "ndhD"}


def load_complete(gene):
    """List of dicts {node: 0/1} on fully observed reads (Err -> 0)."""
    with open(RAW / f"raw_{gene}.csv") as fh:
        reader = csv.reader(fh, delimiter=";")
        header = next(reader)
        raw = list(reader)
    cols = [c for c in header if c in SITES[gene]]
    names = {c: f"{gene}_{SITES[gene][c]}" for c in cols}
    rows = []
    for raw_row in raw:
        rec, ok = {}, True
        for c in cols:
            tok = raw_row[header.index(c)].strip()
            if tok == "" or tok == "nan":
                ok = False
                break
            rec[names[c]] = 0 if tok in ("False", "Err") else 1
        if ok:
            rows.append(rec)
    return rows


def count_bic(rows, nodes, edges):
    parents = {n: [] for n in nodes}
    for a, b in edges:
        if a in parents and b in parents:
            parents[b].append(a)
    N = len(rows)
    cards = {n: len({r[n] for r in rows}) for n in nodes}
    logL, k = 0.0, 0
    for node, pa in parents.items():
        if pa:
            tot = Counter()
            cnt = Counter()
            for r in rows:
                pk = tuple(r[p] for p in pa)
                tot[pk] += 1
                cnt[pk + (r[node],)] += 1
            for key, c in cnt.items():
                logL += c * math.log(c / tot[key[:-1]])
            k += (cards[node] - 1) * math.prod(cards[p] for p in pa)
        else:
            vc = Counter(r[node] for r in rows)
            logL += sum(v * math.log(v / N) for v in vc.values())
            k += cards[node] - 1
    return logL, k, logL - 0.5 * math.log(N) * k


def main():
    rows_out = []
    for gene, suffix, key in ra.DATASETS:
        data = load_complete(DATA_FILE[key])
        nodes_all = sorted({k for r in data for k in r})
        for method in ra.METHODS:
            relations = ra.read_causal_relations(
                ra.CI / f"causal_relations_{method}_{gene}{suffix}.csv")
            _, chron_edges = ra.read_adj(
                ra.CI / f"chron_AdjMatrix_{method}_{gene}{suffix}.csv")
            step1 = {(a, b) for a, b, _ in relations}
            nodes = {x for e in (step1 | chron_edges) for x in e}
            variants = [("step1_all", step1), ("step2_chron", chron_edges)]
            for k in (2, 3):
                variants.append((f"top{k}", ra.keep_topk(relations, k)))
            for thr in (0.10, 0.20):
                variants.append((f"thr{thr:.2f}", ra.keep_threshold(relations, thr)))
            for name, edges in variants:
                e = ra.restrict(nodes, edges)
                if not e:
                    continue
                logL, k, bic = count_bic(data, sorted(nodes), e)
                rows_out.append(dict(dataset=key, method=method, variant=name,
                                     n_edges=len(e), n_reads=len(data), params=k,
                                     BIC=round(bic, 1), logL=round(logL, 1)))
                print(f"{key:16s} {method:8s} {name:12s} edges={len(e):3d} "
                      f"BIC={bic:9.1f} logL={logL:9.1f}")
        # reference graph on the same reads
        _, ref_edges = ra.read_adj(ra.REF_DIR / ra.REFERENCE_FILE[key])
        n = sorted({x for e in ref_edges for x in e})
        logL, k, bic = count_bic(data, n, ref_edges)
        rows_out.append(dict(dataset=key, method="Reference", variant="reference",
                             n_edges=len(ref_edges), n_reads=len(data), params=k,
                             BIC=round(bic, 1), logL=round(logL, 1)))
        print(f"{key:16s} {'Reference':8s} {'reference':12s} edges={len(ref_edges):3d} "
              f"BIC={bic:9.1f} logL={logL:9.1f}")

    with open(OUT / "ablation_scores.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)
    print(f"\nWrote {OUT/'ablation_scores.csv'}")


if __name__ == "__main__":
    main()
