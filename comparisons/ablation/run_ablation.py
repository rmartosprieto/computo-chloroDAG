#!/usr/bin/env python3
"""
Ablation study for reviewer point "motivate/justify Step 2" and
"I strongly recommend performing an ablation study of the method without step 2".

Standalone: nothing here is wired into the Quarto manuscript.

Idea
----
The manuscript's pipeline is:
  Step 1  causal discovery -> a DAG (partial order of events)
  Step 2  causal inference -> keep the strongest causal relation per outcome,
          giving a pruned tree ("chronology")
  Step 3  chronology representation

The reviewer asks why Step 2 is needed at all, i.e. why not stop at the partial
order of Step 1. Here we reconstruct the *Step-1 graph* (all candidate causal
relations recorded in `causal_relations_*.csv`, with their estimated effects)
and compare it to the *Step-2 chronology* (`chron_AdjMatrix_*.csv`) against the
Guilcher reference chronology (Fig6).

We also probe sensitivity to the pruning rule: keeping the top-k strongest
causes per outcome (k = 1, 2, 3) and absolute-effect thresholds.

Conventions
-----------
- Edges are oriented `cause -> outcome` (same orientation as the manuscript's
  `causal_relations_*.csv` and `chron_AdjMatrix_*.csv`).
- A graph induces a partial order by transitive reachability. Two order metrics
  are reported:
    * `coverage`      = share of node pairs ordered by the graph (completeness);
    * `agreement_ref` = among the pairs ordered by the *reference* graph and by
                        the candidate graph, share ordered consistently.
  `agreement_ref` on a small intersection is the closest analogue of the
  "trajectory" accuracy, while `coverage` shows how much the graph commits to.
"""
from __future__ import annotations

import csv
import heapq
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PAPER = ROOT
CI = PAPER / "Tables_Causal_Inference"
REF_DIR = ROOT / "original-paper" / "Results" / "graphs"
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)

REFERENCE_FILE = {
    "ndhB":             "AdjacencyMatrix_chron_ndhB_Fig6_without_intron.csv",
    "ndhB_with_intron": "AdjacencyMatrix_chron_ndhB_Fig6.csv",
    "ndhD":             "AdjacencyMatrix_chron_ndhD_Fig6.csv",
}
DATASETS = [
    ("ndhB", "", "ndhB"),
    ("ndhB", "_with_intron", "ndhB_with_intron"),
    ("ndhD", "", "ndhD"),
]
METHODS = ["HC", "PC", "LiNGAM", "NOTEARS"]


# --------------------------------------------------------------------------- #
# io
# --------------------------------------------------------------------------- #
def read_adj(path):
    with open(path) as fh:
        reader = csv.reader(fh)
        header = next(reader)[1:]
        rows = list(reader)
    edges = {(row[0], header[j]) for row in rows
             for j, v in enumerate(row[1:]) if v not in ("", "nan", "0", "0.0")}
    return header, edges


def read_causal_relations(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        eff = r["Actual Effect"]
        eff = float(eff) if eff not in ("", "nan") else float("nan")
        out.append((r["Cause (X)"], r["Outcome (Y)"], eff))
    return out


# --------------------------------------------------------------------------- #
# graph utilities
# --------------------------------------------------------------------------- #
def restrict(nodes, edges):
    nodes = set(nodes)
    return {(a, b) for a, b in edges if a in nodes and b in nodes}


def trans_closure(nodes, edges):
    nodes = list(nodes)
    adj = {n: set() for n in nodes}
    for a, b in edges:
        if a in adj and b in adj:
            adj[a].add(b)
    out = set()
    for s in nodes:
        stack = list(adj[s])
        seen = set(stack)
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        for v in seen:
            out.add((s, v))
    return out


def topo_order(nodes, edges):
    nodes = [n for n in nodes]
    indeg = {n: 0 for n in nodes}
    adj = {n: [] for n in nodes}
    for a, b in restrict(nodes, edges):
        adj[a].append(b)
        indeg[b] += 1
    available = [n for n in nodes if indeg[n] == 0]
    heapq.heapify(available)
    order = []
    while available:
        u = heapq.heappop(available)
        order.append(u)
        for v in adj[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                heapq.heappush(available, v)
    return order if len(order) == len(nodes) else None


def metrics(nodes, edges, ref_pairs):
    """Return coverage, agreement_ref, n_ordered_pairs."""
    n = len(nodes)
    if n < 2:
        return 0.0, float("nan"), 0
    order = topo_order(nodes, edges)
    closure = trans_closure(nodes, edges)
    coverage = len(closure) / (n * (n - 1) / 2)
    if order is None:
        return coverage, float("nan"), len(closure)
    pos = {x: i for i, x in enumerate(order)}
    common = [(a, b) for (a, b) in ref_pairs if a in pos and b in pos]
    agreement = (sum(1 for a, b in common if pos[a] < pos[b]) / len(common)
                 if common else float("nan"))
    return coverage, agreement, len(closure)


# --------------------------------------------------------------------------- #
# pruning variants
# --------------------------------------------------------------------------- #
def keep_topk(relations, k):
    """Keep the k strongest causes per outcome (largest Actual Effect)."""
    by_outcome = {}
    for a, b, eff in relations:
        by_outcome.setdefault(b, []).append((eff, a))
    edges = set()
    for b, lst in by_outcome.items():
        lst = sorted(lst, key=lambda t: (-(t[0] if t[0] == t[0] else -1.0), t[1]))
        for e, a in lst[:k]:
            edges.add((a, b))
    return edges


def keep_threshold(relations, thr):
    return {(a, b) for a, b, eff in relations if eff == eff and eff >= thr}


def describe(edges):
    ns = sorted({x for e in edges for x in e})
    return len(ns), len(edges)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    rows = []
    for gene, suffix, key in DATASETS:
        ref_nodes, ref_edges = read_adj(REF_DIR / REFERENCE_FILE[key])
        ref_pairs = trans_closure(ref_nodes, ref_edges)
        for method in METHODS:
            relations = read_causal_relations(CI / f"causal_relations_{method}_{gene}{suffix}.csv")
            _, chron_edges = read_adj(CI / f"chron_AdjMatrix_{method}_{gene}{suffix}.csv")
            step1 = {(a, b) for a, b, _ in relations}
            nodes = [n for n in ref_nodes if n in {x for e in (step1 | chron_edges) for x in e}]

            variants = [("step1_all", step1), ("step2_chron", chron_edges)]
            for k in (1, 2, 3):
                variants.append((f"top{k}", keep_topk(relations, k)))
            for thr in (0.05, 0.10, 0.20):
                variants.append((f"thr{thr:.2f}", keep_threshold(relations, thr)))

            for name, edges in variants:
                e = restrict(nodes, edges)
                n_nodes, n_edges = describe(e)
                coverage, agreement, n_ordered = metrics(nodes, e, ref_pairs)
                chron_recall = (len(e & chron_edges) / len(chron_edges)
                                if chron_edges else float("nan"))
                rows.append(dict(
                    dataset=key, method=method, variant=name,
                    n_nodes=n_nodes, n_edges=n_edges,
                    n_ordered_pairs=n_ordered, coverage=round(coverage, 3),
                    n_ref_pairs=len([p for p in ref_pairs if p[0] in nodes and p[1] in nodes]),
                    agreement_ref=round(agreement, 3),
                    step2_edge_recall=round(chron_recall, 3),
                ))
                print(f"{key:16s} {method:8s} {name:12s} edges={n_edges:3d} "
                      f"cov={coverage:4.2f} agree_ref={agreement:5.3f} "
                      f"step2_recall={chron_recall:4.2f}")

    fields = list(rows[0].keys())
    with open(OUT / "ablation_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    write_summary(rows)
    print(f"\nWrote {OUT/'ablation_metrics.csv'} and {OUT/'summary.md'}")


def write_summary(rows):
    def sel(**kw):
        out = rows
        for k, v in kw.items():
            out = [r for r in out if r[k] == v]
        return out

    lines = ["# Ablation of Step 2 (structure only)", "",
             "Generated by `comparisons/ablation/run_ablation.py`.",
             "`agreement_ref` = share of reference-ordered pairs ordered consistently;",
             "`coverage` = share of node pairs ordered by the candidate graph.",
             "", "## Step 1 (all candidate edges) vs Step 2 (pruned chronology)", "",
             "| dataset | method | variant | edges | coverage | agreement_ref |",
             "|---|---|---|---|---|---|"]
    for gene, _, key in DATASETS:
        for method in METHODS:
            for variant in ("step1_all", "step2_chron"):
                r = sel(dataset=key, method=method, variant=variant)[0]
                lines.append(f"| {key} | {method} | {variant} | {r['n_edges']} "
                             f"| {r['coverage']} | {r['agreement_ref']} |")
    lines += ["", "## Sensitivity to the pruning rule (top-k / threshold)", "",
              "| dataset | method | variant | edges | coverage | agreement_ref | step2_edge_recall |",
              "|---|---|---|---|---|---|---|"]
    for gene, _, key in DATASETS:
        for method in METHODS:
            for r in sel(dataset=key, method=method):
                if r["variant"] in ("step1_all", "step2_chron"):
                    continue
                lines.append(f"| {key} | {method} | {r['variant']} | {r['n_edges']} "
                             f"| {r['coverage']} | {r['agreement_ref']} | {r['step2_edge_recall']} |")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
