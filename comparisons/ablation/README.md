# Ablation of Step 2 (causal inference / pruning)

Standalone experiment answering the reviewer's request:

> *"Motivate and justify step 2 (causal inference) in your framework. Why do we
> not simply take the partial order of the DAG learnt in step 1? Why only keep
> the strongest effects and why do we want a tree? I strongly recommend
> performing an ablation study of the method without step 2."*

The appendix of the manuscript (`@sec-ablation-step2`, `@tbl-ablation-step2`)
reads the CSVs produced here; the rest is not wired into the Quarto build.

## What is compared

For each gene (`ndhB` without/with intron, `ndhD`) and each method
(HC, PC, LiNGAM, NOTEARS) we build graphs from the manuscript's own tables:

- **Step 1** — the *unpruned* graph: every candidate causal relation recorded in
  `Tables_Causal_Inference/causal_relations_*.csv`, with its estimated
  `Actual Effect` (before keeping only the strongest cause per outcome).
- **Step 2** — the pruned chronology stored in
  `Tables_Causal_Inference/chron_AdjMatrix_*.csv`.
- **Reference** — the Guilcher chronology (Fig6), read from
  `original-paper/Results/graphs/AdjacencyMatrix_chron_*_Fig6*.csv`.

We validate the reconstruction: keeping the single strongest cause per outcome
reproduces the stored chronology exactly for all 12 models, once the
**NOTEARS/ndhD** file is corrected. Originally that file contained a 3-edge graph
inconsistent with both the top-1 rule and the scored BIC in `Tables_Scoring`
(−11 127,3); it has been replaced by the correct 4-edge graph. A dedicated
script, `make_scoring_tables.py`, regenerates the corrected `Tables_Scoring/*`
on a common variable set.

Sensitivity variants are also built: `top2`, `top3` (keep the 2 or 3 strongest
causes per outcome) and `thr0.10`, `thr0.20` (keep all causes with
`Actual Effect >= threshold`).

## Metrics

- **`coverage`** — share of node pairs ordered by the graph (completeness of the
  partial order). Step 2 has lower coverage by construction.
- **`agreement_ref`** — among the pairs ordered by the *Guilcher reference
  chronology* and by the candidate graph, share ordered consistently. Read it
  together with `coverage` (it is computed on the intersection).
- **`BIC` / `logL`** — discrete Bayesian-network scores on the **exact imputed
  working datasets** used by the paper (`original-paper/Data/*_imputed_*.csv`),
  computed **by direct counting** (Python standard library only):
  `logL = sum_i log P(x_i)` (MLE) and `BIC = logL - 0.5 ln(N) k`. These values
  are identical to `pgmpy`'s `structure_score(..., "bic-d")` /
  `log_likelihood_score` (verified: maximum absolute difference 0.0), so `pgmpy`
  is not needed. Every DAG is scored on a **fixed, common node set** (isolated
  events included), which is required for a fair comparison.
- `run_scores.py` is the complete-case variant (N = 117 for ndhB, 930 for ndhD).

## Run

Python 3.11+, standard library only (no pandas, no pgmpy — runs in
`ComputoPython` as well as anywhere else):

```bash
python3 comparisons/ablation/run_ablation.py          # structure, coverage, agreement
python3 comparisons/ablation/run_scores_imputed.py    # BIC / logL on imputed data
python3 comparisons/ablation/run_scores.py            # complete-case fallback
python3 comparisons/ablation/make_scoring_tables.py   # regenerates Tables_Scoring/*
```

Outputs: `results/ablation_metrics.csv`, `results/ablation_scores_imputed.csv`,
`results/ablation_scores.csv`, `results/summary.md`. See `SYNTHESE.md` for the
reading.

## Caveats / findings

- **Isolated nodes were dropped by the paper's original scoring.** The
  manuscript built the Bayesian network from `DiscreteBayesianNetwork(edges)`, so
  events with no edge were omitted; models were then scored on different variable
  sets, which advantages the model with fewer variables. Example: the *ndhB*
  reference (12 events, 4 isolated) was scored on 8 variables in the old
  `Tables_Scoring` (−5 023,78); on 12 variables it gives −6 682,37. On a
  **common node set** the *ndhB* conclusion changes (the reference beats all four
  causal models). `run_scores_imputed.py` uses a common node set on purpose.
- **NOTEARS/ndhD** stored graph and scored model disagreed; the graph has been
  corrected (see above).
- The structural `agreement_ref` is data-independent; the BIC signal is not:
  complete cases favour Step 2 (11/12), the full imputed data favour Step 1
  (8/12). Do not over-interpret the fit comparison.
