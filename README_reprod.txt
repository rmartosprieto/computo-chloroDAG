How to compile these sources
============================

This folder is a self-contained Quarto project (the Computo extension is
bundled in _extensions/).

1. Create the Python environment and a Jupyter kernel named "computopython":

     conda env create -f environment.yml          # or micromamba
     conda activate ComputoPython
     python -m ipykernel install --user --name computopython --display-name ComputoPython

   (environment.yml provides jupyter, matplotlib, seaborn, numpy, pandas, networkx.)

2. Compile the manuscript:

     quarto render submission-chloroDAG.qmd --to computo-pdf

   The appendix table is computed at render time from
   comparisons/ablation/results/*.csv (included).


3. Regenerate the Step-2 ablation and the Bayesian scores (optional). Python
   3.11+, standard library only (no conda environment needed). Run from this
   folder:

     python3 comparisons/ablation/run_ablation.py         # structure, coverage, agreement
     python3 comparisons/ablation/run_scores_imputed.py   # BIC / logL on the imputed data
     python3 comparisons/ablation/run_scores.py           # complete-case variant
     python3 comparisons/ablation/make_scoring_tables.py  # rewrites Tables_Scoring/*

   Inputs: Tables_Causal_Inference/, Tables_Raw_Data/, the imputed datasets in
   original-paper/Data/ and the reference chronologies in
   original-paper/Results/graphs/ (all included). Outputs are written to
   comparisons/ablation/results/; the last script overwrites Tables_Scoring/*.
   The files shipped here are the outputs of these commands, so re-running them
   changes nothing. See comparisons/ablation/README.md (method, metrics) and
   SYNTHESE.md (reading of the results, in French).
