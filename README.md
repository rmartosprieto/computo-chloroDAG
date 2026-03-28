# Causal inference of post-transcriptional regulation timelines from
long-read sequencing in Arabidopsis thaliana
Rubén Martos, Christophe Ambroise, Guillem Rigaill
2026-03-28

*We present a framework for reconstructing the chronology of chloroplast
gene regulation in Arabidopsis thaliana using causal inference,
addressing missing data challenges, and generating testable models for
ndhB and ndhD genes.*

[![build and
publish](https://github.com/rmartosprieto/computo-chloroDAG/actions/workflows/build.yml/badge.svg)](https://github.com/rmartosprieto/computo-chloroDAG/actions/workflows/build.yml)
[![Creative Commons
License](https://i.creativecommons.org/l/by/4.0/80x15.png)](http://creativecommons.org/licenses/by/4.0/)

### Authors

- [Rubén Martos]() (Université d’Évry Paris-Saclay)
- [Christophe Ambroise]() (Université d’Évry Paris-Saclay)
- [Guillem Rigaill]() (Université Paris-Saclay / CNRS / INRAE /
  Université d’Évry)

### Abstract

We propose a novel framework for reconstructing the chronology of
genetic regulation using causal inference based on Pearl’s theory. The
approach proceeds in three main stages: causal discovery, causal
inference, and chronology construction. We apply it to the and genes of
the chloroplast in , generating four alternative maturation timeline
models per gene, each derived from a different causal discovery
algorithm (HC, PC, LiNGAM, or NOTEARS). Two methodological challenges
are addressed: the presence of missing data, handled via an EM algorithm
that jointly imputes missing values and estimates the Bayesian network,
and the selection of the $\ell_1$-regularization parameter in NOTEARS,
for which we introduce a stability selection strategy. The resulting
causal models consistently outperform reference chronologies in terms of
both reliability and model fit. Moreover, by combining causal reasoning
with domain expertise, the framework enables the formulation of testable
hypotheses and the design of targeted experimental interventions
grounded in theoretical predictions.
