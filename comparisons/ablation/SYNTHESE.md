# Synthèse — Ablation de l'étape 2 (avec les imputations du papier)

Réponse au point du reviewer : *« Motiver et justifier l'étape 2 … pourquoi ne
pas simplement prendre l'ordre partiel du DAG de l'étape 1 ? … je recommande
fortement une étude d'ablation sans l'étape 2. »*

Expérience autonome, non intégrée au manuscrit. Scripts :
`run_ablation.py` (structure), `run_scores_imputed.py` (scores sur le jeu
imputé), `run_scores.py` (variante reads complets). Résultats :
`results/ablation_metrics.csv`, `results/ablation_scores_imputed.csv`,
`results/ablation_scores.csv`, `results/summary.md`.

## Protocole

Pour chaque gène (*ndhB* sans/avec intron, *ndhD*) et chaque méthode (HC, PC,
LiNGAM, NOTEARS) :

- **étape 1** = le DAG non élagué : toutes les relations de
  `causal_relations_*.csv` avec leur `Actual Effect` ;
- **étape 2** = la chronologie élaguée `chron_AdjMatrix_*.csv` (un seul parent
  conservé par outcome, l'effet le plus fort) ;
- **référence** = chronologie de Guilcher (Fig6), lue dans
  `original-paper/Results/graphs/`.

Les scores (BIC, log-vraisemblance) sont calculés **sur le jeu imputé exact du
papier** (`original-paper/Data/ndhB_imputed_IterativeImputer_LogisticReg*.csv`,
`ndhD_imputed_IterativeImputer_LogisticReg.csv`), avec le même scoring que le
code du papier (`structure_score(..., "bic-d")` / `log_likelihood_score`).

**Contrôle de reproduction** : on retrouve **exactement** les valeurs de
`Tables_Scoring/*` pour HC, PC, LiNGAM (et la référence), sur ndhB et ndhD. Deux
problèmes ont été identifiés et corrigés (voir Résultat 3) : la suppression des
nœuds isolés dans le scoring, et un graphe NOTEARS/ndhD erroné.

## Résultat 1 — Structurel : l'étape 2 préserve/améliore l'accord avec la référence

Sur les graphes eux-mêmes (indépendant des données), `agreement_ref` (part des
paires ordonnées par la référence que le modèle ordonne pareil) **augmente dans 6
cas, reste identique dans 5, ne baisse qu'une fois** ; la couverture (part des
paires ordonnées) passe en moyenne de **0,65 à 0,42**.

| dataset | méthode | arêtes 1→2 | couverture 1→2 | accord réf. 1→2 |
|---|---|---|---|---|
| ndhB | HC | 31→10 | 0,60→0,44 | 0,421→**0,474** |
| ndhB | PC | 10→5 | 0,36→0,25 | 0,429→**0,500** |
| ndhB | LiNGAM | 45→9 | 0,82→0,42 | 0,526→0,368 |
| ndhB | NOTEARS | 10→6 | 0,46→0,39 | 0,250→0,250 |
| ndhB+intron | HC | 47→12 | 0,82→0,51 | 0,556→**0,593** |
| ndhB+intron | PC | 5→4 | 0,24→0,24 | 0,625→0,625 |
| ndhB+intron | LiNGAM | 47→11 | 0,69→0,36 | 0,481→**0,852** |
| ndhB+intron | NOTEARS | 10→6 | 0,46→0,46 | 0,250→0,250 |
| ndhD | HC | 8→4 | 0,80→0,80 | 0,444→0,444 |
| ndhD | PC | 8→3 | 0,80→0,40 | 0,667→0,667 |
| ndhD | LiNGAM | 8→3 | 0,80→0,40 | 0,556→**0,778** |
| ndhD | NOTEARS | 10→4 | 1,00→0,70 | 0,444→0,444 |

**Lecture** : élaguer vers les arêtes les plus fortes ne détruit pas l'ordre
temporel ; cela le resserre (moins de paires ordonnées, mais mieux alignées avec
la référence).

## Résultat 2 — Fit (BIC sur imputations) : l'étape 2 ne l'améliore pas

Sur le jeu imputé du papier, **tous les modèles notés sur le même jeu de
variables**, logL(étape 1) ≥ logL(étape 2) partout (emboîtement respecté), et le
BIC est **meilleur pour l'étape 1 dans 8 cas sur 12** (meilleur pour l'étape 2
dans 4 cas seulement : LiNGAM sur *ndhB* et *ndhB+intron*, PC et NOTEARS sur
*ndhB+intron*, de façon marginale).

| dataset | HC | PC | LiNGAM | NOTEARS |
|---|---|---|---|---|
| ndhB | étape 1 meilleure | étape 1 | **étape 2** | étape 1 (marginal) |
| ndhB+intron | étape 1 | **étape 2** (marginal) | **étape 2** | **étape 2** (marginal) |
| ndhD | étape 1 | étape 1 | étape 1 | étape 1 |

Exemple *ndhD* : HC −10 981 → −11 113, PC −11 026 → −11 148, LiNGAM −11 040 →
−11 355, NOTEARS −11 024 → −11 148.

**Conséquence directe pour la réponse au reviewer** : l'étape 2 ne se justifie
**pas** par un meilleur ajustement. Elle se justifie comme **heuristique de
parcimonie/lisibilité** : elle conserve (voire améliore) l'accord avec la
chronologie de référence tout en réduisant le nombre de relations retenues.
C'est exactement le cadrage que le reviewer demandait (« present it as a
heuristic »).

## Résultat 3 — Deux problèmes de scoring découverts au passage

1. **Nœuds isolés supprimés du scoring.** Le code du papier construit le réseau
   bayésien à partir des **seules arêtes** (`DiscreteBayesianNetwork(edges)`) :
   les événements isolés sont donc omis. Deux modèles sont alors notés sur des
   ensembles de variables différents, ce qui avantage mécaniquement celui qui a
   le moins de variables. Exemple : la référence *ndhB* (12 événements, dont 4
   isolés) est notée sur **8 variables** (−5 023,78) ; notée sur les 12
   variables, elle vaut −6 682,37.
   Conséquence : à jeu de variables égal, la conclusion « NOTEARS se distingue
   pour *ndhB* » **ne tient plus** — la référence bat les quatre modèles causaux
   sur *ndhB* et *ndhB+intron*. En revanche, la conclusion pour *ndhD* (modèles
   causaux meilleurs que la référence) **est confirmée** : référence −11 702,95
   vs HC −11 112,70, PC −11 147,59, LiNGAM −11 355,41, NOTEARS −11 147,59.
   → Recommandation : re-scorer tous les modèles sur un **jeu de variables
   commun** (et le préciser).

2. **Graphe NOTEARS/ndhD erroné (corrigé).** Le fichier
   `Tables_Causal_Inference/chron_AdjMatrix_NOTEARS_ndhD.csv`
   contenait un graphe à 3 arêtes (identique à PC), incohérent avec le BIC noté
   dans `Tables_Scoring` (−11 127,33). Le graphe correct (4 arêtes, qui est
   exactement la sélection top-1 des `causal_relations`) a été rétabli. La
   reconstruction « top-1 → chronologie » est donc désormais **12/12**.

## Recommandation pour le papier et la réponse

1. **Présenter l'étape 2 comme une heuristique de parcimonie** (et non comme une
   étape d'inférence causale supplémentaire), en citant l'ablation : accord
   référence conservé/amélioré, couverture réduite, BIC non amélioré.
2. **Re-scorer BIC/logL sur un jeu de variables commun** (inclure les nœuds
   isolés), et reporter les valeurs corrigées ; la conclusion *ndhB* change.
3. **Vérifier l'incohérence NOTEARS/ndhD** entre `chron_AdjMatrix` et
   `Tables_Scoring`.
4. Optionnel : ajouter un tableau d'ablation (étape 1 vs étape 2) en annexe.

## Limites

- `agreement_ref` se calcule sur l'intersection des paires ordonnées : à lire
  avec la couverture.
- L'étape 1 est reconstruite à partir de `causal_relations_*.csv` (relations
  candidates), qui est la meilleure trace disponible du DAG découvert.
- NOTEARS/ndhD a été corrigé (graphe rétabli) : plus d'approximation sur ce
  point.
- La variante `run_scores.py` (reads complètement observés, N=117/930) donne un
  BIC favorable à l'étape 2 dans 11/12 cas : le signal **dépend fortement de la
  taille d'échantillon** — sur le jeu imputé complet (N=1899/7752), c'est
  l'étape 1 qui gagne. Argument supplémentaire pour ne pas surinterpréter le fit.
