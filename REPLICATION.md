# Reproducing the preregistered replication

This file lists how to rerun **Mixing and Switching Between Models of the Iowa Gambling Task**, preregistered at
<https://doi.org/10.17605/OSF.IO/U86SY>, from the raw data to every number reported. The scripts for the earlier
work (scripts 01–42) are described under *How to run* in the README.

## What was committed when

| Commit | What | When |
|---|---|---|
| `ebdddfe` | Data, simulation and cross-validation scripts (43–45), committed before they were run | Nov 20, 2026 |
| `8c4dcd2` | Stacking, switching and scoring scripts (46–48), committed before they were run | Nov 21, 2026 |
| `a59e7b5` | Results of the preregistered analyses (`output/prereg/`) | Nov 21, 2026 |
| `65d4d8d`, `cf7f329` | Exploratory: simulated switchers and one-model players (44b), weight recovery (49) | Nov 21, 2026 |
| `a0fcc39` | Exploratory: results of the simulations | Nov 22, 2026 |

Deviation from the registration: the registration planned to run the cross-validation on Baylor's Kodiak cluster.
It was run on an Apple silicon MacBook Air instead, because the cluster could not be reached off campus. The code
and seeds are the same; `kodiak/prereg_cv.sh` still runs it on a PBS cluster.

## Setup

1. Install R 4.6.1, open the project and run `renv::restore()`. The analyses used **loo 2.10.1**. Use the
   versions in `renv.lock`: each person's stacking weights are fit to about 90 choices, and other loo versions or
   platforms can stop the optimizer at slightly different weights. Scripts 46 print and save the versions they ran with.
2. Download `IGTdata.rdata` from <https://osf.io/8t7rm> (Steingroever et al., 2015; CC BY-SA 4.0) and put it in
   `data-raw/`. The data are not in this repository.

Run every command from the project folder.

## The preregistered analyses

| Step | Script | What it does | Time* |
|---|---|---|---|
| 1 | `44_null_simulation.R` | Simulated players who mix the eight models and never switch (N = 324, seed 2027) | minutes |
| 2 | `45_prereg_cv.R` | Ten-fold cross-validation of each model, one model, fold and quarter of the participants per run (280 runs per data set) | about 7.5 h (replication) and 6 h (simulated), 6 runs at a time |
| 3 | `46_prereg_stacking.R` | Held-out predictions; shared, personal and shrinkage mixes; hedged single models | about 2 min per data set |
| 4 | `47_prereg_switching.R` | Switching, shuffled, drifting weights, personal switching | about 3–10 min per data set |
| 5 | `48_prereg_score.R` | Scores every prediction, then the robustness checks and exploratory descriptions | seconds |

*On an Apple silicon MacBook Air.

Steps 1–2, all runs in parallel and resumable (a run whose output already exists is skipped):

```
Rscript scripts/44_null_simulation.R
for d in replication simulated; do for m in orl variant orl_alpha variant_alpha vse vse_split eef; do
  for f in 01 02 03 04 05 06 07 08 09 10; do for p in 1 2 3 4; do
    [ -f output/prereg/cv/$d/${m}_fold${f}_part$p.rds ] || echo "$d $m $f $p"
  done; done; done; done | xargs -P 6 -n 4 sh -c 'Rscript scripts/45_prereg_cv.R $0 $1 $2 $3 4'
```

On a PBS cluster: `qsub -v DATA=replication kodiak/prereg_cv.sh` (and `DATA=simulated`).

Steps 3–5:

```
for d in replication simulated replication342; do Rscript scripts/46_prereg_stacking.R $d; done
for d in replication simulated replication342; do Rscript scripts/47_prereg_switching.R $d; done
Rscript scripts/48_prereg_score.R
```

`replication342` (robustness check 5) reuses the replication sample's cross-validation fits. The results are in
`output/prereg/`; `results.rds` holds every test.

Built-in checks: script 46 prints the largest difference between the refit held-out predictions and the totals saved
by script 45, and the difference between the λ = 0 shrinkage mix and the shared mix (both should be 0); script 47
prints its optimizer codes (0 = converged). Running scripts 46 and 47 with `horstmann` reproduces the exploratory
Horstmann results of scripts 36–38, 42 and 42b.

## Exploratory analyses

| Script | What it does | Time* |
|---|---|---|
| `44b_switch_simulation.R` | Simulated players with known switching: one switching at the Horstmann rate (`sim_switch`), one keeping one model all session (`sim_types`); N = 162 each | under a minute |
| `45`–`47` with `sim_switch` and `sim_types` | The same pipeline on those players | about 6 h, then minutes |
| `49_weight_recovery.R` | Compares the personal weights estimated for the no-switching players with the true weights | seconds |
| `51_switch_triggers.R` | When the switching model infers a switch (after gains vs. losses), and whether a random-choice state explains part of the switching gain | minutes |
| `50_hybrid_model.R` | The hybrid model (ORL + α with a starting value, decay and separate bending); not yet fit | — |

```
Rscript scripts/44b_switch_simulation.R
# then the cross-validation loop above with: for d in sim_switch sim_types
for d in sim_switch sim_types; do Rscript scripts/46_prereg_stacking.R $d; Rscript scripts/47_prereg_switching.R $d; done
Rscript scripts/49_weight_recovery.R
for d in replication sim_switch; do Rscript scripts/51_switch_triggers.R $d; done
```
