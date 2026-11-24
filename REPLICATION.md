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
| `50_hybrid_model.R` | The hybrid model: ORL + α with a starting value, decay and (in one version) separate bending | — |
| `50b_hybrid_cv.R` | Ten-fold cross-validation of the four hybrid versions, as script 45 (plus ORL + α's fit as one more start) | about 3–5 h (Horstmann), 6 runs at a time |
| `50c_hybrid_score.R` | The hybrid against ORL + α, the other models and the mixes | minutes |
| `50d_nova_replication.R` | NOVA (the hybrid) on the replication sample: four tests written down first (`predictions/nova_replication.md`) | minutes, after script 50b on `replication` (about 3 h) |
| `54_pooled_cv.R` | Every model refit with a group prior (empirical Bayes, expectation-maximization with a Laplace approximation), same folds | about 2 h (Horstmann), 6 runs at a time |
| `54b_pooled_score.R` | Pooled against unpooled fits, NOVA against ORL + α both pooled, and the mixes | minutes |
| `54c_pooled_predictions.R` | The pooled fits on the replication sample: five tests written down first (`predictions/pooled_replication.md`) | seconds, after script 54b |
| `56_identity_check.R` | VSE's exploration term equals ORL's perseverance with a negative weight, on every trial of both samples | seconds |
| `57_forward_scores.R` | Forward predictions for the Centaur comparison: each choice from a fit to the earlier trials only, under the other sample's group prior | about 5 seconds per person and model (NOVA about 10); everyone in both samples, four at a time, an hour or two |
| `57b_forward_compare.R` | Our models, their shared mixes and Centaur, forward, on the same people and choices | seconds |
| `58_duplicates.R` | Checks that the only repeated records in the data are Kjome 2/17 and Wood 298/299 (found when matching to Psych-101), and runs the registered tests on the 438 left without the second record of each pair, next to the 440 as registered | seconds, after scripts 46 and 47 on `replication438` |
| `58_duplicates_chain.sh` | Every step that combines participants, rerun without the repeated records (46, 47, 58, 50d, 54–54c, 57, 57b), one after another and resumable: `caffeinate -i bash scripts/58_duplicates_chain.sh` | about 6–7 h |
| `55_make_stan.py` | Writes the hierarchical Stan models (`stan/hier/`) from one copy of the likelihood functions | seconds |
| `55a_stan_check.R` | Every Stan model compiles, and its likelihood equals the R functions' on every trial | minutes |
| `55b_stan_recovery.R` | Parameter recovery for one Stan model (60 simulated players) | about 1 h per model |
| `55_stan_fit.R` | Hierarchical Bayesian fit of one model: all trials (`full`) or one cross-validation fold | about 1 h per fit (Horstmann) |
| `55_run_all.sh` | All 88 hierarchical fits of one data set on one machine (Mac or a rented computer), several at once | about 2 days (Horstmann, 8 cores) |
| `55c_stan_cv_score.R` | The hierarchical cross-validation against the unpooled and pooled fits, NOVA against ORL + α, the mixes, and each fit's diagnostics | about 10 min |
| `kodiak/stan_setup.sh`, `kodiak/stan_cv.sh` | The hierarchical fits on Kodiak: one-time setup, then 88 jobs per data set | about a day |
| `centaur/run_centaur.sh`, `kodiak/centaur.sh` | The Centaur benchmark from start to finish on Kodiak's GPUs or a rented GPU (`centaur/RENTING.md`) | about 1 h |
| `52_model_subsets.R` | The shared mix with fewer models: without one model at a time, and two sets of five | minutes |

```
Rscript scripts/44b_switch_simulation.R
# then the cross-validation loop above with: for d in sim_switch sim_types
for d in sim_switch sim_types; do Rscript scripts/46_prereg_stacking.R $d; Rscript scripts/47_prereg_switching.R $d; done
Rscript scripts/49_weight_recovery.R
for d in replication sim_switch; do Rscript scripts/51_switch_triggers.R $d; done
```

The hybrid model, on the discovery sample (script 46 with `horstmann` first, for the eight models' held-out predictions):

```
Rscript scripts/46_prereg_stacking.R horstmann
for v in orl_alpha_start orl_alpha_decay hybrid hybrid_sep; do for f in 01 02 03 04 05 06 07 08 09 10; do
  for p in 1 2; do [ -f output/hybrid/cv/horstmann/${v}_fold${f}_part$p.rds ] || echo "$v $f $p"; done
done; done | xargs -P 6 -n 3 sh -c 'Rscript scripts/50b_hybrid_cv.R horstmann $0 $1 $2 2 > /dev/null'
Rscript scripts/50c_hybrid_score.R horstmann
for d in horstmann replication; do Rscript scripts/52_model_subsets.R $d; done
```

NOVA on the replication sample (predictions in `predictions/nova_replication.md`, committed before running):

```
for f in 01 02 03 04 05 06 07 08 09 10; do for p in 1 2 3 4; do
  [ -f output/hybrid/cv/replication/hybrid_fold${f}_part$p.rds ] || echo "hybrid $f $p"
done; done | xargs -P 6 -n 3 sh -c 'Rscript scripts/50b_hybrid_cv.R replication $0 $1 $2 4 > /dev/null'
Rscript scripts/50d_nova_replication.R
```

Pooled fits (script 54), after NOVA's held-out predictions exist for that sample (script 50c or 50d):

```
for m in orl variant orl_alpha variant_alpha vse vse_split eef nova; do for f in 01 02 03 04 05 06 07 08 09 10; do
  [ -f output/pooled/cv/horstmann/${m}_fold$f.rds ] || echo "$m $f"
done; done | xargs -P 6 -n 2 sh -c 'Rscript scripts/54_pooled_cv.R horstmann $0 $1 > /dev/null'
Rscript scripts/54b_pooled_score.R horstmann
```

Pooled fits on the replication sample (predictions in `predictions/pooled_replication.md`, committed before running):

```
for m in orl variant orl_alpha variant_alpha vse vse_split eef nova; do for f in 01 02 03 04 05 06 07 08 09 10; do
  [ -f output/pooled/cv/replication/${m}_fold$f.rds ] || echo "$m $f"
done; done | xargs -P 6 -n 2 sh -c 'Rscript scripts/54_pooled_cv.R replication $0 $1 > /dev/null'
Rscript scripts/54b_pooled_score.R replication
Rscript scripts/54c_pooled_predictions.R replication
```

The correspondence between VSE's exploration and ORL's perseverance:

```
Rscript scripts/56_identity_check.R
```

Forward predictions and the Centaur comparison (after script 54 on both samples, and centaur/score_centaur.py):

```
for d in horstmann replication; do for m in orl variant orl_alpha variant_alpha vse vse_split eef nova wsls; do
  echo "$d $m"
done; done | xargs -P 6 -n 2 sh -c 'Rscript scripts/57_forward_scores.R $0 $1 who=centaur/psych101_igt_match.csv > /dev/null'
Rscript scripts/57b_forward_compare.R horstmann centaur/centaur70b_scores.csv
Rscript scripts/57b_forward_compare.R replication centaur/centaur70b_scores.csv
```

Hierarchical Bayesian versions (cmdstanr; the Stan files are written by script 55_make_stan.py):

```
python3 scripts/55_make_stan.py .
Rscript scripts/55a_stan_check.R
for m in orl orl_split orl_alpha orl_split_alpha nova vse vse_split eef; do Rscript scripts/55b_stan_recovery.R $m; done
Rscript scripts/55_stan_fit.R horstmann nova full
```

All hierarchical fits and their cross-validation, on Kodiak (or on one machine with `bash scripts/55_run_all.sh <data> <jobs>`):

```
bash kodiak/stan_setup.sh                  # once, on the login node
qsub -v DATA=horstmann kodiak/stan_cv.sh
qsub -v DATA=replication kodiak/stan_cv.sh
Rscript scripts/55c_stan_cv_score.R horstmann
Rscript scripts/55c_stan_cv_score.R replication
```

The Centaur scores (on Kodiak: `qsub kodiak/centaur.sh`; on a rented GPU: `centaur/RENTING.md`):

```
bash centaur/run_centaur.sh marcelbinz/Llama-3.1-Centaur-70B centaur/centaur70b_scores.csv
```
