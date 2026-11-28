# Hierarchical fits of the eight models: predictions

Written November 27, 2026, before any of these fits was run on real data, and committed together with the scripts that
run and score them (`scripts/55_stan_common.R`, `55_stan_fit.R`, `55c_stan_cv_score.R`, `55d_stan_predictions.R`,
`55e_stan_jobs.R`, `55_queue.sh`, `kodiak/stan_cv.sh`, `kodiak/stan_submit.sh`). Registered on OSF as an addendum to
the preregistered replication (<https://doi.org/10.17605/OSF.IO/U86SY>). Nothing registered in U86SY changes.

## What the hierarchical fits are

In U86SY each person was fit alone, by maximum likelihood ("unpooled"). Script 54 fit each person with a group prior
estimated by empirical Bayes ("pooled"). Here each model is fit to all participants at once in Stan, in the style of
hBayesDM (Ahn, Haines & Zhang, 2017): each person's parameters are drawn from a group distribution, and the group's
mean and spread are estimated along with them (`stan/hier/`, written by `scripts/55_make_stan.py`).

- Models: ORL, the ORL variant, ORL + α, variant + α, VSE, VSE-split, EEF and NOVA, with the equations and limits of
  the separate fits. WSLS is not refit.
- Priors: group means ~ Normal(0, 1); group spreads ~ half-Normal(0, 0.2) for bounded parameters, half-Cauchy(0, 1)
  for unbounded ones; non-centered; bounded parameters through Phi_approx.
- Sampler: 4 chains, 1,000 warm-up and 1,000 sampling iterations, adapt_delta 0.95, seed 2055 + fold.
- Folds: U86SY's ten. Held-out outcomes still update the model; held-out choices are left out of the likelihood, so
  the group distribution is estimated from training choices only. Each held-out choice is scored by its posterior
  predictive probability (the mean over draws of p(choice)).
- Samples: the replication sample (N = 440) decides; Horstmann (N = 162) is reported beside it.

## Predictions

Each test is a difference in held-out log-likelihood, summed over participants; SE = standard deviation of the
per-participant differences × √n; z = difference / SE. Each prediction is directional, so each test is one-sided at
α = .05 (the rule of U86SY), with no correction for multiple tests.

| | Prediction | Met if |
|---|---|---|
| S1a | Hierarchical NOVA beats unpooled NOVA | z > 1.645 |
| S1b | Hierarchical ORL + α beats unpooled ORL + α | z > 1.645 |
| S2 | Hierarchical NOVA beats hierarchical ORL + α | z > 1.645 |
| S3a | Hierarchical NOVA beats pooled NOVA | z > 1.645 |
| S3b | Hierarchical ORL + α beats pooled ORL + α | z > 1.645 |
| S4 | The shared mix of the eight with hierarchical members is worse than the shared mix of the eight unpooled | z < −1.645 |
| S5 | Hierarchical NOVA alone does worse than the shared mix of the eight unpooled | z < −1.645 |

The eight in S4 are the seven learning models plus WSLS, U86SY's mix; with hierarchical members, the seven are swapped
for their hierarchical fits and WSLS stays as it is. If S4 or S5 fails because the hierarchical side does better, that
is reported as a result, not hidden as a failed prediction.

Described, not tested: the same comparisons for the other six models, the group parameters (fits to all trials), and
the replication tests without the second record of each repeated pair (N = 438; nothing is refit, those two people's
scores are dropped).

## Convergence rule

A fit passes with no divergent transitions and no R-hat above 1.05 (group and person parameters).

1. A fit that does not pass is rerun once with 2,000 warm-up iterations and adapt_delta 0.99.
2. If the rerun still fails only on some people's parameters, the tests are reported twice: with everyone, and without
   the people with any R-hat above 1.05 in any fit of the models in that test. The version with everyone decides.
3. If, after the rerun, any fit of a model has a group-level R-hat above 1.05 or more than 1% divergent transitions
   (40 of 4,000 draws), that model's tests count as not tested: neither met nor failed. S4 counts as not tested if
   this happens to any of the seven models in its mix.

## What was known before

No fit of these Stan models had been run on real data. Known: every unpooled result of U86SY; the pooled fits and their
five predictions (`predictions/pooled_replication.md`, all met); the recovery of these Stan models on 60 simulated
people (`output/stan/recovery/`: Spearman .34–.92, no divergent transitions, largest R-hat 1.498 for NOVA and 1.528
for EEF); and paper 1's hierarchical fits to Horstmann, including the cross-validation of hierarchical ORL and the
variant, in which ORL did not fully converge in several folds (`notes/lab_notebook.md`).

## How to run

On the Mac, after this file and the scripts are committed (the queue refuses to start otherwise):

```
bash scripts/55_queue.sh start        # background; status, pause, resume, stop, kill
```

On Kodiak (after `git pull`; once: `bash kodiak/stan_setup.sh`), with any fits already done on the Mac copied into
`output/stan/` first:

```
bash kodiak/stan_submit.sh horstmann
bash kodiak/stan_submit.sh replication
```

Then, on the Mac with every fit in `output/stan/`:

```
for d in replication horstmann replication438; do
  Rscript scripts/55c_stan_cv_score.R $d > output/stan/log55c_$d.txt 2>&1
  Rscript scripts/55d_stan_predictions.R $d > output/stan/log55d_$d.txt 2>&1
done
```

Any change to the scripts or these predictions after the commit will be reported as a deviation.
