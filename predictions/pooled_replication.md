# Pooled fits on the replication sample: predictions

Written November 23, 2026, before the pooled fits were run on the replication sample, and committed together with the
script that scores them (`scripts/54c_pooled_predictions.R`). This is an exploratory follow-up to the preregistered
replication (<https://doi.org/10.17605/OSF.IO/U86SY>), not part of that registration.

## What pooling is

Script 54 refits every model with a group prior: empirical Bayes by expectation-maximization with a Laplace
approximation. Each person's parameters are pulled toward the group in proportion to how little their own choices pin
them down. The folds, held-out choices and scoring are those of every other analysis; the prior is estimated from
training trials only. Script 54b scores the pooled fits against the unpooled ones and builds the shared mixes
(script 46's method) from the pooled members.

## What the discovery sample showed

Horstmann, N = 162 (`output/pooled/log54b_horstmann.txt`):

- Pooling improved every model. NOVA: +454.6 (SE 140.2, z = 3.2), better for 98 of 162.
- Pooled NOVA against pooled ORL + α: +228.8 (SE 49.8, z = 4.6), better for 95. Unpooled, NOVA had lost to ORL + α
  (−55.8, SE 128.7), because a few people's fits collapsed.
- The shared mix of the eight pooled models against the shared mix of the eight unpooled: −126.8 (SE 24.8, z = −5.1).
  Pooled members made a worse mix.
- Five pooled models with NOVA (NOVA, ORL + split + α, VSE-split, EEF, with WSLS unpooled) against the eight unpooled:
  −42.1 (SE 30.5, z = −1.4).
- Pooled NOVA alone against the shared mix of the eight unpooled: −184.6 (SE 49.6, z = −3.7), better for 62. The best
  single model, even pooled, did not match the mix.

## Predictions

Replication sample, N = 440. Each test is a difference in held-out log-likelihood, summed over participants;
SE = standard deviation of the per-participant differences × √n; z = difference / SE. Each prediction is directional,
so each test is one-sided at α = .05 (the rule of the preregistered tests, script 48), with no correction for multiple
tests.

| | Prediction | Met if |
|---|---|---|
| Q1 | Pooling improves NOVA (pooled against unpooled NOVA) | z > 1.645 |
| Q2 | Pooled NOVA beats pooled ORL + α | z > 1.645 |
| Q3 | The shared mix of the eight pooled models is worse than the shared mix of the eight unpooled | z < −1.645 |
| Q4 | Five pooled models with NOVA do not beat the shared mix of the eight unpooled | z < 1.645 |
| Q5 | Pooled NOVA alone does worse than the shared mix of the eight unpooled | z < −1.645 |

Described, not tested: pooled NOVA against the pooled five, each model's gain from pooling, and the group parameters.

## What follows from the results

- Q1 and Q2 met: NOVA, fit with a group prior, is the best single model, and it is reported as such in paper 2.
- Q2 not met: NOVA's advantage when pooled was specific to the discovery sample.
- Q3 met: pooling and mixing do not add up. Better single models made a worse mix, and the shared mix stays built from
  unpooled fits.
- Q5 met: one model, even the best one fit with a group prior, does not predict as well as the mix, so the mix's
  advantage is not only noise in single fits. This is the question left open at the end of paper 1.
- Q4 predicts no improvement. If the pooled five do beat the unpooled eight (z ≥ 1.645), that is reported as a result,
  not hidden as a failed prediction.

## How to run

From the project folder, after this file and the script are committed (about 5 hours):

```
for m in orl variant orl_alpha variant_alpha vse vse_split eef nova; do for f in 01 02 03 04 05 06 07 08 09 10; do
  [ -f output/pooled/cv/replication/${m}_fold$f.rds ] || echo "$m $f"
done; done | caffeinate -i xargs -P 6 -n 2 sh -c 'Rscript scripts/54_pooled_cv.R replication $0 $1 > /dev/null'
Rscript scripts/54b_pooled_score.R replication > output/pooled/log54b_replication.txt 2>&1
Rscript scripts/54c_pooled_predictions.R replication > output/pooled/log54c_replication.txt 2>&1
```

Any change to the script or these predictions after the commit will be reported as a deviation.
