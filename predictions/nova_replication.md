# NOVA on the replication sample: predictions

Written November 23, 2026, before the analysis was run, and committed together with the script that scores it
(`scripts/50d_nova_replication.R`). This is an exploratory follow-up to the preregistered replication
(<https://doi.org/10.17605/OSF.IO/U86SY>), not part of that registration.

## What NOVA is

NOVA is the hybrid model of script 50 (version `hybrid`): ORL + α, with a starting value for every deck (Q0) and
decay of the unchosen decks' values toward 0 (D). Script 50b fits it with the same ten folds, held-out choices and
fitting as every other model in the replication: 20 random starts (seed 2050 + fold), plus ORL + α's own fit as one
more start.

## What the discovery sample showed

Horstmann, N = 162 (`output/hybrid/log50c_horstmann.txt`):

- NOVA alone against ORL + α: −55.8 (SE 128.7). No better. It helped the median person slightly (+0.6), but a few
  people's predictions collapsed.
- Adding NOVA to the shared mix of the eight models: +82.1 (SE 20.6, z = 4.0). NOVA got the largest weight (.355).
- Five models with NOVA (NOVA, ORL + split + α, VSE-split, EEF, win-stay/lose-shift) against all eight: +75.7 (SE 22.7,
  z = 3.3). Against the same five with ORL + α in place of NOVA: +84.6 (SE 22.7, z = 3.7).

The five were chosen after seeing script 52's results, so these are claims to test, not findings.

## Predictions

Replication sample, N = 440. Each test is a difference in held-out log-likelihood, summed over participants;
SE = standard deviation of the per-participant differences × √n; z = difference / SE. Each prediction is directional, so
each test is one-sided at α = .05 (z > 1.645), the rule of the preregistered tests (script 48), with no correction for
multiple tests. Mixes are fit as script 46's
shared mix: one set of weights for everyone, fit on the other nine folds and scored on the held-out fold.

| | Prediction | Met if |
|---|---|---|
| P1 | NOVA alone does not beat ORL + α | z < 1.645 |
| P2 | Adding NOVA to the eight improves the shared mix | z > 1.645 |
| P3 | Five models with NOVA beat the shared mix of all eight | z > 1.645 |
| P4 | Five models with NOVA beat the same five with ORL + α | z > 1.645 |

Described, not tested: NOVA's weight in the mix of nine, the median per-person difference from ORL + α, the five
largest losses, and the fitted Q0 and D.

## What follows from the results

- P2, P3 and P4 met: the five-model Seldon with NOVA is carried into later work.
- P2 met, P3 or P4 not: NOVA joins the mix, but the eight (plus NOVA) stay.
- P2 not met: NOVA's gain in the mix was specific to the discovery sample, and the eight stay.
- P1 predicts no improvement. If NOVA alone does beat ORL + α (z ≥ 1.645), that is reported as a result in NOVA's
  favor, not hidden as a failed prediction.

## How to run

From the project folder, after this file and the script are committed:

```
for f in 01 02 03 04 05 06 07 08 09 10; do for p in 1 2 3 4; do
  [ -f output/hybrid/cv/replication/hybrid_fold${f}_part$p.rds ] || echo "hybrid $f $p"
done; done | caffeinate -i xargs -P 6 -n 3 sh -c 'Rscript scripts/50b_hybrid_cv.R replication $0 $1 $2 4 > /dev/null'
Rscript scripts/50d_nova_replication.R > output/hybrid/log50d_replication.txt 2>&1
```

Any change to the script or these predictions after the commit will be reported as a deviation.
