# Forward predictions: predictions

Written November 24, 2026, before any forward predictions were run on either sample. The scripts that make and score
them (`scripts/57_forward_scores.R`, `scripts/57b_forward_compare.R`) were committed the day before. This is an
exploratory follow-up to the preregistered replication (<https://doi.org/10.17605/OSF.IO/U86SY>), not part of that
registration.

## What forward prediction is

Every other score in this project fits each person to 90% of their trials, later ones included, and scores the rest.
Here each choice is predicted only from the trials before it, the way a person meets the task and the way Centaur is
scored. Each model is refit after every trial under a group prior estimated from the other sample (script 54's pooled
fits); the shared mixes use stacking weights fit in the other sample. Nothing about the person being scored, other
than their earlier choices, informs the prediction.

## Predictions

Each sample separately: Horstmann (N = 162, prior and weights from the replication sample) and the replication sample
(N = 440, prior and weights from Horstmann). Each test is a difference in log-probability summed over trials 2 to 100
and over participants; SE = standard deviation of the per-participant differences × √n; z = difference / SE.
One-sided at α = .05. Script 57b reports every model and mix against the mix of the eight.

| | Prediction | Met if |
|---|---|---|
| F1 | The shared mix of the eight beats ORL + α, the best single model in the preregistered analysis | ORL + α against the mix: z < −1.645, in each sample |

Why it could fail: here every model is fit under a group prior, and pooled members made a worse mix in the pooled
replication (`predictions/pooled_replication.md`, Q3: z = −9.7).

Described, not tested: every single model's forward score and which is best, the mix of five with NOVA, and NOVA
alone against the mix.

## What follows from the results

- F1 met in both samples: the advantage of mixing does not depend on fitting with later trials; it holds when each
  choice is predicted from the past alone.
- F1 not met: the mix's edge comes partly from seeing a person's later choices, and paper 2 says so.

## How to run

From the project folder (about 5 seconds per person and model on one core, NOVA about 10):

    for d in horstmann replication; do for m in orl variant orl_alpha variant_alpha vse vse_split eef nova wsls; do
      for p in 1 2 3 4; do echo "$d $m $p"; done; done; done |
      caffeinate -i xargs -P 4 -n 3 sh -c 'Rscript scripts/57_forward_scores.R "$0" "$1" part=$2 parts=4 > "logs/forward_$0_$1_$2.log" 2>&1'
    Rscript scripts/57b_forward_compare.R horstmann > output/forward/log57b_horstmann.txt 2>&1
    Rscript scripts/57b_forward_compare.R replication > output/forward/log57b_replication.txt 2>&1
