#!/bin/bash
# Runs every hierarchical fit of scripts/55_stan_fit.R for one data set (8 models x (all trials + 10 folds) = 88 fits,
# plus the reruns of the convergence rule) on one machine, in the foreground: a rented many-core computer or an
# interactive cluster node. On the Mac use scripts/55_queue.sh instead, which runs in the background and can be paused.
# (On Kodiak's batch system use kodiak/stan_submit.sh.) Each fit runs 4 chains on 4 cores, so set JOBS to about
# cores / 4. Saved fits are skipped, so it can be stopped and started again.
# Usage (from the project folder):  bash scripts/55_run_all.sh <data> [JOBS]
#   e.g. bash scripts/55_run_all.sh replication 16   (a rented 64-core computer: about a day)
# Then: Rscript scripts/55c_stan_cv_score.R <data> and Rscript scripts/55d_stan_predictions.R <data>
set -e
DATA=$1
JOBS=${2:-2}
[ -n "$DATA" ] || { echo "usage: bash scripts/55_run_all.sh <horstmann|replication> [JOBS]"; exit 1; }
bash scripts/55_queue.sh check
bash scripts/55_queue.sh _run "$JOBS" "$DATA"
tail -n 3 logs/stan_queue.log
