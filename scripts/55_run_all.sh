#!/bin/bash
# Runs every hierarchical Stan fit of scripts/55_stan_fit.R for one data set, 8 models x (all trials + 10 folds) = 88
# fits, on a single machine: your Mac, a rented many-core computer, or an interactive cluster node. (On Kodiak's batch
# system use kodiak/stan_cv.sh instead.) Each fit runs 4 chains on 4 cores, so set JOBS to about cores / 4. Fits that
# are already saved are skipped, so it can be stopped and started again.
# Usage (from the project folder):  bash scripts/55_run_all.sh <data> [JOBS]
#   e.g. bash scripts/55_run_all.sh horstmann 2      (an 8-core laptop: about 2 days for Horstmann)
#        bash scripts/55_run_all.sh replication 16   (a rented 64-core computer: about a day)
# Then: Rscript scripts/55c_stan_cv_score.R <data>
set -e
export DATA=$1
JOBS=${2:-2}
[ -n "$DATA" ] || { echo "usage: bash scripts/55_run_all.sh <horstmann|replication> [JOBS]"; exit 1; }
mkdir -p logs
Rscript -e 'library(cmdstanr); for (f in Sys.glob("stan/hier/igt_*.stan")) invisible(cmdstan_model(f))'   # compile once
for m in orl orl_split orl_alpha orl_split_alpha nova vse vse_split eef; do
  [ -f output/stan/$DATA/${m}_full.rds ] || echo "$m full"
  for f in 1 2 3 4 5 6 7 8 9 10; do
    [ -f output/stan/cv/$DATA/${m}_fold$(printf %02d $f).rds ] || echo "$m $f"
  done
done | xargs -P "$JOBS" -n 2 sh -c 'Rscript scripts/55_stan_fit.R "$DATA" "$0" "$1" > "logs/stan_${DATA}_$0_$1.log" 2>&1'
echo "all fits done"
