#!/bin/bash
# The analyses without the two repeated records (scripts/43_prereg_data.R, duplicate_rows; data set "replication438"),
# one step after another. Every model is fit to each person on their own, so those fits are reused; the steps that
# combine participants are rerun: the shared mix and switching (46, 47), the registered tests (58), NOVA's tests (50d),
# the pooled fits, whose group prior is estimated from everyone (54, 54b, 54c), and the forward predictions for
# Horstmann, which use the replication sample's group prior (57, 57b). The replication sample's forward predictions use
# the Horstmann prior, so 57b only drops the two records from them.
# Resumable: pooled and forward runs whose output already exists are skipped, so after an interruption run it again.
# From the project folder:   caffeinate -i bash scripts/58_duplicates_chain.sh
# About 6 to 7 hours on an Apple silicon MacBook Air: steps 1-4 about 20 minutes, step 5 (the pooled refit) about 5 to
# 6 hours, steps 6-8 about 30 minutes. Progress: logs/chain_duplicates.log; results: the log files named below.
set -u
cd "$(dirname "$0")/.."
mkdir -p logs output/prereg output/hybrid output/pooled output/forward
step() { echo "$(date '+%Y-%m-%d %H:%M') $*" | tee -a logs/chain_duplicates.log; }
run() { local out=$1; shift; "$@" > "$out" 2>&1 || { step "FAILED: $* (see $out)"; exit 1; }; }

step "1. shared mix without the repeated records (46)"
run output/prereg/log46_replication438.txt Rscript scripts/46_prereg_stacking.R replication438
step "2. switching (47)"
run output/prereg/log47_replication438.txt Rscript scripts/47_prereg_switching.R replication438
step "3. the registered tests on 438 next to 440 (58)"
run output/prereg/log58_duplicates.txt Rscript scripts/58_duplicates.R
step "4. NOVA's tests (50d)"
run output/hybrid/log50d_replication438.txt Rscript scripts/50d_nova_replication.R replication438

step "5. pooled refit: 8 models x 10 folds, 6 at a time (54)"
for m in orl variant orl_alpha variant_alpha vse vse_split eef nova; do for f in 01 02 03 04 05 06 07 08 09 10; do
  [ -f "output/pooled/cv/replication438/${m}_fold$f.rds" ] || echo "$m $f"
done; done | xargs -P 6 -n 2 sh -c 'Rscript scripts/54_pooled_cv.R replication438 "$0" "$1" > "logs/pooled_replication438_$0_$1.log" 2>&1'
n=$(ls output/pooled/cv/replication438/*_fold??.rds 2>/dev/null | wc -l | tr -d ' ')
[ "$n" -eq 80 ] || { step "FAILED: $n of 80 pooled fits (see logs/pooled_replication438_*.log)"; exit 1; }
step "6. pooled scores (54b, 54c)"
run output/pooled/log54b_replication438.txt Rscript scripts/54b_pooled_score.R replication438
run output/pooled/log54c_replication438.txt Rscript scripts/54c_pooled_predictions.R replication438

step "7. forward predictions for Horstmann with the replication438 prior: 9 models x 4 parts, 4 at a time (57)"
for m in orl variant orl_alpha variant_alpha vse vse_split eef nova wsls; do for p in 1 2 3 4; do
  [ -f "output/forward/horstmann_prior-replication438/${m}_part$p.rds" ] || echo "$m $p"
done; done | xargs -P 4 -n 2 sh -c 'Rscript scripts/57_forward_scores.R horstmann "$0" part=$1 parts=4 prior=replication438 > "logs/forward_horstmann_prior438_$0_$1.log" 2>&1'
n=$(ls output/forward/horstmann_prior-replication438/*_part?.rds 2>/dev/null | wc -l | tr -d ' ')
[ "$n" -eq 36 ] || { step "FAILED: $n of 36 forward runs (see logs/forward_horstmann_prior438_*.log)"; exit 1; }
step "8. forward comparisons (57b)"
run output/forward/log57b_horstmann_prior-replication438.txt Rscript scripts/57b_forward_compare.R horstmann prior=replication438
run output/forward/log57b_replication438.txt Rscript scripts/57b_forward_compare.R replication438
step "done"
