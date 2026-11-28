#!/bin/bash
# Submits the hierarchical fits of one data set on Kodiak as two chained job arrays: every fit (kodiak/stan_cv.sh), then,
# once that array has ended, the reruns of the fits that did not pass the convergence rule. Each job stays under the
# 12-hour limit, and Kodiak runs as many at once as the queue allows.
# Usage (from the project folder, after git pull and, once, bash kodiak/stan_setup.sh):
#   bash kodiak/stan_submit.sh replication
#   bash kodiak/stan_submit.sh horstmann
# If the Mac has already done some fits, copy them here first so they are skipped, from the Mac's project folder:
#   scp -r output/stan <you>@<kodiak address>:<path to>/IGT-Modeling/output/
# Check on the jobs with: qstat -u $USER   (and qstat -t for each array's jobs)
set -e
DATA=$1
[ -n "$DATA" ] || { echo "usage: bash kodiak/stan_submit.sh <horstmann|replication>"; exit 1; }
FIRST=$(qsub -v DATA=$DATA kodiak/stan_cv.sh)
echo "pass 1 (every fit): $FIRST"
SECOND=$(qsub -W depend=afterany:$FIRST -v DATA=$DATA,LONG=long kodiak/stan_cv.sh)
echo "pass 2 (reruns, starts when pass 1 has ended): $SECOND"
