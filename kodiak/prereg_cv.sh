#!/bin/bash
#PBS -N igt_prereg_cv
#PBS -J 0-279
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -j oe

# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), step 2: held-out predictions (script 45).
# 7 models x 10 folds x 4 parts of the participants = 280 jobs. Choose the data set when submitting:
#   qsub -v DATA=replication kodiak/prereg_cv.sh
#   qsub -v DATA=simulated kodiak/prereg_cv.sh      (after script 44 has made the simulated data)
# For a quick test (2 participants per job), add QUICK=quick: qsub -v DATA=replication,QUICK=quick kodiak/prereg_cv.sh
cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

PARTS=4
MODELS=(orl variant orl_alpha variant_alpha vse vse_split eef)
MODEL=${MODELS[$((PBS_ARRAY_INDEX / (10 * PARTS)))]}
FOLD=$(( (PBS_ARRAY_INDEX / PARTS) % 10 + 1 ))
PART=$(( PBS_ARRAY_INDEX % PARTS + 1 ))

Rscript scripts/45_prereg_cv.R $DATA $MODEL $FOLD $PART $PARTS $QUICK \
  > logs/prereg_cv_${DATA}_${MODEL}_fold${FOLD}_part${PART}${QUICK:+_quick}.log 2>&1
