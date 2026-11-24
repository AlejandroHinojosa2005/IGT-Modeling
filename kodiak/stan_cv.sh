#!/bin/bash
#PBS -N igt_stan
#PBS -J 0-87
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=24:00:00
#PBS -j oe

# Hierarchical Bayesian fits (scripts/55_stan_fit.R): 8 models x (all trials + 10 cross-validation folds) = 88 jobs per
# data set, each running 4 chains on 4 cores. Rough times: about 1 h per fit for Horstmann (N = 162), about 3 h for the
# replication sample (N = 440); NOVA and the VSE models may take longer.
# First, once: bash kodiak/stan_setup.sh. Then choose the data set when submitting:
#   qsub -v DATA=horstmann kodiak/stan_cv.sh
#   qsub -v DATA=replication kodiak/stan_cv.sh
# Quick test (the first 20 participants, short chains): qsub -v DATA=horstmann,QUICK=quick kodiak/stan_cv.sh
# Afterwards (any machine, with the outputs copied back): Rscript scripts/55c_stan_cv_score.R <data>
# Results come back to the Mac by copying output/stan/; nothing is committed or pushed from Kodiak.
cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

MODELS=(orl orl_split orl_alpha orl_split_alpha nova vse vse_split eef)
MODEL=${MODELS[$((PBS_ARRAY_INDEX / 11))]}
F=$((PBS_ARRAY_INDEX % 11))
if [ "$F" -eq 0 ]; then FOLD=full; else FOLD=$F; fi

Rscript scripts/55_stan_fit.R $DATA $MODEL $FOLD $QUICK > logs/stan_${DATA}_${MODEL}_${FOLD}${QUICK:+_quick}.log 2>&1
