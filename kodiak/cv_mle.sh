#!/bin/bash
#PBS -N igt_cvmle
#PBS -J 0-69
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

MODELS=(orl variant orl_alpha variant_alpha vse vse_split eef)
MODEL=${MODELS[$((PBS_ARRAY_INDEX / 10))]}
FOLD=$((PBS_ARRAY_INDEX % 10 + 1))

Rscript scripts/35_cv_mle.R $MODEL $FOLD $QUICK > logs/cv_mle_${MODEL}_fold${FOLD}${QUICK:+_quick}.log 2>&1
