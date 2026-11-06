#!/bin/bash
#PBS -N orl_mix
#PBS -J 0-1
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

FOLDS=(0 8)
FOLD=${FOLDS[$PBS_ARRAY_INDEX]}

Rscript scripts/24_orl_mixing.R $FOLD $QUICK > logs/orl_mixing_fold${FOLD}${QUICK:+_quick}.log 2>&1
