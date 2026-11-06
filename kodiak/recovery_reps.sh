#!/bin/bash
#PBS -N igt_recov
#PBS -J 1-10
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

Rscript scripts/26_recovery_reps.R $PBS_ARRAY_INDEX $QUICK > logs/recovery_rep${PBS_ARRAY_INDEX}${QUICK:+_quick}.log 2>&1
