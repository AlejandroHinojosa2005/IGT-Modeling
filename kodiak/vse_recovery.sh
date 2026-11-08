#!/bin/bash
#PBS -N igt_vse_rec
#PBS -J 0-1
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

GENS=(vse vse_split)
GEN=${GENS[$PBS_ARRAY_INDEX]}

Rscript scripts/33_vse_recovery.R $GEN $QUICK > logs/vse_recovery_${GEN}${QUICK:+_quick}.log 2>&1
