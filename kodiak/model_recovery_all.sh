#!/bin/bash
#PBS -N igt_recovery_all
#PBS -J 0-7
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=24:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

GENS=(orl variant orl_alpha variant_alpha vse vse_split eef wsls)
GEN=${GENS[$PBS_ARRAY_INDEX]}

Rscript scripts/39_model_recovery_all.R $GEN $QUICK > logs/model_recovery_${GEN}${QUICK:+_quick}.log 2>&1
