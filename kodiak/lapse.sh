#!/bin/bash
#PBS -N igt_lapse
#PBS -J 0-1
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

MODELS=(orl wl)
MODEL=${MODELS[$PBS_ARRAY_INDEX]}

Rscript scripts/25_lapse.R $MODEL $QUICK > logs/lapse_${MODEL}${QUICK:+_quick}.log 2>&1
