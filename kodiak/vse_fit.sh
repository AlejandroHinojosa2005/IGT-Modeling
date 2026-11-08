#!/bin/bash
#PBS -N igt_vse
#PBS -J 0-2
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

MODELS=(vse eef wsls)
MODEL=${MODELS[$PBS_ARRAY_INDEX]}

Rscript scripts/32_fit_vse.R $MODEL $QUICK > logs/vse_fit_${MODEL}${QUICK:+_quick}.log 2>&1
