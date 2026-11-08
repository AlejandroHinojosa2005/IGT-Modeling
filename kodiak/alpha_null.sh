#!/bin/bash
#PBS -N igt_alpha_null
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=12:00:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

Rscript scripts/29b_alpha_null.R $QUICK > logs/alpha_null${QUICK:+_quick}.log 2>&1