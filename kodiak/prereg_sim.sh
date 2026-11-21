#!/bin/bash
#PBS -N igt_prereg_sim
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=1:00:00
#PBS -j oe

# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), H6: makes the no-switching data (script 44).
# Takes under a minute. Then submit the cross-validation of the simulated data (kodiak/prereg_cv.sh, DATA=simulated).
cd $PBS_O_WORKDIR
module load R/4.5.2
mkdir -p logs

Rscript scripts/44_null_simulation.R > logs/prereg_sim.log 2>&1
