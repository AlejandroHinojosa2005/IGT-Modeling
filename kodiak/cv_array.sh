   #!/bin/bash
   #PBS -N igt_cv
   #PBS -J 0-21
   #PBS -l select=1:ncpus=4:mem=16gb
   #PBS -l walltime=12:00:00
   #PBS -j oe

   cd $PBS_O_WORKDIR
   module load R/4.5.2
   mkdir -p logs

   MODELS=(orl wl)
   MODEL=${MODELS[$((PBS_ARRAY_INDEX / 11))]}
   FOLD=$((PBS_ARRAY_INDEX % 11))

   Rscript scripts/23_cv_fold.R $MODEL $FOLD $QUICK > logs/cv_${MODEL}_fold${FOLD}${QUICK:+_quick}.log 2>&1






