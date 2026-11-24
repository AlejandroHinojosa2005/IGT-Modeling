#!/bin/bash
#PBS -N igt_centaur
#PBS -l select=1:ncpus=8:ngpus=2:mem=200gb
#PBS -l walltime=08:00:00
#PBS -j oe

# The Centaur benchmark on Kodiak's GPUs (centaur/run_centaur.sh). Centaur 70B in 16-bit needs about 141 GB of GPU
# memory: two 80 GB GPUs (A100 or H100). The scoring itself takes minutes; downloading the weights takes longest.
# To check with Baylor IT before the first run (the lines marked CHECK):
#   - the GPU queue (add a line  #PBS -q <queue>  if GPUs have their own) and how to ask for 80 GB GPUs;
#   - the Python and CUDA module names (module avail python cuda);
#   - where scratch space is (the weights need about 150 GB) and whether compute nodes can reach the internet.
# Once, on the login node, from the project folder:
#   module load python/3.11                       # CHECK the name
#   python -m venv ~/centaur-env && source ~/centaur-env/bin/activate
#   pip install -r centaur/requirements.txt
#   huggingface-cli login                         # type your token here; it is never written into the project
#   copy data-raw/igt_trials.csv from your Mac (scp)
#   python centaur/download_psych101.py && python centaur/match_psych101.py
#   If compute nodes cannot reach the internet, download the weights here first:
#   HF_HOME=/scratch/$USER/hf_cache huggingface-cli download marcelbinz/Llama-3.1-Centaur-70B
# Then submit:  qsub kodiak/centaur.sh
# Pilot with the small model (one GPU is enough; change ngpus=2 to ngpus=1 above):
#   qsub -v MODEL=marcelbinz/Llama-3.1-Minitaur-8B,OUT=centaur/minitaur8b_scores.csv,LIMIT=3 kodiak/centaur.sh
# Copy the scores back to the Mac (centaur/*.csv); nothing is committed or pushed from Kodiak.
cd $PBS_O_WORKDIR
module load python/3.11 cuda                      # CHECK the names
source ~/centaur-env/bin/activate
mkdir -p logs
export HF_HOME=${HF_HOME:-/scratch/$USER/hf_cache}  # CHECK the scratch path

bash centaur/run_centaur.sh "${MODEL:-marcelbinz/Llama-3.1-Centaur-70B}" "${OUT:-centaur/centaur70b_scores.csv}" "${LIMIT:-0}" \
  > logs/centaur_$(basename "${MODEL:-Centaur-70B}").log 2>&1
