#!/bin/bash
# The Centaur benchmark from start to finish, on any machine with a suitable GPU: Kodiak (kodiak/centaur.sh submits
# this script) or a rented GPU (centaur/RENTING.md). From the project folder:
#   bash centaur/run_centaur.sh [model] [output] [limit]
#     model:  marcelbinz/Llama-3.1-Centaur-70B (the default; two 80 GB GPUs)
#             marcelbinz/Llama-3.1-Minitaur-8B (the small version; one GPU with 24 GB or more, for a pilot)
#     output: where the scores go (default centaur/centaur70b_scores.csv; centaur/*.csv is never committed)
#     limit:  score only the first few participants, for a test
# Needs: the packages in centaur/requirements.txt; a Hugging Face login on this machine (huggingface-cli login, which
# asks for your token in the terminal: never put the token in a file or a message); the terms of
# marcelbinz/Psych-101-test accepted on the Hugging Face website; and data-raw/igt_trials.csv (from
# scripts/53_export_trials.R on the Mac; data-raw is never pushed, so copy the file over by hand).
# Steps it skips when their output already exists: downloading Psych-101 (about 1 GB) and matching our participants.
# The model weights (about 141 GB for the 70B) are downloaded to HF_HOME; point HF_HOME at a disk with room.
set -euo pipefail
MODEL=${1:-marcelbinz/Llama-3.1-Centaur-70B}
OUT=${2:-centaur/centaur70b_scores.csv}
LIMIT=${3:-0}
export HF_HOME=${HF_HOME:-$PWD/hf_cache}
# A login made before HF_HOME was set is saved in the default place; use it if there is none in HF_HOME
if [ -z "${HF_TOKEN_PATH:-}" ] && [ ! -f "$HF_HOME/token" ] && [ -f "$HOME/.cache/huggingface/token" ]; then
  export HF_TOKEN_PATH=$HOME/.cache/huggingface/token
fi
[ -f data-raw/igt_trials.csv ] || { echo "data-raw/igt_trials.csv is missing: copy it from the Mac"; exit 1; }
if [ ! -f data-raw/psych101/train.jsonl ] || [ ! -f data-raw/psych101/test.jsonl ]; then
  python centaur/download_psych101.py
fi
[ -f centaur/psych101_igt_match.csv ] || python centaur/match_psych101.py
python -c "import torch; print('GPUs:', torch.cuda.device_count(), [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())])"
python centaur/score_centaur.py --model "$MODEL" --out "$OUT" --limit "$LIMIT"
echo "scores written to $OUT"
