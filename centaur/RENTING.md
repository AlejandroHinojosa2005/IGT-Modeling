# Running the Centaur benchmark on a rented GPU

The same script that runs on Kodiak (`centaur/run_centaur.sh`) runs on a rented GPU. Nothing in the code changes.

## What to rent

- **Centaur 70B:** two 80 GB GPUs (2 × A100 80GB or 2 × H100 80GB), with at least 250 GB of disk (the weights
  are about 141 GB). On RunPod, a "PyTorch" template on a 2 × A100 80GB pod.
- **Minitaur 8B (pilot):** one GPU with 24 GB or more, 50 GB of disk.

Time: downloading and loading the 70B takes the longest (roughly 15 to 45 minutes); scoring the participants Centaur
did not train on takes minutes, because each participant's 100 choices are scored in one pass. At roughly $1.20 to
$2.20 per GPU per hour, the whole run should cost well under $20. Billing runs until the pod is stopped and deleted.

## Steps

1. On Hugging Face, accept the terms of `marcelbinz/Psych-101-test` (once, on the website).
2. Start the pod and open its terminal. Then:
   ```
   git clone https://github.com/AlejandroHinojosa2005/IGT-Modeling.git
   cd IGT-Modeling
   pip install -r centaur/requirements.txt
   huggingface-cli login            # paste your token here, in the pod's terminal only
   mkdir -p data-raw
   ```
3. Upload `data-raw/igt_trials.csv` from the Mac into the pod's `IGT-Modeling/data-raw/` (the pod's file upload,
   or `scp`). It is not on GitHub: `data-raw/` is never pushed.
4. Put the weights on the big disk and run (a three-person test first, then everyone):
   ```
   export HF_HOME=/workspace/hf_cache
   bash centaur/run_centaur.sh marcelbinz/Llama-3.1-Centaur-70B centaur/centaur70b_test.csv 3
   bash centaur/run_centaur.sh
   ```
5. Download `centaur/centaur70b_scores.csv` (and `centaur/psych101_igt_match.csv`) to the Mac's `centaur/` folder.
6. Stop and delete the pod, so its disk (with the Psych-101 test set and the weights) is gone and billing ends.

## Then, on the Mac

Forward predictions of our models for the same participants, and the comparison (see `REPLICATION.md`):
```
for d in horstmann replication; do for m in orl variant orl_alpha variant_alpha vse vse_split eef nova wsls; do
  echo "$d $m"
done; done | xargs -P 6 -n 2 sh -c 'Rscript scripts/57_forward_scores.R $0 $1 who=centaur/psych101_igt_match.csv > /dev/null'
Rscript scripts/57b_forward_compare.R horstmann centaur/centaur70b_scores.csv
Rscript scripts/57b_forward_compare.R replication centaur/centaur70b_scores.csv
```

## The hierarchical cross-validation on a rented computer

It needs no GPU. On a rented computer with many CPU cores (64 cores runs 16 fits at once), with R and the packages
from `renv.lock` installed:
```
bash scripts/55_run_all.sh horstmann 16
bash scripts/55_run_all.sh replication 16
Rscript scripts/55c_stan_cv_score.R horstmann
```
