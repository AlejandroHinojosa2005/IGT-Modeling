"""
download_psych101.py - Saves Psych-101's training and held-out sets to data-raw/psych101/ as .jsonl files.

The held-out set (marcelbinz/Psych-101-test) is gated: accept its terms on Hugging Face, then log in on this
computer once with  huggingface-cli login  (it asks for your token in the terminal; never paste it anywhere else).
Needs about 1 GB of disk space. Run from the project folder:  python centaur/download_psych101.py
"""
import json
import os

from datasets import load_dataset

os.makedirs("data-raw/psych101", exist_ok=True)
for name, out in (("marcelbinz/Psych-101", "train"), ("marcelbinz/Psych-101-test", "test")):
    ds = load_dataset(name)
    split = list(ds.keys())[0]                                 # each repository holds one split
    with open(f"data-raw/psych101/{out}.jsonl", "w") as f:
        for row in ds[split]:
            f.write(json.dumps({k: row[k] for k in ("text", "experiment", "participant") if k in row}) + "\n")
    print(name, split, len(ds[split]), "rows ->", f"data-raw/psych101/{out}.jsonl")
