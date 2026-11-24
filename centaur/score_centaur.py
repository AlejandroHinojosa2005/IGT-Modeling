"""
score_centaur.py - Centaur's probability for every choice of the participants it never trained on.

Centaur predicts each choice from the transcript of the trials before it, so every score here is a forward
prediction (choice t from trials 1 to t - 1). Two numbers per choice:
  lp_token  log-probability of the chosen deck's letter token (Centaur's own metric, Binz et al., 2025)
  lp_4      the same, renormalized over the four deck letters: a probability over the four decks, comparable with
            our models' predictions
Participants scored (from centaur/psych101_igt_match.csv, made by match_psych101.py):
  test    held out from Centaur's training: scored on their own Psych-101 transcript
  absent  never in Psych-101: a transcript in the same format is written from our data (the instructions are
          copied from a Psych-101 IGT transcript; deck letters drawn at random, seed 2026)
Participants Centaur trained on are skipped. Only the first 100 choices are scored.

Run from the project folder, after match_psych101.py:
  python centaur/score_centaur.py --model marcelbinz/Llama-3.1-Minitaur-8B --out centaur/minitaur8b_scores.csv
  add --limit 3 for a quick test on three participants
Memory: the 8B model needs about 16 GB in 16-bit (a GPU, or a Mac with more memory than that); the 70B needs about
140 GB (two 80 GB GPUs).
"""
import argparse
import csv
import json
import random
import re
from collections import defaultdict

import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

CHOICE = re.compile(r"<<\s*([A-Za-z])\s*>>")
LABELS = re.compile(r"labeled\s+([A-Za-z]),\s*([A-Za-z]),\s*([A-Za-z]),?\s*and\s+([A-Za-z])")


def load_ours(path="data-raw/igt_trials.csv"):
    people = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            people[(r["sample"], int(r["subj"]))].append(
                (int(r["trial"]), int(r["choice"]), float(r["win"]), abs(float(r["loss"]))))
    return {k: [x[1:] for x in sorted(v)] for k, v in people.items()}


def write_transcript(template, seq, rng):
    """A Psych-101-style transcript of our data: Psych-101's instructions and trial-line wording, new deck letters."""
    head = template[:template.index("You press")]
    assert LABELS.search(head), "no deck labels in the instructions"
    first = template[template.index("You press"):].split("\n")[0]          # e.g. You press <<H>>. You win 100.0$ and lose 200.0$.
    nums = re.findall(r"-?\d+(?:\.(\d+))?(?=\s*\$)", first)
    digits = len(nums[0]) if nums and nums[0] else 0                    # decimals used for amounts
    line = CHOICE.sub("<<{L}>>", first, count=1)
    line = re.sub(r"-?\d+(?:\.\d+)?(?=\s*\$)", "{N}", line)
    assert line.count("{N}") == 2, "unexpected trial-line wording: " + first
    a, b = line.split("{N}", 2)[0], line.split("{N}", 2)[1]
    rest = line.split("{N}", 2)[2]
    new = rng.sample([c for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"], 4)
    head = LABELS.sub(f"labeled {new[0]}, {new[1]}, {new[2]}, and {new[3]}", head, count=1)
    lines = [a.replace("{L}", new[c - 1]) + f"{w:.{digits}f}" + b + f"{l:.{digits}f}" + rest for c, w, l in seq[:100]]
    return head + "\n".join(lines) + "\n", new


def score(model, tok, text, letters, device):
    """Per choice: (letter, lp_token, lp_4), for the first 100 choices."""
    marks = list(CHOICE.finditer(text))[:100]
    text = text[:marks[-1].end()]                              # later trials cannot change earlier predictions
    enc = tok(text, return_offsets_mapping=True, return_tensors="pt")
    ids = enc["input_ids"].to(device)
    starts = [s for s, _ in enc["offset_mapping"][0].tolist()]
    letter_ids = [tok.encode(L, add_special_tokens=False) for L in letters]
    assert all(len(x) == 1 for x in letter_ids), "a deck letter is not a single token"
    letter_ids = [x[0] for x in letter_ids]
    with torch.no_grad():
        logits = model(ids).logits[0].float()
    out = []
    for m in marks:
        pos = m.start(1)
        t = starts.index(pos) if pos in starts else None
        assert t is not None and tok.decode(ids[0, t]) == m.group(1), \
            "the chosen letter is not its own token here; check the tokenization before scoring"
        lp = torch.log_softmax(logits[t - 1], -1)
        lp4 = torch.log_softmax(logits[t - 1, letter_ids], -1)
        k = letters.index(m.group(1))
        out.append((m.group(1), lp[ids[0, t]].item(), lp4[k].item()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="marcelbinz/Llama-3.1-Minitaur-8B")
    ap.add_argument("--out", default="centaur/minitaur8b_scores.csv")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
    dtype = torch.bfloat16 if device == "cuda" else torch.float16 if device == "mps" else torch.float32
    tok = AutoTokenizer.from_pretrained(a.model)
    model = AutoModelForCausalLM.from_pretrained(a.model, dtype=dtype,
                                                 device_map="auto" if device == "cuda" else None)
    if device != "cuda":
        model.to(device)
    model.eval()

    match = list(csv.DictReader(open("centaur/psych101_igt_match.csv")))
    test_rows = {}
    template = None
    with open("data-raw/psych101/test.jsonl") as f:
        for line in f:
            row = json.loads(line)
            if "steingroever" in str(row.get("experiment", "")).lower():
                test_rows[(row["experiment"], str(row["participant"]))] = row["text"]
                template = template or row["text"]
    ours = load_ours()
    rng = random.Random(2026)
    todo = [m for m in match if m["split"] in ("test", "absent")]
    if a.limit:
        todo = todo[:a.limit]

    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "subj", "split", "trial", "deck", "lp_token", "lp_4"])
        for i, m in enumerate(todo):
            key = (m["sample"], int(m["subj"]))
            if m["split"] == "test":
                text = test_rows[(m["experiment"], m["participant"])]
                letters = list(LABELS.search(text).groups())
            else:
                text, letters = write_transcript(template, ours[key], rng)
            for trial, ((letter, lp_tok, lp4), (deck, _, _)) in enumerate(zip(score(model, tok, text, letters, device), ours[key]), 1):
                w.writerow([*key, m["split"], trial, deck, f"{lp_tok:.6f}", f"{lp4:.6f}"])
            print(f"{i + 1}/{len(todo)}", key, m["split"], flush=True)


if __name__ == "__main__":
    main()
