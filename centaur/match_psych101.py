"""
match_psych101.py - Which of our participants did Centaur see in training?

Centaur (Binz et al., 2025) was fine-tuned on Psych-101, whose Iowa gambling task data are the Steingroever et al.
(2015) pool (reference 116 of the Centaur paper), split by participant: about 90% for training, 10% held out.
This script finds each of our participants (scripts 43 and 53: Horstmann, N = 162; replication, N = 440) among the
Psych-101 IGT transcripts, by the sequence of outcomes on their first 100 trials, and records the split.

Inputs (all local; nothing here is pushed):
  data-raw/igt_trials.csv            from: Rscript scripts/53_export_trials.R
  data-raw/psych101/train.jsonl      the Psych-101 training set (Hugging Face: marcelbinz/Psych-101)
  data-raw/psych101/test.jsonl       the held-out set (Hugging Face: marcelbinz/Psych-101-test; accept its terms first)
Each line of the .jsonl files is one participant: {"text": ..., "experiment": ..., "participant": ...}.
Download them on the Mac with download_psych101.py (it saves the two files above).

Output: centaur/psych101_igt_match.csv, one row per participant of ours:
  sample, subj, study, split (train / test / absent), experiment, participant (Psych-101's ids)

A transcript line looks like:  You press <<H>>. You win 100.0$ and lose 200.0$.
Deck letters differ between participants; they are mapped to our decks by the choices, and the mapping must be
one-to-one for a match to count.

Run from the project folder:  python centaur/match_psych101.py
"""
import csv
import json
import re
from collections import defaultdict

TRIAL = re.compile(r"<<\s*([A-Za-z])\s*>>[^\n]*?win\s+(-?[\d.]+)\s*\$[^\n]*?lose\s+(-?[\d.]+)\s*\$")


def parse_transcript(text):
    """Letters chosen, wins and losses (as positive amounts), in order."""
    return [(m.group(1), float(m.group(2)), abs(float(m.group(3)))) for m in TRIAL.finditer(text)]


def load_ours(path="data-raw/igt_trials.csv"):
    people = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            people[(r["sample"], int(r["subj"]), r["study"])].append(
                (int(r["trial"]), int(r["choice"]), float(r["win"]), abs(float(r["loss"]))))
    out = {}
    for key, rows in people.items():
        rows.sort()
        out[key] = [(c, w, l) for _, c, w, l in rows]
    return out


def outcome_key(seq):
    return tuple((round(w, 2), round(l, 2)) for _, w, l in seq[:100])


def same_choices(ours, theirs):
    """True if the letters map one-to-one onto our deck numbers over the first 100 trials."""
    to_deck, to_letter = {}, {}
    for (deck, _, _), (letter, _, _) in zip(ours[:100], theirs[:100]):
        if to_deck.setdefault(letter, deck) != deck or to_letter.setdefault(deck, letter) != letter:
            return False
    return True


def main():
    ours = load_ours()
    index = defaultdict(list)
    for key, seq in ours.items():
        index[outcome_key(seq)].append(key)

    found = {}
    counts = defaultdict(int)
    for split in ("train", "test"):
        with open(f"data-raw/psych101/{split}.jsonl") as f:
            for line in f:
                row = json.loads(line)
                if "steingroever" not in str(row.get("experiment", "")).lower():
                    continue
                counts[split] += 1
                seq = parse_transcript(row["text"])
                if len(seq) < 100:
                    continue                                   # the 95-trial study is not in our samples
                for key in index.get(outcome_key(seq), []):
                    if same_choices(ours[key], seq):
                        if key in found and found[key][0] != split:
                            raise RuntimeError(f"{key} matched in both splits")
                        found[key] = (split, row.get("experiment"), row.get("participant"))

    with open("centaur/psych101_igt_match.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "subj", "study", "split", "experiment", "participant"])
        for key in sorted(ours):
            split, exp, part = found.get(key, ("absent", "", ""))
            w.writerow([*key, split, exp, part])

    print("Psych-101 IGT transcripts read:", dict(counts))
    table = defaultdict(lambda: defaultdict(int))
    for key in ours:
        table[key[0]][found.get(key, ("absent",))[0]] += 1
    for sample, row in table.items():
        print(sample, dict(row))                               # how many Centaur trained on, held out, never saw


if __name__ == "__main__":
    main()
