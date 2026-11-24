# Notes for paper 2

Things to remember when writing paper 2, in the order they came up.

## Two repeated records in the data pool (found November 24, 2026)

**What prompted it.** For the Centaur comparison we needed to know which of our participants Centaur was trained on.
Centaur was fine-tuned on Psych-101, whose Iowa Gambling Task data are this same many-labs pool, split by participant
into a training set and a held-out set. `centaur/match_psych101.py` finds each of our participants among the Psych-101
transcripts by their sequence of choices and outcomes. Every match was one-to-one except two: two Psych-101 training
transcripts each matched two of our participants.

**What it was.** Comparing all 602 of our records with each other (162 Horstmann, 440 replication) found exactly two
pairs that are identical on all 100 trials, choices, wins and losses:

- Kjome: participants 2 and 17 (rows of `load_set("replication")`)
- Wood: participants 298 and 299

These are almost certainly the same person entered twice in the pool. Psych-101 keeps one record of each pair. One
lookalike is not a repeat: Steingroever2011 participants 350 and 399 made identical choices (deck C on all 100 trials)
but their outcomes differ on 52 trials, so they are different people.

**What it affects.**

- Paper 1: nothing. The repeated records are in Kjome and Wood, not in Horstmann or the held-out sample (SteingroverInPrep).
- Paper 2: the preregistered replication (N = 440) counts each of the two people twice. The registered results stand as
  registered. The same analyses without the second record of each pair (N = 438, data set `replication438`) are an
  exploratory check, run with `scripts/58_duplicates_chain.sh`: the shared mix and switching (46, 47), the registered
  tests (58), NOVA's tests (50d), the pooled fits (54, 54b, 54c), and the forward predictions (57, 57b). Each model is
  fit to each person on their own, so those fits are reused; only the steps that combine participants are rerun.
- The Centaur comparison: nothing. All four records are in Centaur's training set, so none of them are scored.

**Results (to fill in after the run).** Registered tests on 438 vs. 440: `output/prereg/log58_duplicates.txt`.
NOVA: `output/hybrid/log50d_replication438.txt`. Pooled: `output/pooled/log54c_replication438.txt`.
Forward: `output/forward/log57b_replication438.txt` and `output/forward/log57b_horstmann_prior-replication438.txt`.

**For the paper.** One or two sentences in Methods (the pool contains two repeated records; Psych-101 keeps one of
each), the registered results as registered, and the N = 438 results as a robustness check next to robustness check 5.

## Also from the Psych-101 match (November 24, 2026)

- Psych-101 has 577 IGT records: the whole pool (617) except Maia's study (40), which is not in Psych-101 at all.
- Of our 602: 498 in Centaur's training set, 64 in its held-out set (17 Horstmann, 47 replication), and Maia's 40 absent.
  Centaur is scored on the 104 it never saw (`centaur/score_centaur.py` scores the held-out and the absent).
