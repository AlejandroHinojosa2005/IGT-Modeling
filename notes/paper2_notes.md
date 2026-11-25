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

**Results (run finished November 24, 2026, 21:24).** Every test that was met with 440 is met with 438, and almost
nothing moves. z for 440 → 438:

- Registered tests (`output/prereg/log58_duplicates.txt`): H1 12.5 → 12.5, H1b 11.6 → 11.6, H2 7.3 → 7.4, H3 6.6 → 6.7,
  H4 6.4 → 5.4, H5 4.8 → 4.7, H7 6.5 → 6.3. H4 is the only one that moves by more than 0.3. Most likely the reason is
  that its test uses random draws (the seeded shuffles and groupings in script 47), and with two fewer people the same
  seeds give different draws. Worth checking before the paper says so.
- NOVA's tests (`output/hybrid/log50d_replication438.txt`): P1 −1.0, P2 5.4, P3 3.9, P4 4.9, all met as before.
  NOVA's weight in the shared mix: .330.
- Pooled (`output/pooled/log54c_replication438.txt`): Q1 4.3 → 4.3 (+1184.2), Q2 6.4 → 6.4 (+458.1),
  Q3 −9.7 → −9.6 (−339.5), Q4 −4.1 → −4.0 (−168.9), Q5 −8.3 → −8.2 (−634.9). NOVA's group estimates: Q0 2.08 (was
  2.09), D 0.96 (same).
- Forward, F1, the mix of eight beats ORL + α (`output/forward/log57b_replication438.txt` and
  `output/forward/log57b_horstmann_prior-replication438.txt`): Horstmann −477.1, SE 64.8, z −7.4 (was −475.2, z −7.4);
  replication −1215.8, SE 79.4, z −15.3 (was −1226.6, z −15.4). The mix of five plus NOVA still beats the mix of eight
  (z 4.7 and 5.0, unchanged), and NOVA is still the best single model.

**For the paper.** One or two sentences in Methods (the pool contains two repeated records; Psych-101 keeps one of
each), the registered results as registered, and the N = 438 results as a robustness check next to robustness check 5.

## Also from the Psych-101 match (November 24, 2026)

- Psych-101 has 577 IGT records: the whole pool (617) except Maia's study (40), which is not in Psych-101 at all.
- Of our 602: 498 in Centaur's training set, 64 in its held-out set (17 Horstmann, 47 replication), and Maia's 40 absent.
  Centaur is scored on the 104 it never saw (`centaur/score_centaur.py` scores the held-out and the absent).

## Centaur against the models (November 24, 2026)

Centaur 70B (Binz et al., 2025) run on a rented RunPod pod (2 × A100 80GB, about 30 minutes, under $5), scoring every
choice of the 104 participants it never trained on (`centaur/run_centaur.sh`; scores in `centaur/centaur70b_scores.csv`,
never committed). Compared with our models' forward predictions on the same people and choices (`scripts/57b`, with
the scores file). Exploratory: no prediction was written for this comparison. Probability given to the choice actually
made, trials 2-100:

| Participants | Centaur | Mix of eight | NOVA | ORL + α | Mix vs. Centaur |
|---|---|---|---|---|---|
| Horstmann, held out of Psych-101 (17) | .439 | .384 | .380 | .368 | z = −5.8; mix better for 0 of 17 |
| Replication, held out (47) | .400 | .359 | .356 | .350 | z = −6.7; 7 of 47 |
| Maia, a study not in Psych-101 at all (40) | .415 | .377 | .374 | .363 | z = −8.7; 1 of 40 |

Centaur predicts the next choice better than every model and both mixes, for nearly everyone, including a whole study
it never saw. Points for the paper: Centaur learned from about 500 other people doing this exact task (same decks and
payoff schedules), while our models know the task only through their equations and a group prior from the other sample.
So the gap measures how much predictable structure the models miss, not which explanation is right. NOVA and the mixes
come closest. This is the case for the small-RNN benchmark (Ji-An, Benna & Mattar, 2025): find out what Centaur captures
that the models do not.

(Numbers computed in Claude's workspace with the repository's scripts; rerun on the Mac for the record:
`Rscript scripts/57b_forward_compare.R horstmann centaur/centaur70b_scores.csv` and the same with `replication`.)
