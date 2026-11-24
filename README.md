# Reinforcement-learning models of the Iowa Gambling Task

This started as a class project for a cognition and AI course. I wanted to see whether I could rebuild the
Outcome-Representation Learning model (ORL; Haines, Vassileva & Ahn, 2018) by hand, and whether anything about it could
be improved. I implemented three reinforcement-learning models of the Iowa Gambling Task (PVL-Delta, VPP, and ORL) in R
from their published equations and fit them to 162 healthy adults by maximum likelihood. Each answer raised the next
question, and the project grew into two papers.

**Status:** Personal research project, ongoing (last updated November 23, 2026). A plain-language version of the main
results, with a few interactive figures, is on my website: [alejandrohinojosa.com](https://alejandrohinojosa.com).

## What I found

**1. Outcome curvature (paper 1).** Of the original three models, ORL fit best, but it under-produced win-stay/lose-shift.
I built a variant whose perseveration depends on whether the last outcome was a gain or a loss, and it fixed that. Then
I added diminishing sensitivity (α, the utility curvature PVL-Delta uses and ORL dropped): most people feel large wins
and losses less than in proportion to their size. ORL + α improved the fit more than my own split did, and once α was in,
the split no longer helped predict held-out choices. In a cross-validated comparison of eight models, ORL + α predicted
held-out choices best. On a held-out sample of 70 people, all five predictions I wrote down before fitting were met, and
in the 440 participants of the other studies in the data pool it was again the best of the eight single models.

**2. People mix strategies (paper 2).** No single model was best for everyone, so I asked whether people mix strategies.
In the first sample, a weighted mix of the eight models predicted held-out choices far better than ORL + α. I
preregistered ten predictions from those results ([OSF](https://doi.org/10.17605/OSF.IO/U86SY)) and tested them in 440 new
participants from eight other studies. All ten were met. The shared mix beat ORL + α for 328 of 440 people (z = 12.5),
and personal weights pulled toward the shared ones did better still (z = 7.3). A model in which one strategy is in charge
for a run of trials and then hands off beat a fixed mix (z = 6.6), and it needed the choices in their real order
(z = 6.4). On simulated players who mix strategies but never switch, the method found no switching. I call the whole
procedure (held-out choices, mixing, shrinkage toward the group, and the switching tests) Seldon.

**3. Two models, one term.** VSE (Ligneul, 2019) and ORL are treated as rival accounts of the same behavior: VSE says
people explore decks they haven't chosen in a while, and ORL says they persevere on decks they chose recently. They are
the same term. Once their parameters are mapped onto each other, VSE's exploration bonus is exactly ORL's perseveration
with the sign reversed, a pull away from recently chosen decks. Script 56 checks this on every trial of all 602
participants (the differences are about 10⁻¹²) and shows the one caveat, about starting values. What really separates
the two families is how they learn from outcomes, not how they explore.

**4. NOVA.** Putting the pieces together gave NOVA: **N**ovelty (an optimistic starting value for decks not yet tried),
**O**utcome curvature (α), **V**alue decay (for decks not chosen), and **A**lternation (VSE's exploration, written as
negative perseveration). Fit one person at a time, NOVA alone did not beat ORL + α, because a few people's fits collapsed,
but it became the most useful member of the mix. With predictions written in advance, adding it improved the mix in the
440 (z = 5.4), and a mix of just five models built around it beat all eight (z = 3.9). When each person's parameters are
pulled toward the group, NOVA beats ORL + α on its own in the first sample (z = 4.6). The same test on the 440 is running
now, with its predictions committed first (`predictions/pooled_replication.md`).

## What's coming

- **The papers.** Paper 1 (ORL + α) is going up on PsyArXiv. Paper 2 (mixing, switching and NOVA) is being written.
- **Hierarchical versions of every model.** hBayesDM has four Iowa Gambling Task models; `stan/hier/` has Stan versions
  of all eight models here, seven of them new. Their likelihoods match the R functions to about 10⁻¹¹, and parameter
  recovery is running.
- **A benchmark against Centaur.** Centaur (Binz et al., 2025) is a language model fine-tuned on about 10 million choices
  from 160 experiments, including this data pool. Scored only on people it never trained on, against my
  models scored the same way (each choice predicted only from the choices before it), it gives a rough upper bound on how
  predictable these choices are. The scripts are ready; it needs GPUs (`centaur/`).
- **Further out:** hierarchical model selection (Piray et al., 2019), small recurrent neural networks that can be read
  like models, tasks built to separate exploration from perseveration (like the Horizon task), and EEG.

## Background

Back in 1994, Bechara et al. wanted to figure out a way to assess decision making
deficits in people with ventromedial prefrontal damage who despite normal intelligence,
displayed poor decision making. To do this they developed the Iowa Gambling Task.
The scheme goes as follows:

The participant would be given $2000 of facsimile money. They would then be shown four
decks of cards; they choose one deck and pick one card from said deck at a time. Each
card gives a reward (rew), with some also carrying a loss/punishment (pun). Their
instructions are to make as much money as possible while not being told deck properties.
Over every ten cards, decks A & B both have a -250 expected outcome (losing decks),
while decks C & D have expected outcome of +250 (winning decks). The paired decks varied
in both how often pun was doled out & how much any given pun is worth. For losing decks,
rew was fixed at +100 while pun varied between -150 to -350 for deck A and was fixed at
-1250 for deck B. For winning decks, rew was fixed at +50 while pun varied between -25
to -75 for deck C and was fixed at -250 for deck D. Every card has a rew, but decks A & C
carry a pun on 5 of every 10 cards, while B & D carry one on only 1 of every 10. Participants
are assessed on ability to favor winning decks over the losing ones.

Summary scores like the net score (C+D-A-B) can't tell apart people who reach the same score
through different means, e.g. faster vs. slower learning or more versus less consistent choosing.
To address this, computational modeling is used. These models assume that following each trial,
participants update their expectations about the decks. The models I coded build on the Expectancy
Valence model (Busemeyer & Stout, 2002) and Prospect Valence Learning (PVL; Ahn et al., 2008).

## What's in this repository

| Where | What |
|---|---|
| `scripts/01`–`42` | Paper 1: the three models, fits, recovery, the variant, α, the newer models, cross-validation and the held-out test |
| `scripts/43`–`49`, `51`, `predictions/` | The preregistered replication and its exploratory follow-ups |
| `scripts/50`–`52` | NOVA and the smaller mixes |
| `scripts/54` | Pooled fits (each person pulled toward the group) |
| `scripts/55`, `stan/hier/` | Hierarchical Stan versions of the eight models |
| `scripts/56` | The check that VSE's exploration is ORL's perseveration |
| `scripts/53`, `57`, `centaur/` | The Centaur benchmark |
| `kodiak/` | Job files for Baylor's Kodiak cluster |
| `output/` | Saved results and logs, so later scripts don't refit |
| `notes/lab_notebook.md` | The full notebook of the paper 1 work, with every number (this README until November 23, 2026) |

## How to run

R 4.6.1. Package versions are recorded with renv: open the project in RStudio and run `renv::restore()`. The data aren't
included: download `IGTdata.rdata` from <https://osf.io/8t7rm> (Steingroever et al., 2015) and put it in `data-raw/`.
Scripts 01–42 are listed, with run times, under *How to run* in `notes/lab_notebook.md`. Everything from script 43 on is
in `REPLICATION.md`, in the order it was run, along with the commits that show what was written down before it was run.

## How this project was built

Goals:
- Teach myself how to code in R
- Understand how scientists use computational methods in research relating to Computational Psychology

Assistance:
I used Anthropic's Claude (Opus 5.5) as a tutor: it explained the models' equations and R syntax,
reviewed my code for errors, and helped edit this write-up. I implemented the three models from their published equations
(Ahn et al., 2008; Worthy, Pang & Byrne, 2013; Haines, Vassileva & Ahn, 2018), following hBayesDM's outcome scaling and consulting
its ORL code (Ahn, Haines & Zhang, 2017). I wrote and ran the code myself in RStudio.

## References

- Ahn, W.-Y., Busemeyer, J. R., Wagenmakers, E.-J., & Stout, J. C. (2008). Comparison of decision learning models using the generalization criterion method. *Cognitive Science, 32*(8), 1376–1402. https://doi.org/10.1080/03640210802352992
- Ahn, W.-Y., Haines, N., & Zhang, L. (2017). Revealing neurocomputational mechanisms of reinforcement learning and decision-making with the hBayesDM package. *Computational Psychiatry, 1*, 24–57. https://doi.org/10.1162/CPSY_a_00002
- Bechara, A., Damasio, A. R., Damasio, H., & Anderson, S. W. (1994). Insensitivity to future consequences following damage to human prefrontal cortex. *Cognition, 50*(1–3), 7–15. https://doi.org/10.1016/0010-0277(94)90018-3
- Binz, M., et al. (2025). A foundation model to predict and capture human cognition. *Nature*. https://doi.org/10.1038/s41586-025-09215-4
- Busemeyer, J. R., & Stout, J. C. (2002). A contribution of cognitive decision models to clinical assessment: Decomposing performance on the Bechara gambling task. *Psychological Assessment, 14*(3), 253–262. https://doi.org/10.1037/1040-3590.14.3.253
- Haines, N., Vassileva, J., & Ahn, W.-Y. (2018). The Outcome-Representation Learning model: A novel reinforcement learning model of the Iowa Gambling Task. *Cognitive Science, 42*(8), 2534–2561. https://doi.org/10.1111/cogs.12688
- Ligneul, R. (2019). Sequential exploration in the Iowa gambling task: Validation of a new computational model in a large dataset of young and old healthy participants. *PLOS Computational Biology, 15*(6), e1006989. https://doi.org/10.1371/journal.pcbi.1006989
- Piray, P., Dezfouli, A., Heskes, T., Frank, M. J., & Daw, N. D. (2019). Hierarchical Bayesian inference for concurrent model fitting and comparison for group studies. *PLOS Computational Biology, 15*(6), e1007043. https://doi.org/10.1371/journal.pcbi.1007043
- Steingroever, H., et al. (2015). Data from 617 healthy participants performing the Iowa Gambling Task: A "many labs" collaboration. *Journal of Open Psychology Data, 3*(1), e5. https://doi.org/10.5334/jopd.ak
- Worthy, D. A., Pang, B., & Byrne, K. A. (2013). Decomposing the roles of perseveration and expected value representation in models of the Iowa Gambling Task. *Frontiers in Psychology, 4*, 640. https://doi.org/10.3389/fpsyg.2013.00640
- Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using stacking to average Bayesian predictive distributions (with discussion). *Bayesian Analysis, 13*(3), 917–1007. https://doi.org/10.1214/17-BA1091

The full reference list for the paper 1 work is in `notes/lab_notebook.md`.

## License

MIT License: see `LICENSE`. The exception is the Stan files in `stan/`: they are modified versions of hBayesDM's
`igt_orl.stan` and are distributed under the GNU General Public License v3 (GPL-3), as their headers state.
The data are not included and keep their own OSF license (CC BY-SA 4.0).
