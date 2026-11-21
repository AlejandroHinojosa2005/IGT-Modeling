# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), H6: data with no switching.
# Each Horstmann participant gives two simulated players (N = 324, seed 2027). Every model uses that participant's
# full-data maximum-likelihood parameters (the main fits in Hinojosa, 2026: scripts 05, 11, 28 and 32); WSLS uses
# their stay-after-win and shift-after-loss proportions. On each trial:
#   1. Each of the eight models gives its probability for each deck, given the simulated history so far.
#   2. The choice is drawn from the mix of those probabilities, weighted by that participant's shrinkage weights:
#      0.7 x their own stacking weights + 0.3 x the shared weights (script 36, from all held-out choices).
#   3. The outcome is a (win, loss) card drawn at random from the chosen deck's cards in the Horstmann data
#      (scheme 2; deck_draws, as in the simulators of scripts 31 and 39).
#   4. All eight models update on that choice and outcome.
# The mix is drawn fresh on every trial, so no model stays in charge (delta = 1).
# Simulated player 2i - 1 and 2i both come from Horstmann participant i.
# Usage (from the project folder): Rscript scripts/44_null_simulation.R
source("scripts/43_prereg_data.R")
source("scripts/43b_prereg_models.R")

# Each participant's full-data parameters, in the order each model's functions take them
fit_files <- list(orl = "output/orl_fit_horstmann.rds", variant = "output/orl_wl_fit_horstmann.rds",
                  orl_alpha = "output/alpha/alpha_orl.rds", variant_alpha = "output/alpha/alpha_wl.rds",
                  vse = "output/vse/fit_vse.rds", vse_split = "output/vse/fit_vse_split.rds",
                  eef = "output/vse/fit_eef.rds", wsls = "output/vse/fit_wsls.rds")
par_cols <- c(lapply(models, `[[`, "names"), list(wsls = c("p_stay_win", "p_shift_loss")))
pars <- lapply(model_names, function(m) as.matrix(readRDS(fit_files[[m]])[, par_cols[[m]]]))
names(pars) <- model_names
stopifnot(all(sapply(pars, nrow) == 162), !anyNA(unlist(pars)))

# Shrinkage weights: 0.7 x own + 0.3 x shared, from script 36
s <- readRDS("output/stacking.rds")
stopifnot(identical(colnames(s$w_person), model_names), identical(names(s$w_group), model_names))
w_mix <- 0.7 * s$w_person + 0.3 * matrix(s$w_group, 162, 8, byrow = TRUE)
stopifnot(all(abs(rowSums(w_mix) - 1) < 1e-8))

simulate_player <- function(i) {
  st <- lapply(model_names, function(m) steps_of[[m]](pars[[m]][i, ]))
  state <- lapply(st, function(x) x$start())
  choice <- win <- loss <- rep(NA, 100)
  for (t in 1:100) {
    P <- sapply(seq_along(st), function(k) st[[k]]$probs(state[[k]]))   # 4 decks x 8 models
    p <- as.vector(P %*% w_mix[i, ])
    j <- sample(1:4, size = 1, prob = p)
    k <- sample(nrow(deck_draws[[j]]), size = 1)
    choice[t] <- j; win[t] <- deck_draws[[j]][k, 1]; loss[t] <- deck_draws[[j]][k, 2]
    state <- lapply(seq_along(st), function(m) st[[m]]$update(state[[m]], j, win[t], loss[t]))
  }
  list(choice = choice, win = win, loss = loss)
}

set.seed(2027)
sims <- list()
for (i in 1:162) for (r in 1:2) sims[[length(sims) + 1]] <- simulate_player(i)
d <- list(choice = t(sapply(sims, `[[`, "choice")), win = t(sapply(sims, `[[`, "win")),
          loss = t(sapply(sims, `[[`, "loss")), study = rep("simulated", 324), scheme = rep(2, 324),
          source = rep(1:162, each = 2), w_mix = w_mix)

# Description only: share of choices from the good decks (C and D), simulated vs. real Horstmann
h <- load_set("horstmann")
round(c(simulated = mean(d$choice >= 3), horstmann = mean(h$choice >= 3)), 3)

dir.create("output/prereg", showWarnings = FALSE)
saveRDS(d, "output/prereg/null_sim_data.rds")
