# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), exploratory: data WITH switching.
# Script 44 made players who mix the eight models on every choice and never switch; the switching model found no
# switching in them (delta = 1 in every group). This script makes two sets of players who do have a model in charge,
# to check that the same pipeline (scripts 45-47) finds switching when it is there, and at the right rate:
#   "sim_switch": one model is in charge at a time. On each trial, with probability delta the model in charge is
#                 redrawn from the player's weights (it can be redrawn as the same model); otherwise it stays.
#                 delta is the rate fit to all 162 Horstmann participants by script 37 (0.053, a redraw about every
#                 19 trials).
#   "sim_types":  one model is in charge for the whole session (delta = 0), drawn once from the player's weights.
#                 This is the second simulation listed under exploratory analyses in the registration.
# Everything else is as in script 44: each Horstmann participant's full-data parameters for every model, the
# shrinkage weights 0.7 x own + 0.3 x shared as the player's weights, outcomes drawn from the chosen deck's cards in
# the Horstmann data, and all eight models updating on every choice. The model in charge gives the choice.
# One player per Horstmann participant in each set (N = 162 each, the size of the sample the rates come from).
# The model in charge on every trial is saved (in_charge), so later scripts can compare what they infer with the truth.
# Usage (from the project folder): Rscript scripts/44b_switch_simulation.R
source("scripts/43_prereg_data.R")
source("scripts/43b_prereg_models.R")

fit_files <- list(orl = "output/orl_fit_horstmann.rds", variant = "output/orl_wl_fit_horstmann.rds",
                  orl_alpha = "output/alpha/alpha_orl.rds", variant_alpha = "output/alpha/alpha_wl.rds",
                  vse = "output/vse/fit_vse.rds", vse_split = "output/vse/fit_vse_split.rds",
                  eef = "output/vse/fit_eef.rds", wsls = "output/vse/fit_wsls.rds")
par_cols <- c(lapply(models, `[[`, "names"), list(wsls = c("p_stay_win", "p_shift_loss")))
pars <- lapply(model_names, function(m) as.matrix(readRDS(fit_files[[m]])[, par_cols[[m]]]))
names(pars) <- model_names
stopifnot(all(sapply(pars, nrow) == 162), !anyNA(unlist(pars)))

s <- readRDS("output/stacking.rds")
stopifnot(identical(colnames(s$w_person), model_names), identical(names(s$w_group), model_names))
w_mix <- 0.7 * s$w_person + 0.3 * matrix(s$w_group, 162, 8, byrow = TRUE)
delta_h <- readRDS("output/switching.rds")$fit_all$delta     # switching rate fit to all Horstmann participants
delta_h

simulate_player <- function(i, delta) {
  st <- lapply(model_names, function(m) steps_of[[m]](pars[[m]][i, ]))
  state <- lapply(st, function(x) x$start())
  choice <- win <- loss <- z <- rep(NA, 100)
  z_t <- sample(1:8, size = 1, prob = w_mix[i, ])
  for (t in 1:100) {
    if (t > 1 && runif(1) < delta) z_t <- sample(1:8, size = 1, prob = w_mix[i, ])
    z[t] <- z_t
    p <- st[[z_t]]$probs(state[[z_t]])
    j <- sample(1:4, size = 1, prob = p)
    k <- sample(nrow(deck_draws[[j]]), size = 1)
    choice[t] <- j; win[t] <- deck_draws[[j]][k, 1]; loss[t] <- deck_draws[[j]][k, 2]
    state <- lapply(seq_along(st), function(m) st[[m]]$update(state[[m]], j, win[t], loss[t]))
  }
  list(choice = choice, win = win, loss = loss, z = z)
}

make_set <- function(name, delta, seed) {
  set.seed(seed)
  sims <- lapply(1:162, simulate_player, delta = delta)
  d <- list(choice = t(sapply(sims, `[[`, "choice")), win = t(sapply(sims, `[[`, "win")),
            loss = t(sapply(sims, `[[`, "loss")), study = rep(name, 162), scheme = rep(2, 162),
            source = 1:162, w_mix = w_mix, delta = delta, in_charge = t(sapply(sims, `[[`, "z")))
  saveRDS(d, sprintf("output/prereg/%s_data.rds", name))
  d
}
sw <- make_set("sim_switch", delta_h, 2028)
ty <- make_set("sim_types", 0, 2029)

# Description only: redraws that changed the model in charge, and the share of good-deck choices
c(changes_per_player = mean(rowSums(sw$in_charge[, -1] != sw$in_charge[, -100])),
  types_changes = sum(ty$in_charge[, -1] != ty$in_charge[, -100]))
h <- load_set("horstmann")
round(c(sim_switch = mean(sw$choice >= 3), sim_types = mean(ty$choice >= 3), horstmann = mean(h$choice >= 3)), 3)
