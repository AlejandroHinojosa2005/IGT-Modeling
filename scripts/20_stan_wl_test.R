# Validation of the Stan (hierarchical) version of the ORL variant before any real run:
# simulate 30 participants with known parameters (R simulator), fit them with the
# Stan model, and check that the parameters come back.
library(cmdstanr)
source("scripts/10_orl_variant.R")   # sim_orl_wl, Horstmann data
results_orl_wl <- readRDS("output/orl_wl_fit_horstmann.rds")
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2], net_h[choice_h == 3], net_h[choice_h == 4])

# True parameters: the MLE variant fits of 30 randomly chosen Horstmann participants
set.seed(2026)
n_sim <- 30
idx <- sample(1:nrow(results_orl_wl), n_sim)
wl_names <- c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")
true_wl <- as.matrix(results_orl_wl[idx, wl_names])

# Simulate one 100-trial dataset per participant
choice_sim <- matrix(NA, n_sim, 100)
outcome_sim <- matrix(NA, n_sim, 100)
for (k in 1:n_sim) {
  sim <- sim_orl_wl(true_wl[k, ], deck_outcomes)
  choice_sim[k, ] <- sim$choice
  outcome_sim[k, ] <- sim$outcome
}

# Data in the format the Stan model expects (same as hBayesDM's preprocessing)
stan_data <- list(
  N = n_sim, T = 100, Tsubj = rep(100, n_sim),
  choice = matrix(as.integer(choice_sim), n_sim, 100),
  outcome = outcome_sim / 100,
  sign_out = sign(outcome_sim)
)

mod <- cmdstan_model("stan/igt_orl_wl.stan")
fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000,
                  adapt_delta = 0.95, seed = 2026, refresh = 500)

# Posterior mean of each participant's parameters (Stan's K is on the same 0-5 scale as Kprime)
stan_names <- c("Arew", "Apun", "K", "betaF", "betaPw", "betaPl")
rec_wl <- matrix(NA, n_sim, 6, dimnames = list(NULL, wl_names))
for (p in 1:6) {
  rec_wl[, p] <- colMeans(fit$draws(variables = stan_names[p], format = "draws_matrix"))
}
saveRDS(list(idx = idx, true = true_wl, recovered = rec_wl), "output/stan_wl_test.rds")

# Checks
fit$diagnostic_summary()                                    # divergences, treedepth problems
max(fit$summary(variables = stan_names)$rhat)               # should be <= 1.1
round(diag(cor(true_wl, rec_wl, method = "spearman")), 2)   # recovery, per parameter
true_diff <- true_wl[, "betaPw"] - true_wl[, "betaPl"]
rec_diff <- rec_wl[, "betaPw"] - rec_wl[, "betaPl"]
round(cor(true_diff, rec_diff, method = "spearman"), 2)     # the split
mean(sign(true_diff) == sign(rec_diff))                     # direction recovered






