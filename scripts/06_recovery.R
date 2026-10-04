# Parameter recovery: simulate participants with known ("true") parameters,
# refit them with the same multi-start procedure, & compare.
# True values = the 162 fitted Horstmann estimates, so the simulated
# participants have realistic parameter combinations.

source("scripts/02_pvl_delta.R")
results_pvl <- readRDS("output/pvl_delta_fit_horstmann.rds")

# Each deck's pool of observed net outcomes, for the simulator draw from
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2], net_h[choice_h == 3], net_h[choice_h == 4])

# PVL-Delta recovery
true_par <- as.matrix(results_pvl[, c("A", "alpha", "cons", "lambda")])
n_sim <- nrow(true_par)
recovered <- matrix(NA, nrow=n_sim, ncol = 4)
colnames(recovered) <- colnames(true_par)

set.seed(2026)
for (k in 1:n_sim) {
  sim <- sim_pvl_delta(true_par[k, ], deck_outcomes)
  f <- fit_pvl_delta(sim$choice, sim$outcome)
  recovered[k, ] <- f$par
}

saveRDS(list(true = true_par, recovered = recovered),
"output/pvl_delta_recovery.rds")

# How well does each parameter come back?
round(diag(cor(true_par, recovered)), 2)

par(mfrow = c(2,2))
for (k in 1:4) {
  plot(true_par[, k], recovered[, k], xlab = "True", ylab = "Recovered",
       main=colnames(true_par)[k])
  abline(0,1)
}