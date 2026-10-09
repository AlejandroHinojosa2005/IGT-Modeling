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

# VPP recovery (slow: roughly an 20min with 50 starts)
source("scripts/03_vpp.R")
results_vpp <- readRDS("output/vpp_fit_horstmann.rds")

vpp_names <- c("A", "alpha", "cons", "lambda", "epP", "epN", "K", "w")
true_vpp <- as.matrix(results_vpp[, vpp_names])
recovered_vpp <- matrix(NA, nrow = nrow(true_vpp), ncol = 8)
colnames(recovered_vpp) <- vpp_names

set.seed(2026)
for (k in 1:nrow(true_vpp)) {
  sim <- sim_vpp(true_vpp[k, ], deck_outcomes)
  f <- fit_vpp(sim$choice, sim$outcome, n_starts = 50)
  recovered_vpp[k, ] <- f$par
}

saveRDS(list(true = true_vpp, recovered = recovered_vpp), "output/vpp_recovery.rds")

round(diag(cor(true_vpp, recovered_vpp)), 2)
round(diag(cor(true_vpp, recovered_vpp, method = "spearman")), 2)

par(mfrow = c(2, 4), mar = c(4,4,2,1))
for (k in 1:8) {
  plot(true_vpp[, k], recovered_vpp[, k], xlab = "True", ylab = "Recovered", main = vpp_names[k])
  abline(0,1)
}

# How many estimates sit at a parameter limit? (out of 162)
# True = real-data VPP fits: recovered = fits to the simulated

cat("True values at a limit: \n")
for (k in 1:8) {
  at_limit <- true_vpp[, k] <= lower_vpp[k] + 0.001 |
    true_vpp[, k] >= upper_vpp[k] - 0.001
  cat(vpp_names[k], sum(at_limit), "\n")
}

cat("Recovered vlaues at a limit: \n")
lower_vpp <- c(0.01, 0.01, 0.01, 0.01, -5, -5, 0.01, 0.01)
upper_vpp <- c(1, 2, 5, 10, 5, 5, 1, 1)
for (k in 1:8) {
  at_limit <- recovered_vpp[, k] <= lower_vpp[k] + 0.001 |
    recovered_vpp[, k] >= upper_vpp[k] - 0.001
  cat(vpp_names[k], sum(at_limit), "\n")
}