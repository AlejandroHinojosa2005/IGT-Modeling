# ORL variant parameter recovery: simulate from the fitted variant parameters,
# refit with the same procedure as the real fit (ORL fit as one start + 20 starts)
source("scripts/10_orl_variant.R")
results_orl_wl <- readRDS("output/orl_wl_fit_horstmann.rds")

deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2], net_h[choice_h == 3], net_h[choice_h == 4])

wl_names <- c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")
true_wl <- as.matrix(results_orl_wl[, wl_names])
recovered_wl <- matrix(NA, nrow = nrow(true_wl), ncol = 6)

colnames(recovered_wl) <- wl_names

set.seed(2026)
for (k in 1:nrow(true_wl)) {
  sim <- sim_orl_wl(true_wl[k, ], deck_outcomes)
  f_orl <- fit_orl(sim$choice, sim$outcome)
  f <- fit_orl_wl(sim$choice, sim$outcome, start_par = c(f_orl$par, f_orl$par[5]))
  recovered_wl[k, ] <- f$par
}

saveRDS(list(true = true_wl, recovered = recovered_wl), "output/orl_wl_recovery.rds")

round(diag(cor(true_wl, recovered_wl)), 2)
round(diag(cor(true_wl, recovered_wl, method = "spearman")), 2)

par(mfrow = c(2,3), mar = c(4,4,2,1))
for (k in 1:6) {
  plot(true_wl[, k], recovered_wl[, k], xlab = "True", ylab = "Recovered", main = wl_names[k])
  abline(0,1)
}

# How many estimates sit at a limit? (out of 162)
# True = real-data variant fits; recovered = fits to the simulated data
lower_wl <- c(.01,.01,0,-10,-10,-10)
upper_wl <- c(1,1,5,10,10,10)

cat("True values at a limit:\n")
for (k in 1:6) {
  at_limit <- true_wl[, k] <= lower_wl[k] + .001 | true_wl[, k] >= upper_wl[k] - .001
  cat(wl_names[k], sum(at_limit), "\n")
}

cat("Recovered values at a limit:\n")
for (k in 1:6) {
  at_limit <- recovered_wl[, k] <= lower_wl[k] + .001 | recovered_wl[, k] >= upper_wl[k] - .001
  cat(wl_names[k], sum(at_limit), "\n")
}

# Can the split itself be recovered? (prediction 2 depends on it)
true_diff <- true_wl[, "betaPw"] - true_wl[, "betaPl"]
rec_diff <- recovered_wl[, "betaPw"] - recovered_wl[, "betaPl"]
round(cor(true_diff, rec_diff), 2)
round(cor(true_diff, rec_diff, method = "spearman"), 2)
mean(sign(true_diff) == sign(rec_diff))


# Does the fitting tilt the difference?
mean(rec_diff - true_diff)
c(mean(true_diff > 0), mean(rec_diff > 0))

















