# Parameter recovery for the variant + alpha: simulate from each participant's fitted
# parameters (script 28), then refit with the same procedure as the real fits
# (ORL, then the variant starting from ORL, then variant + alpha starting from the variant).
# Usage (from the project folder):  Rscript scripts/29_alpha_recovery.R [quick]
args <- commandArgs(trailingOnly = TRUE)
quick <- length(args) >= 1 && args[1] == "quick"
tag <- if (quick) "_quick" else ""

source("scripts/27_alpha_model.R")
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2],
                      net_h[choice_h == 3], net_h[choice_h == 4])

names_alpha <- c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl", "alpha")
fits <- readRDS(sprintf("output/alpha/alpha_wl%s.rds", tag))
true_par <- as.matrix(fits[, names_alpha])
recovered <- true_par * NA

set.seed(2026)
for (i in 1:nrow(true_par)) {
  sim <- sim_orl_wl_alpha(true_par[i, ], deck_outcomes)
  f_orl <- fit_orl(sim$choice, sim$outcome)
  f_wl <- fit_orl_wl(sim$choice, sim$outcome, start_par = c(f_orl$par, f_orl$par[5]))
  f <- fit_alpha(orl_wl_alpha_nll, sim$choice, sim$outcome, lower_wl_alpha, upper_wl_alpha,
                 start_par = c(f_wl$par, 1))
  recovered[i, ] <- f$par
}

saveRDS(list(true = true_par, recovered = recovered), sprintf("output/alpha/alpha_recovery%s.rds", tag))
round(diag(cor(true_par, recovered)), 2)
