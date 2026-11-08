# Null check for prediction 1: if nobody's alpha were below 1, how often would the fits say it is?
# Simulate from each participant's fitted variant (alpha = 1 built in), refit with the same
# procedure as the real fits (ORL -> variant -> variant + alpha), and count alpha < 1.
# Usage (from the project folder):  Rscript scripts/29b_alpha_null.R [quick]
args <- commandArgs(trailingOnly = TRUE)
quick <- length(args) >= 1 && args[1] == "quick"
tag <- if (quick) "_quick" else ""

source("scripts/27_alpha_model.R")
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2],
                      net_h[choice_h == 3], net_h[choice_h == 4])

wl <- readRDS("output/orl_wl_fit_horstmann.rds")
true_par <- as.matrix(wl[, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")])
n <- if (quick) 3 else nrow(true_par)
alpha_hat <- rep(NA, n)

set.seed(2026)
for (i in 1:n) {
  sim <- sim_orl_wl(true_par[i, ], deck_outcomes)
  f_orl <- fit_orl(sim$choice, sim$outcome)
  f_wl <- fit_orl_wl(sim$choice, sim$outcome, start_par = c(f_orl$par, f_orl$par[5]))
  f <- fit_alpha(orl_wl_alpha_nll, sim$choice, sim$outcome, lower_wl_alpha, upper_wl_alpha,
                 start_par = c(f_wl$par, 1))
  alpha_hat[i] <- f$par[7]
}

saveRDS(alpha_hat, sprintf("output/alpha/alpha_null%s.rds", tag))
summary(alpha_hat)
sum(alpha_hat < 1)          # real data: 130 of 162
sum(alpha_hat <= 0.011)     # real data: 39 at the lower limit