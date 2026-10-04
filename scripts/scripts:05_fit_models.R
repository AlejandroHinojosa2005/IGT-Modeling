# Fit models to all Horstmann participants (multi-startmaximum likelihood)
# urrently: PVL-Delta. VPP & ORL will be added once their scripts exist.

source("scripts/02_pvl_delta.R")

# PVL-Delta: fit every participant (take a few minutes)
set.seed(2026)
n_subj <- nrow(choice_h)

results_pvl <- data.frame(
  subj = 1:n_subj,
  A = NA,
  alpha = NA,
  cons = NA,
  lambda = NA,
  nll = NA,
  conv = NA
)

for (i in 1:n_subj) {
  f <- fit_pvl_delta(choice_h[i, ], net_h[i, ])
  results_pvl[i, c("A", "alpha", "cons", "lambda")] <- f$par
  results_pvl$nll[i] <- f$value
  results_pvl$conv[i] <- f$convergence
}

saveRDS(results_pvl, "output/pvl_delta_fit_horstmann.rds")