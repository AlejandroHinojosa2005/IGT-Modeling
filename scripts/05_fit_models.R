# Fit models to all Horstmann participants (multi-startmaximum likelihood)
# urrently: PVL-Delta, VPP. ORL will be added once it's script exist.

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

# VPP: fit every participant
# 50 starts: VPP's starts rarely agree, so 20 can miss best fit
# Slow: roughly 45-60 min.
source("scripts/03_vpp.R")
n_subj <- nrow(choice_h)

set.seed(2026)
results_vpp <- data.frame(
  subj = 1:n_subj,
  A=NA, alpha=NA, cons=NA, lambda=NA, epP=NA, epN=NA, K=NA, w=NA, nll=NA, conv=NA
)

for (i in 1:n_subj) {
  f <- fit_vpp(choice_h[i, ], net_h[i, ], n_starts = 50)
  results_vpp[i, c("A", "alpha", "cons", "lambda", "epP", "epN", "K", "w")] <- f$par
  results_vpp$nll[i] <- f$value
  results_vpp$conv[i] <- f$convergence
}

saveRDS(results_vpp, "output/vpp_fit_horstmann.rds")

# ORL: fit every participant (20 starts; starts agreed for participant 1)
source("scripts/04_orl.R")
n_subj <- nrow(choice_h)

set.seed(2026)
results_orl <- data.frame(
  subj = 1:n_subj,
  Arew = NA, Apun = NA, Kprime = NA, betaF = NA, betaP = NA, nll = NA, conv = NA
)

for (i in 1:n_subj) {
  f <- fit_orl(choice_h[i, ], net_h[i, ])
  results_orl[i, c("Arew", "Apun", "Kprime", "betaF", "betaP")] <- f$par
  results_orl$nll[i] <- f$value
  results_orl$conv[i] <- f$convergence
}

saveRDS(results_orl, "output/orl_fit_horstmann.rds")

# Check: is the betaP limit (+/-10) costing fit?
# Refit only participants whose betaP hit the limit, with limits
source("scripts/04_orl.R")
results_orl <- readRDS("output/orl_fit_horstmann.rds")
pinned <- which(abs(results_orl$betaP) >= 9.999)
nll_wide <- rep(NA, length(pinned))

set.seed(2026)
for (n in 1:length(pinned)) {
  i <- pinned[n]
  nll_wide[n] <- fit_orl(choice_h[i, ], net_h[i, ], beta_lim = 30)$value
}

improve <- results_orl$nll[pinned] - nll_wide # positive = wider limit fits better
summary(improve)
sum(improve > 1)