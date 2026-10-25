# ORL variant copied into this script: fit every participant
# (20 starts, one of them at that participant's ORL fit)

source("scripts/10_orl_variant.R")
n_subj <- nrow(choice_h)

set.seed(2026)
orl <- readRDS("output/orl_fit_horstmann.rds")
results_orl_wl <- data.frame(
  subj = 1:n_subj,
  Arew = NA, Apun = NA, Kprime = NA, betaF = NA, betaPw = NA, betaPl = NA, nll = NA, conv = NA
)

for (i in 1:n_subj) {
  p <- unlist(orl[i, c("Arew", "Apun", "Kprime", "betaF", "betaP", "betaP")])
  f <- fit_orl_wl(choice_h[i, ], net_h[i, ], start_par = p)
  results_orl_wl[i, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")] <- f$par
  results_orl_wl$nll[i] <- f$value
  results_orl_wl$conv[i] <- f$convergence
}

# Check: 3 fits ended with conv = 52 (both betaPs near the -10 limit, K' near 0).
# Refit those with 50 starts & keep the refits: nll improved by at most 0.08
# (participatn 75), enough to move 75 past the AIC bar (83 -> 84 variant wins).

bad <- which(results_orl_wl$conv != 0)
set.seed(2026)
for (i in bad) {
  p <- unlist(orl[i, c("Arew","Apun","Kprime","betaF","betaP","betaP")])
  f <- fit_orl_wl(choice_h[i, ], net_h[i, ], n_starts = 50, start_par = p)
  results_orl_wl[i, c("Arew","Apun","Kprime","betaF","betaPw","betaPl")] <- f$par
  results_orl_wl$nll[i] <- f$value
  results_orl_wl$conv[i] <- f$convergence
}

saveRDS(results_orl_wl, "output/orl_wl_fit_horstmann.rds")


# Checks & comparison with ORL and VPP
all(results_orl_wl$nll <= orl$nll + 1e-8)   # never worse than ORL
table(results_orl_wl$conv)

vpp <- readRDS("output/vpp_fit_horstmann.rds")
c(orl_aic = sum(2 * orl$nll + 10), wl_aic = sum(2 * results_orl_wl$nll + 12),
  vpp_aic = sum(2 * vpp$nll + 16))
c(orl_bic = sum(2 * orl$nll + 5 * log(100)), wl_bic = sum(2 * results_orl_wl$nll + 6 * log(100)),
  vpp_bic = sum(2 * vpp$nll + 8 * log(100)))

sum(orl$nll - results_orl_wl$nll > 1)              # variant wins vs. ORL, AIC
sum(orl$nll - results_orl_wl$nll > log(100) / 2)   # variant wins vs. ORL, BIC
sum(results_orl_wl$betaPw > results_orl_wl$betaPl) # prediction 2

# Test ran that asked if weather spiting perspiration leaves ORL's other parameters unchanged
# Shows that A_pun is absorbing the switching after losses
p <- c("Arew", "Apun", "Kprime", "betaF")
round(diag(cor(orl[, p], wl[, p])), 2)          # same people: does each parameter keep its ranking?
round(colMeans(wl[, p]) - colMeans(orl[, p]), 2) # average shift

















