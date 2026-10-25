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






