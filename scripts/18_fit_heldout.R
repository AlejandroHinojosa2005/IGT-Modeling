# Held-out test: fit ORL, then the ORL variant, to the SteingroverInPrep sample (N = 70).
# Same procedure as Horstmann (scripts 05 & 11), as specified in the held-out test plan.
source("scripts/10_orl_variant.R")    # ORL & variant functions
source("scripts/17_load_heldout.R")   # choice_ho, net_ho
n_ho <- nrow(choice_ho)

# ORL: 20 starts per participant
set.seed(2026)
orl_ho <- data.frame(subj = 1:n_ho, Arew = NA, Apun = NA, Kprime = NA, betaF = NA, betaP = NA, nll = NA, conv = NA)
for (i in 1:n_ho) {
  f <- fit_orl(choice_ho[i, ], net_ho[i, ])
  orl_ho[i, c("Arew", "Apun", "Kprime", "betaF", "betaP")] <- f$par
  orl_ho$nll[i] <- f$value
  orl_ho$conv[i] <- f$convergence
}
saveRDS(orl_ho, "output/orl_fit_heldout.rds")

# Variant: 20 starts plus one at that participant's ORL fit
set.seed(2026)
wl_ho <- data.frame(subj = 1:n_ho, Arew = NA, Apun = NA, Kprime = NA, betaF = NA,
                    betaPw = NA, betaPl = NA, nll = NA, conv = NA)
for (i in 1:n_ho) {
  p <- unlist(orl_ho[i, c("Arew", "Apun", "Kprime", "betaF", "betaP", "betaP")])
  f <- fit_orl_wl(choice_ho[i, ], net_ho[i, ], start_par = p)
  wl_ho[i, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")] <- f$par
  wl_ho$nll[i] <- f$value
  wl_ho$conv[i] <- f$convergence
}

# Same rule as script 11: refit variant fits that didn't converge with 50 starts (keep the better fit)
bad <- which(wl_ho$conv != 0)
set.seed(2026)
for (i in bad) {
  p <- unlist(orl_ho[i, c("Arew", "Apun", "Kprime", "betaF", "betaP", "betaP")])
  f <- fit_orl_wl(choice_ho[i, ], net_ho[i, ], n_starts = 50, start_par = p)
  if (f$value < wl_ho$nll[i]) {
    wl_ho[i, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")] <- f$par
    wl_ho$nll[i] <- f$value
    wl_ho$conv[i] <- f$convergence
  }
}
saveRDS(wl_ho, "output/orl_wl_fit_heldout.rds")








