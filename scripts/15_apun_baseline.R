# Baseline for the A_pun finding: simulate from ORL (no outcome-dependent perseveration),
# fit ORL & the variant, and see how much the shared parameters reorder by chance.
source("scripts/10_orl_variant.R")
orl <- readRDS("output/orl_fit_horstmann.rds")
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2], net_h[choice_h == 3], net_h[choice_h == 4])

true_orl <- as.matrix(orl[, c("Arew", "Apun", "Kprime", "betaF", "betaP")])
rec_orl <- matrix(NA, nrow = nrow(true_orl), ncol = 5)
rec_wl <- matrix(NA, nrow = nrow(true_orl), ncol = 6)
colnames(rec_orl) <- c("Arew", "Apun", "Kprime", "betaF", "betaP")
colnames(rec_wl) <- c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")

set.seed(2026)
for (k in 1:nrow(true_orl)) {
  sim <- sim_orl(true_orl[k, ], deck_outcomes)
  f_orl <- fit_orl(sim$choice, sim$outcome)
  f_wl <- fit_orl_wl(sim$choice, sim$outcome, start_par = c(f_orl$par, f_orl$par[5]))
  rec_orl[k, ] <- f_orl$par
  rec_wl[k, ] <- f_wl$par
}

saveRDS(list(true = true_orl, orl = rec_orl, wl = rec_wl), "output/apun_baseline.rds")

# Baseline vs. real data. Same-person r between ORL & the variant, and the mechanism check.
# Results (Oct 30, 2026): A_pun r = .91 baseline vs. .59 real, so the reordering is beyond chance;
# mechanism check .23 baseline vs. .26 real, so that check doesn't test the mechanism.
p <- c("Arew", "Apun", "Kprime", "betaF")
wl <- readRDS("output/orl_wl_fit_horstmann.rds")
rbind(baseline = round(diag(cor(rec_orl[, p], rec_wl[, p])), 2),
      real     = round(diag(cor(orl[, p], wl[, p])), 2))
c(baseline = cor(rec_wl[, "Apun"] - rec_orl[, "Apun"], rec_wl[, "betaPl"]),
  real     = cor(wl$Apun - orl$Apun, wl$betaPl))






