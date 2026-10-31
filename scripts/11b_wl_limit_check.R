# Sensitivity: does the +/-10 limit cost the variant fit? Refit participants with a
# betaP at the limit using +/-30 (the variant's version of ORL's check in 05, which gained 30.6 nll).
# Each refit starts at that participant's +/-10 fit, so it can't come out worse.
source("scripts/10_orl_variant.R")
results_orl_wl <- readRDS("output/orl_wl_fit_horstmann.rds")
pinned_wl <- which(abs(results_orl_wl$betaPw) >= 9.999 | abs(results_orl_wl$betaPl) >= 9.999)
nll_wide_wl <- rep(NA, length(pinned_wl))

set.seed(2026)
for (n in 1:length(pinned_wl)) {
  i <- pinned_wl[n]
  p <- unlist(results_orl_wl[i, c("Arew","Apun","Kprime","betaF","betaPw","betaPl")])
  nll_wide_wl[n] <- fit_orl_wl(choice_h[i, ], net_h[i, ], beta_lim = 30, start_par = p)$value
}

improve_wl <- results_orl_wl$nll[pinned_wl] - nll_wide_wl  # positive = wider limit fits better
saveRDS(list(pinned = pinned_wl, improve = improve_wl), "output/orl_wl_limit_check.rds")

length(pinned_wl)         # participants at the limit
sum(improve_wl)           # total nll gain (ORL's was 30.6)
sum(improve_wl < -1e-8)   # refits that did worse (should be 0)

# Does summed BIC flip with both models at +/-30? Variant minus ORL was 41.18 at +/-10.
# Negative = BIC now favors the variant. ORL's 30.6 slightly understates its gain
# (4 of its refits were search misses), so this comparison leans slightly toward the variant.
41.18 - 2 * (sum(improve_wl) - 30.6)





