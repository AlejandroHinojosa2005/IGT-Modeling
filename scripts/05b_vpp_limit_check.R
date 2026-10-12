# check: is VPP's epsilon limit (+/- 5) consting fit?
# Refit participants whose epP or epN hit the limit, with limits of +/- 15
#Mirror of ORL betaP check in 05_fit_models.R
source("scripts/03_vpp.R")
results_vpp <- readRDS("output/vpp_fit_horstmann.rds")
pinned_vpp <- which(abs(results_vpp$epP) >= 4.999 | abs(results_vpp$epN) >= 4.999)
nll_wide_vpp <- rep(NA, length(pinned_vpp))

set.seed(2026)
for (n in 1:length(pinned_vpp)) {
  i <- pinned_vpp[n]
  nll_wide_vpp[n] <- fit_vpp(choice_h[i, ], net_h[i, ], n_starts=50,eps_lim=15)$value
}

improve_vpp <- results_vpp$nll[pinned_vpp] - nll_wide_vpp
saveRDS(list(pinned = pinned_vpp,improve = improve_vpp), "output/vpp_limit_check.rds")