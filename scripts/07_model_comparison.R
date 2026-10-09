# Model Comparison: PVL-Delta vs VPP, per participant
# AIC = 2*nll + 2*k (k = # of free parameters)
# BIC = 2*nll + k*log(n_trials) (stricter penalty per parameter)
# Lower is better for both.

results_pvl <- readRDS("output/pvl_delta_fit_horstmann.rds")
results_vpp <- readRDS("output/vpp_fit_horstmann.rds")
n_trials <- 100

comp <- data.frame(
  subj = results_pvl$subj,
  aic_pvl = 2 * results_pvl$nll + 2 * 4,
  aic_vpp = 2 * results_vpp$nll + 2 * 8,
  bic_pvl = 2 * results_pvl$nll + 4 * log(n_trials),
  bic_vpp = 2 * results_vpp$nll + 8 * log(n_trials)
)

# How many participants are better described by VPP?
sum(comp$aic_vpp < comp$aic_pvl)
sum(comp$bic_vpp < comp$bic_pvl)

# Total across all participants (lower total = better overall)
colSums(comp[, c("aic_pvl", "aic_vpp", "bic_pvl", "bic_vpp")])

saveRDS(comp, "output/model_comparison_pvl_vpp.rds")