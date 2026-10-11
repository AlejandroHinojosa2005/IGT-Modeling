# Model Comparison: PVL-Delta vs VPP, per participant
# AIC = 2*nll + 2*k (k = # of free parameters)
# BIC = 2*nll + k*log(n_trials) (stricter penalty per parameter)
# Lower is better for both.

results_pvl <- readRDS("output/pvl_delta_fit_horstmann.rds")
results_vpp <- readRDS("output/vpp_fit_horstmann.rds")
results_orl <- readRDS("output/orl_fit_horstmann.rds")
n_trials <- 100

comp <- data.frame(
  subj = results_pvl$subj,
  aic_pvl = 2 * results_pvl$nll + 2 * 4,
  aic_vpp = 2 * results_vpp$nll + 2 * 8,
  aic_orl = 2 * results_orl$nll + 2 * 5,
  bic_pvl = 2 * results_pvl$nll + 4 * log(n_trials),
  bic_vpp = 2 * results_vpp$nll + 8 * log(n_trials),
  bic_orl = 2 * results_orl$nll + 5 * log(n_trials)
)

# WHich model wins for each participant? (1=PVL-Delta, 2= VPP, 3=ORL)
aic <- comp[, c("aic_pvl", "aic_vpp", "aic_orl")]
bic <- comp[, c("bic_pvl", "bic_vpp", "bic_orl")]
table(apply(aic, 1, which.min))
table(apply(bic, 1, which.min))

# Total across participants (lower = better)
colSums(aic)
colSums(bic)

saveRDS(comp, "output/model_comparison_all.rds")