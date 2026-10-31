# Hierarchical Bayesian fits (hBayesDM 2.0, via CmdStan) of PVL-Delta, VPP & ORL
# to the Horstmann data, as a cross-check on the maximum-likelihood fits.
# Slow: run as a background job.
library(hBayesDM)
source("scripts/01_load_data.R")

# One row per trial: subjID, choice, gain, loss
dat_h <- data.frame(
  subjID = rep(1:nrow(choice_h), each = ncol(choice_h)),
  choice = as.vector(t(choice_h)),
  gain = as.vector(t(wi_h)),
  loss = as.vector(t(lo_h))
)

# Save each model's summaries as soon as it finishes, so a later error can't lose earlier fits.
save_hb <- function(fit, name) {
  saveRDS(list(ind = fit$all_ind_pars,
               rhat = rhat(fit),
               looic = extract_ic(fit, ic = "looic")),
          paste0("output/hbayesdm_", name, "_horstmann.rds"))
}

fit_orl_hb <- igt_orl(data = dat_h, niter = 4000, nwarmup = 1000, nchain = 4, ncore = 4, seed = 2026)
save_hb(fit_orl_hb, "orl")
fit_pvl_hb <- igt_pvl_delta(data = dat_h, niter = 4000, nwarmup = 1000, nchain = 4, ncore = 4, seed = 2026)
save_hb(fit_pvl_hb, "pvl")
fit_vpp_hb <- igt_vpp(data = dat_h, niter = 4000, nwarmup = 1000, nchain = 4, ncore = 4, seed = 2026)
save_hb(fit_vpp_hb, "vpp")

# LOOIC comparison of all three (lower = better)
looic_all <- print_fit(fit_pvl_hb, fit_vpp_hb, fit_orl_hb, ic = "looic")
looic_all
saveRDS(looic_all, "output/hbayesdm_looic_horstmann.rds")

library(cmdstanr)   # Stan backend for hBayesDM 2.0 (listed so renv records it)

