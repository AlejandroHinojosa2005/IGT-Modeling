# Hierarchical Bayesian fit of the ORL variant (Stan, via cmdstanr) to the 162 Horstmann participants.
# Same settings as the hBayesDM run (script 16): 4 chains, 1,000 warmup + 3,000 sampling.
# ~1.5-3 hours: run as a background job. Everything needed later is saved, since the fit itself isn't.
library(cmdstanr)
source("scripts/01_load_data.R")

stan_data <- list(
  N = nrow(choice_h), T = 100, Tsubj = rep(100, nrow(choice_h)),
  choice = matrix(as.integer(choice_h), nrow(choice_h), 100),
  outcome = net_h / 100,
  sign_out = sign(net_h)
)

mod <- cmdstan_model("stan/igt_orl_wl.stan")
fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 3000,
                  adapt_delta = 0.95, seed = 2026, refresh = 500)

# Posterior mean of each participant's parameters (Stan's K is on the same 0-5 scale as Kprime)
stan_names <- c("Arew", "Apun", "K", "betaF", "betaPw", "betaPl")
wl_names <- c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")
ind <- matrix(NA, nrow(choice_h), 6, dimnames = list(NULL, wl_names))
for (p in 1:6) {
  ind[, p] <- colMeans(fit$draws(variables = stan_names[p], format = "draws_matrix"))
}

# Group-level draws, for the population-level asymmetry test
group_draws <- fit$draws(variables = c("mu_Arew", "mu_Apun", "mu_K", "mu_betaF", "mu_betaPw", "mu_betaPl"),
                         format = "draws_df")

saveRDS(list(
  ind = ind,                                                          # posterior means, per participant
  group_draws = group_draws,                                          # group-level draws
  diagnostics = fit$diagnostic_summary(),                             # divergences, treedepth, E-BFMI
  summary = fit$summary(variables = c("mu_pr", "sigma", stan_names))  # Rhat & ESS for every parameter
), "output/stan_wl_horstmann.rds")





