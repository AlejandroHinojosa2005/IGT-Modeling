# Model Recovery: can AIC/BIC tell ORL & ORL_variant apart?
# For each generating model, simulate one synthetic participant per real participant
# (using that participant's fitted parameters), both models, record the scores.

source("scripts/10_orl_variant.R")

deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2],
                      net_h[choice_h == 3], net_h[choice_h == 4])

orl <- readRDS("output/orl_fit_horstmann.rds")
wl <- readRDS("output/orl_wl_fit_horstmann.rds")

true_par <- list(
  orl = as.matrix(orl[, c("Arew", "Apun", "Kprime", "betaF", "betaP")]),
  wl = as.matrix(wl[, c("Arew","Apun","Kprime","betaF","betaPw","betaPl")])
)
sims <- list(orl = sim_orl, wl = sim_orl_wl)

# Fit all both models to one dataset; return AIC and BIC for each
fit_all <- function(choice, outcome) {
  f_orl <- fit_orl(choice, outcome)
  f_wl <- fit_orl_wl(choice, outcome, start_par = c(f_orl$par, f_orl$par[5]))
  nll <- c(orl = orl_nll(f_orl$par, choice, outcome),
           wl = orl_wl_nll(f_wl$par, choice, outcome))
  k <- c(orl = 5, wl = 6)
  c(aic = 2 * nll + 2 * k, bic = 2 * nll + k * log(100))
}

set.seed(2026)
gens <- c("orl", "wl")
results_mr <- data.frame()
for (gen in gens) {
  for (i in 1:nrow(true_par[[gen]])) {
    sim <- sims[[gen]](true_par[[gen]][i, ], deck_outcomes)
    scores <- fit_all(sim$choice, sim$outcome)
    results_mr <- rbind(results_mr, data.frame(gen = gen, subj = i, t(scores)))
  }
  cat("finished generating model:", gen, "\n")
}

saveRDS(results_mr, "output/model_recovery_wl.rds")