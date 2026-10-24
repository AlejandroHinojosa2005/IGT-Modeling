# Model Recovery: can AIC/BIC tell PVL-Delta, VPP, & ORL apart?
# For each generating model, simulate one synthetic participant per real participant
# (using that participant's fitted parameters), fit all three models, record the scorse.

source("scripts/01_load_data.R")
source("scripts/02_pvl_delta.R")
source("scripts/03_vpp.R")
source("scripts/04_orl.R")

deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2],
                      net_h[choice_h == 3], net_h[choice_h == 4])

pvl <- readRDS("output/pvl_delta_fit_horstmann.rds")
vpp <- readRDS("output/vpp_fit_horstmann.rds")
orl <- readRDS("output/orl_fit_horstmann.rds")

true_par <- list(
  pvl = as.matrix(pvl[, c("A", "alpha", "cons", "lambda")]),
  vpp = as.matrix(vpp[, c("A", "alpha", "cons", "lambda", "epP", "epN", "K", "w")]),
  orl = as.matrix(orl[, c("Arew", "Apun", "Kprime", "betaF", "betaP")])
)
sims <- list(pvl = sim_pvl_delta, vpp = sim_vpp, orl = sim_orl)

# Fit all three models to one dataset; return AIC and BIC for each
fit_all <- function(choice, outcome) {
  f_pvl <- fit_pvl_delta(choice, outcome)
  f_vpp <- fit_vpp(choice, outcome, n_starts = 50)
  f_orl <- fit_orl(choice, outcome)
  nll <- c(pvl = pvl_delta_nll(f_pvl$par, choice, outcome),
           vpp = vpp_nll(f_vpp$par, choice, outcome),
           orl = orl_nll(f_orl$par, choice, outcome))
  k <- c(pvl = 4, vpp = 8, orl = 5)
  c(aic = 2 * nll + 2 * k, bic = 2 * nll + k * log(100))
}

set.seed(2026)
gens <- c("pvl", "vpp", "orl")
results_mr <- data.frame()
for (gen in gens) {
  for (i in 1:nrow(true_par[[gen]])) {
    sim <- sims[[gen]](true_par[[gen]][i, ], deck_outcomes)
    scores <- fit_all(sim$choice, sim$outcome)
    results_mr <- rbind(results_mr, data.frame(gen = gen, subj = i, t(scores)))
  }
  cat("finished generating model:", gen, "\n")
}

saveRDS(results_mr, "output/model_recovery.rds")















