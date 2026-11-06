# Repeated model and parameter recovery for ORL and the variant.
# One replicate does what scripts 12 and 13 did once: simulate one dataset per participant
# from each model (using that participant's fitted parameters) and fit both models, with
# the same procedure as the real fits. Ten replicates (different seeds) show how much the
# recovery results vary from one simulated dataset to the next.
# Usage (from the project folder):  Rscript scripts/26_recovery_reps.R <replicate> [quick]
#   quick: first 5 participants only, for testing
args <- commandArgs(trailingOnly = TRUE)
rep_id <- as.integer(args[1])
quick <- length(args) >= 2 && args[2] == "quick"

source("scripts/10_orl_variant.R")   # data, ORL and variant functions

deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2],
                      net_h[choice_h == 3], net_h[choice_h == 4])

orl <- readRDS("output/orl_fit_horstmann.rds")
wl <- readRDS("output/orl_wl_fit_horstmann.rds")
n <- if (quick) 5 else nrow(choice_h)
true_par <- list(
  orl = as.matrix(orl[1:n, c("Arew", "Apun", "Kprime", "betaF", "betaP")]),
  wl = as.matrix(wl[1:n, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")])
)
sims <- list(orl = sim_orl, wl = sim_orl_wl)
recovered <- list(orl = true_par$orl * NA, wl = true_par$wl * NA)   # same shape, filled below
scores <- data.frame()

set.seed(2026 + rep_id)
for (gen in c("orl", "wl")) {
  for (i in 1:n) {
    sim <- sims[[gen]](true_par[[gen]][i, ], deck_outcomes)
    f_orl <- fit_orl(sim$choice, sim$outcome)
    f_wl <- fit_orl_wl(sim$choice, sim$outcome, start_par = c(f_orl$par, f_orl$par[5]))
    nll <- c(orl = f_orl$value, wl = f_wl$value)
    k <- c(orl = 5, wl = 6)
    scores <- rbind(scores, data.frame(gen = gen, subj = i, t(c(aic = 2 * nll + 2 * k, bic = 2 * nll + k * log(100)))))
    # parameter recovery: the generating model's own fit to its simulated data
    if (gen == "orl") recovered$orl[i, ] <- f_orl$par else recovered$wl[i, ] <- f_wl$par
  }
  cat("replicate", rep_id, "finished generating model:", gen, "\n")
}

dir.create("output/recovery_reps", showWarnings = FALSE)
saveRDS(list(rep = rep_id, scores = scores, true = true_par, recovered = recovered),
        sprintf("output/recovery_reps/rep_%02d%s.rds", rep_id, if (quick) "_quick" else ""))
