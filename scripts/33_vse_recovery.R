# Recovery for VSE and VSE-split: simulate one dataset per participant from their fitted
# parameters (script 32), then refit both models with the same procedure as the real fits.
#   gen = vse_split: parameter recovery (prediction 6: psi recovers at r >= .70), and how often
#                    the split is picked when it is true
#   gen = vse:       how often the split is wrongly picked when the data come from plain VSE
# Usage (from the project folder):  Rscript scripts/33_vse_recovery.R <vse|vse_split> [quick]
args <- commandArgs(trailingOnly = TRUE)
gen <- args[1]
quick <- length(args) >= 2 && args[2] == "quick"
tag <- if (quick) "_quick" else ""

source("scripts/31_vse_models.R")
names_vse <- c("theta", "Delta", "alpha", "phi", "beta")
names_split <- c(names_vse, "psi")

fits <- readRDS(sprintf("output/vse/fit_%s%s.rds", gen, tag))
true_par <- as.matrix(fits[, if (gen == "vse") names_vse else names_split])
n <- nrow(true_par)
rec_vse <- matrix(NA, n, 5, dimnames = list(NULL, names_vse))
rec_split <- matrix(NA, n, 6, dimnames = list(NULL, names_split))
nll_vse <- nll_split <- rep(NA, n)

set.seed(2026)
for (i in 1:n) {
  sim <- if (gen == "vse") sim_vse(true_par[i, ], deck_draws) else sim_vse_split(true_par[i, ], deck_draws)
  f_vse <- fit_multi(vse_nll, sim$choice, sim$win, sim$loss, lower_vse, upper_vse)
  f_split <- fit_multi(vse_split_nll, sim$choice, sim$win, sim$loss, lower_vse_split, upper_vse_split,
                       start_par = c(f_vse$par, 0))
  rec_vse[i, ] <- f_vse$par
  rec_split[i, ] <- f_split$par
  nll_vse[i] <- f_vse$value
  nll_split[i] <- f_split$value
}

saveRDS(list(gen = gen, true = true_par, rec_vse = rec_vse, rec_split = rec_split,
             nll_vse = nll_vse, nll_split = nll_split),
        sprintf("output/vse/recovery_%s%s.rds", gen, tag))

# How often each criterion picks the split (one extra parameter: AIC needs > 1, BIC > log(100)/2)
c(aic = mean(nll_vse - nll_split > 1), bic = mean(nll_vse - nll_split > log(100) / 2))
# Parameter recovery, generating model's parameters
if (gen == "vse") round(diag(cor(true_par, rec_vse)), 2) else round(diag(cor(true_par, rec_split)), 2)
