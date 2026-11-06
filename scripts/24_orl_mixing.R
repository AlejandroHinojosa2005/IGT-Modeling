# Why does hierarchical ORL mix poorly in some cross-validation folds?
# Refits ORL on one fold with the same data, mask and sampler settings as the final
# cross-validation run (script 23), but saves R-hat for every parameter and the draws
# needed to see where the chains disagree.
# Usage (from the project folder):  Rscript scripts/24_orl_mixing.R <fold> [quick]
args <- commandArgs(trailingOnly = TRUE)
fold <- as.integer(args[1])
quick <- length(args) >= 2 && args[2] == "quick"

library(cmdstanr)
source("scripts/01_load_data.R")

n <- if (quick) 10 else nrow(choice_h)
n_trials <- ncol(choice_h)
train <- matrix(1L, n, n_trials)
if (fold > 0) train[, seq(fold, n_trials, by = 10)] <- 0L

stan_data <- list(
  N = n, T = n_trials, Tsubj = rep(n_trials, n),
  choice = matrix(as.integer(choice_h[1:n, ]), n, n_trials),
  outcome = net_h[1:n, ] / 100,
  sign_out = sign(net_h[1:n, ]),
  train = train
)

n_iter <- if (quick) 200 else 2000
mod <- cmdstan_model("stan/igt_orl_cv.stan")
fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4,
                  iter_warmup = n_iter, iter_sampling = n_iter,
                  adapt_delta = 0.95, seed = 2026, refresh = 500)

# R-hat and ESS for every parameter (K here is K', on the 0-5 scale)
summ <- as.data.frame(fit$summary(variables = c("mu_pr", "sigma", "Arew", "Apun", "K", "betaF", "betaP")))

# Group-level draws, plus which chain each draw came from (draws are stacked chain by chain)
group <- unclass(fit$draws(variables = c("mu_pr", "sigma"), format = "draws_matrix"))
chain <- rep(1:4, each = n_iter)

# Draws for the 5 participants whose parameters mixed worst
ind <- summ[!grepl("mu_pr|sigma", summ$variable), ]
ind <- ind[order(-ind$rhat), ]
worst_subj <- unique(as.integer(sub(".*\\[([0-9]+)\\]", "\\1", ind$variable)))[1:5]
worst_vars <- paste0(rep(c("Arew", "Apun", "K", "betaF", "betaP"), times = 5), "[", rep(worst_subj, each = 5), "]")
worst <- unclass(fit$draws(variables = worst_vars, format = "draws_matrix"))

saveRDS(list(fold = fold, quick = quick, diagnostics = fit$diagnostic_summary(), summary = summ,
             group = group, chain = chain, worst_subj = worst_subj, worst = worst),
        sprintf("output/orl_mixing_fold%02d%s.rds", fold, if (quick) "_quick" else ""))
