# Hierarchical Bayesian fits of the eight models (stan/hier/, written by scripts/55_make_stan.py), with cmdstanr.
#   full:  every trial of every participant; saves group means, each person's posterior means and diagnostics
#   fold:  ten-fold cross-validation within each person, the folds of scripts 35 and 45 (trials f, f + 10, ..., f + 90
#          held out; trial 1 never). Each held-out trial is scored by its posterior predictive probability, the mean
#          over draws of p(choice), so these scores are comparable with the maximum-likelihood ones.
# Check the models first: Rscript scripts/55a_stan_check.R (likelihood identical to R) and scripts/55b_stan_recovery.R.
# Usage (from the project folder):
#   Rscript scripts/55_stan_fit.R <data> <model> <fold> [quick]
#     data: horstmann or replication
#     model: orl, orl_split, orl_alpha, orl_split_alpha, nova, vse, vse_split, eef
#     fold: full, or 1-10
#     quick: the first 20 participants and short chains, for testing
# Runs 4 chains in parallel (4 cores). A full fit of one model to 162 people takes roughly an hour on a laptop; the
# cross-validation needs ten fits per model, so it belongs on a cluster.
library(cmdstanr)
library(posterior)
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
model <- args[2]
fold <- if (args[3] == "full") NA else as.integer(args[3])
quick <- "quick" %in% args

source("scripts/43_prereg_data.R")
d <- load_set(data_name)
who <- if (quick) 1:20 else seq_len(nrow(d$choice))
N <- length(who)

par_names <- list(orl = c("Arew", "Apun", "Kp", "betaF", "betaP"),
                  orl_split = c("Arew", "Apun", "Kp", "betaF", "betaPw", "betaPl"),
                  orl_alpha = c("Arew", "Apun", "Kp", "betaF", "betaP", "alpha"),
                  orl_split_alpha = c("Arew", "Apun", "Kp", "betaF", "betaPw", "betaPl", "alpha"),
                  nova = c("Arew", "Apun", "Kp", "betaF", "betaP", "alpha", "Q0", "D"),
                  vse = c("theta", "Delta", "alpha", "phi", "beta"),
                  vse_split = c("theta", "Delta", "alpha", "phi", "beta", "psi"),
                  eef = c("theta", "lambda", "phi", "beta"))[[model]]
stopifnot(!is.null(par_names))

test <- if (is.na(fold)) 1L else test_trials(fold)     # a full fit scores trial 1 only, as a placeholder
train <- matrix(1L, N, 100)
if (!is.na(fold)) train[, test] <- 0L
stan_data <- list(N = N, T = 100, Tsubj = rep(100L, N),
                  choice = matrix(as.integer(d$choice[who, ]), N, 100), win = d$win[who, ], loss = d$loss[who, ],
                  train = train, n_test = length(test), test_t = matrix(as.integer(test), N, length(test), byrow = TRUE))

mod <- cmdstan_model(sprintf("stan/hier/igt_%s.stan", model))
fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4,
                  iter_warmup = if (quick) 300 else 1000, iter_sampling = if (quick) 300 else 1000,
                  adapt_delta = 0.95, seed = 2055 + if (is.na(fold)) 0 else fold, refresh = 0, show_messages = FALSE)

# Diagnostics: divergent transitions, and the largest R-hat over group and person parameters
diag <- fit$diagnostic_summary(quiet = TRUE)
summ <- fit$summary(variables = c(paste0("mu_", par_names), "sigma", par_names), "mean", "rhat", "ess_bulk")
diagnostics <- c(divergent = sum(diag$num_divergent), max_rhat = max(summ$rhat, na.rm = TRUE),
                 min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE))
print(round(diagnostics, 3))

# Each person's posterior means
person <- data.frame(subj = who)
for (p in par_names) person[[p]] <- summ$mean[match(sprintf("%s[%d]", p, seq_len(N)), summ$variable)]
group <- fit$summary(variables = paste0("mu_", par_names), "mean", ~quantile(.x, c(.025, .975)))

out_dir <- if (is.na(fold)) file.path("output/stan", data_name) else file.path("output/stan/cv", data_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tag <- if (quick) "_quick" else ""
if (is.na(fold)) {
  ll <- fit$draws("log_lik", format = "draws_matrix")
  saveRDS(list(group = group, person = person, diagnostics = diagnostics, log_lik_draws = ll),
          sprintf("%s/%s_full%s.rds", out_dir, model, tag))
} else {
  # Posterior predictive log-probability of each held-out trial: log of the mean over draws of p(choice)
  lt <- fit$draws("lp_test", format = "draws_matrix")
  log_mean_exp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
  heldout_lp <- matrix(NA, N, length(test), dimnames = list(NULL, test))
  for (i in seq_len(N)) for (k in seq_along(test)) heldout_lp[i, k] <- log_mean_exp(lt[, sprintf("lp_test[%d,%d]", i, k)])
  person$heldout_ll <- rowSums(heldout_lp)
  person$conv <- as.integer(diagnostics["divergent"] > 0 || diagnostics["max_rhat"] > 1.05)
  saveRDS(list(fit = person, heldout_lp = heldout_lp, group = group, diagnostics = diagnostics),
          sprintf("%s/%s_fold%02d%s.rds", out_dir, model, fold, tag))
}
group
