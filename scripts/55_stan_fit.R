# Hierarchical Bayesian fits of the eight models (stan/hier/, written by scripts/55_make_stan.py), with cmdstanr.
# The analysis registered in the addendum to U86SY (predictions/hierarchical.md).
#   full:  every trial of every participant; saves group means, each person's posterior means and diagnostics
#          (the log-likelihood draws go to a separate *_draws.rds, which git ignores)
#   fold:  ten-fold cross-validation within each person, the folds of scripts 35 and 45 (trials f, f + 10, ..., f + 90
#          held out; trial 1 never). Each held-out trial is scored by its posterior predictive probability, the mean
#          over draws of p(choice), so these scores are comparable with the maximum-likelihood ones.
# Check the models first: Rscript scripts/55a_stan_check.R (likelihood identical to R) and scripts/55b_stan_recovery.R.
# Usage (from the project folder):
#   Rscript scripts/55_stan_fit.R <data> <model> <fold> [quick] [long] [force]
#     data: horstmann or replication (or simulated, for testing)
#     model: orl, orl_split, orl_alpha, orl_split_alpha, nova, vse, vse_split, eef
#     fold: full, or 1-10
#     quick: the first 20 participants and short chains, for testing (saved as *_quick.rds, which git ignores)
#     long: the rerun of the convergence rule (2,000 warm-up iterations, adapt_delta 0.99), saved as *_long.rds. It
#           runs only if the first fit exists and did not pass (scripts/55_stan_common.R).
#     force: fit again even if the fit is already saved
# A fit that is already saved is skipped, so the same command can be repeated (scripts/55_queue.sh, kodiak/stan_cv.sh).
# Runs 4 chains in parallel (4 cores). A fit to Horstmann (N = 162) takes about an hour on a laptop, one to the
# replication sample (N = 440) about three.
source("scripts/55_stan_common.R")
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
model <- args[2]
fold <- if (args[3] == "full") NA else as.integer(args[3])
quick <- "quick" %in% args
long <- "long" %in% args
force <- "force" %in% args
stopifnot(model %in% stan_models, is.na(fold) || fold %in% 1:10)

out_file <- stan_file(data_name, model, fold, long, quick)
if (file.exists(out_file) && !force) {
  cat("already done:", out_file, "\n")
  quit(save = "no")
}
if (long) {
  first <- stan_file(data_name, model, fold, FALSE, quick)
  if (!file.exists(first)) {
    cat("no rerun yet: the first fit is not done (", first, ")\n")
    quit(save = "no")
  }
  if (stan_pass(readRDS(first)$diagnostics)) {
    cat("no rerun needed: the first fit passed (", first, ")\n")
    quit(save = "no")
  }
}

library(cmdstanr)
library(posterior)
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

test <- if (is.na(fold)) 1L else stan_test_trials(fold)  # a full fit scores trial 1 only, as a placeholder
train <- matrix(1L, N, 100)
if (!is.na(fold)) train[, test] <- 0L
stan_data <- list(N = N, T = 100, Tsubj = rep(100L, N),
                  choice = matrix(as.integer(d$choice[who, ]), N, 100), win = d$win[who, ], loss = d$loss[who, ],
                  train = train, n_test = length(test), test_t = matrix(as.integer(test), N, length(test), byrow = TRUE))

settings <- stan_settings(long, quick)
mod <- cmdstan_model(sprintf("stan/hier/igt_%s.stan", model))
fit <- mod$sample(data = stan_data, chains = settings$chains, parallel_chains = settings$chains,
                  iter_warmup = settings$iter_warmup, iter_sampling = settings$iter_sampling,
                  adapt_delta = settings$adapt_delta, seed = 2055 + if (is.na(fold)) 0 else fold,
                  refresh = 0, show_messages = FALSE)

# Diagnostics: divergent transitions, the largest R-hat over all group and person parameters and over the group
# parameters alone, the smallest bulk effective sample size, and each person's largest R-hat
diag <- fit$diagnostic_summary(quiet = TRUE)
summ <- fit$summary(variables = c(paste0("mu_", par_names), "sigma", par_names), "mean", "rhat", "ess_bulk")
is_group <- summ$variable %in% paste0("mu_", par_names) | startsWith(summ$variable, "sigma[")
minutes <- tryCatch(fit$time()$total / 60, error = function(e) NA)
diagnostics <- c(divergent = sum(diag$num_divergent), max_rhat = max(summ$rhat, na.rm = TRUE),
                 max_rhat_group = max(summ$rhat[is_group], na.rm = TRUE), min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
                 draws = settings$chains * settings$iter_sampling, iter_warmup = settings$iter_warmup,
                 adapt_delta = settings$adapt_delta, minutes = minutes)
print(round(diagnostics, 3))
cat(if (stan_pass(diagnostics)) "passed" else "did not pass", "the convergence rule\n")

# Each person's posterior means and largest R-hat
person <- data.frame(subj = who)
rhat_person <- matrix(NA, N, length(par_names))
for (k in seq_along(par_names)) {
  rows <- match(sprintf("%s[%d]", par_names[k], seq_len(N)), summ$variable)
  person[[par_names[k]]] <- summ$mean[rows]
  rhat_person[, k] <- summ$rhat[rows]
}
person$max_rhat <- apply(rhat_person, 1, max, na.rm = TRUE)
group <- fit$summary(variables = paste0("mu_", par_names), "mean", ~quantile(.x, c(.025, .975)))
info <- list(data = data_name, model = model, fold = if (is.na(fold)) "full" else fold, long = long, quick = quick,
             settings = settings, finished = format(Sys.time()), cmdstan = cmdstan_version())

if (is.na(fold)) {
  save_atomic(list(group = group, person = person, diagnostics = diagnostics, info = info), out_file)
  save_atomic(list(log_lik_draws = fit$draws("log_lik", format = "draws_matrix"), info = info),
              sub("\\.rds$", "_draws.rds", out_file))
} else {
  # Posterior predictive log-probability of each held-out trial: log of the mean over draws of p(choice)
  lt <- fit$draws("lp_test", format = "draws_matrix")
  log_mean_exp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
  heldout_lp <- matrix(NA, N, length(test), dimnames = list(NULL, test))
  for (i in seq_len(N)) for (k in seq_along(test)) heldout_lp[i, k] <- log_mean_exp(lt[, sprintf("lp_test[%d,%d]", i, k)])
  person$heldout_ll <- rowSums(heldout_lp)
  person$conv <- as.integer(!stan_pass(diagnostics))
  save_atomic(list(fit = person, heldout_lp = heldout_lp, group = group, diagnostics = diagnostics, info = info), out_file)
}
cat("saved", out_file, "\n")
group
