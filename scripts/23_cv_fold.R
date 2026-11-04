# Cross-validation of hierarchical ORL vs. the variant, by held-out TRIALS.
# Usage (from the project folder):  Rscript scripts/23_cv_fold.R <model> <fold> [quick]
#   model: "orl" or "wl"
#   fold:  0 = fit all trials (per-trial PSIS-LOO); 1-10 = hold out every 10th trial, starting at trial <fold>
#   quick: optional; 10 people & short chains, for testing
# Held-out trials still update the model's values (their outcomes were observed);
# only their choices are left out of the fit, then predicted.
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
fold <- as.integer(args[2])
quick <- length(args) >= 3 && args[3] == "quick"

library(cmdstanr)
source("scripts/01_load_data.R")

n <- if (quick) 10 else nrow(choice_h)
n_trials <- ncol(choice_h)
train <- matrix(1L, n, n_trials)
if (fold > 0) train[, seq(fold, n_trials, by = 10)] <- 0L   # hold out trials fold, fold+10, ...

stan_data <- list(
  N = n, T = n_trials, Tsubj = rep(n_trials, n),
  choice = matrix(as.integer(choice_h[1:n, ]), n, n_trials),
  outcome = net_h[1:n, ] / 100,
  sign_out = sign(net_h[1:n, ]),
  train = train
)

stan_file <- if (model == "orl") "stan/igt_orl_cv.stan" else "stan/igt_orl_wl_cv.stan"
mod <- cmdstan_model(stan_file)
fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4,
                  iter_warmup = if (quick) 200 else 1000,
                  iter_sampling = if (quick) 200 else 1000,
                  adapt_delta = 0.95, seed = 2026, refresh = 500)

dir.create("output/cv", showWarnings = FALSE)
out <- list(model = model, fold = fold, quick = quick,
            diagnostics = fit$diagnostic_summary(),
            max_rhat = max(fit$summary(variables = c("mu_pr", "sigma"))$rhat))

ll <- fit$draws(variables = "log_lik", format = "draws_matrix")   # one column per (person, trial)
if (fold == 0) {
  out$loo <- loo::loo(ll[, !grepl(",1\\]$", colnames(ll))])       # per-trial LOO, with Pareto k
} else {
  test <- which(train == 0, arr.ind = TRUE)                       # (person, trial) of held-out choices
  cols <- paste0("log_lik[", test[, 1], ",", test[, 2], "]")
  out$test <- data.frame(subj = test[, 1], trial = test[, 2],
                         elpd = apply(ll[, cols], 2, matrixStats::logSumExp) - log(nrow(ll)))
}
saveRDS(out, sprintf("output/cv/cv_%s_fold%02d%s.rds", model, fold, if (quick) "_quick" else ""))






