# Shared by the hierarchical Stan scripts (55_stan_fit.R, 55c, 55d, 55e): where each fit is saved, the order the fits
# are run in, and the convergence rule of the addendum to U86SY (predictions/hierarchical.md). Loads no data.

stan_models <- c("orl", "orl_split", "orl_alpha", "orl_split_alpha", "nova", "vse", "vse_split", "eef")

# Held-out trials of fold f, as test_trials() in scripts/43b_prereg_models.R (trial 1 is never held out). Kept here so
# the fits need none of the maximum-likelihood code; script 55c checks that every fit used the same trials.
stan_test_trials <- function(f) { tt <- seq(f, 100, by = 10); tt[tt != 1] }

# Sampler settings: the registered ones, and the one rerun allowed by the convergence rule
stan_settings <- function(long = FALSE, quick = FALSE) {
  s <- list(chains = 4, iter_warmup = if (long) 2000 else 1000, iter_sampling = 1000,
            adapt_delta = if (long) 0.99 else 0.95)
  if (quick) { s$iter_warmup <- s$iter_warmup / 10 * 3; s$iter_sampling <- 300 }   # quick tests only
  s
}

# Where a fit is saved. fold: NA (or "full") for the fit to all trials, 1-10 for a cross-validation fold.
stan_file <- function(data, model, fold, long = FALSE, quick = FALSE) {
  tag <- paste0(if (long) "_long" else "", if (quick) "_quick" else "")
  if (is.na(fold) || identical(fold, "full")) sprintf("output/stan/%s/%s_full%s.rds", data, model, tag)
  else sprintf("output/stan/cv/%s/%s_fold%02d%s.rds", data, model, as.integer(fold), tag)
}

# The convergence rule. A fit passes with no divergent transitions and no R-hat above 1.05 (group and person
# parameters). A fit that does not pass is rerun once with stan_settings(long = TRUE). After that, a model whose fits
# still include one with a group-level R-hat above 1.05 or more than 1% divergent transitions counts as not tested.
stan_pass <- function(diagnostics) diagnostics[["divergent"]] == 0 && diagnostics[["max_rhat"]] <= 1.05
stan_unusable <- function(diagnostics)
  diagnostics[["max_rhat_group"]] > 1.05 || diagnostics[["divergent"]] > 0.01 * diagnostics[["draws"]]

# The fit to use: the rerun when there is one, otherwise the first fit (NA if neither exists yet)
stan_final_file <- function(data, model, fold, quick = FALSE) {
  for (long in c(TRUE, FALSE)) { f <- stan_file(data, model, fold, long, quick); if (file.exists(f)) return(f) }
  NA_character_
}

# Comparing two columns of script 55c's per-person scores (sc = its saved output). Differences are summed over
# participants; SE = sd of the per-participant differences x sqrt(n), as every other test. converged_only drops the
# people with an R-hat above 1.05 in any fit of the Stan models the two columns use (reported; it does not decide).
stan_learning_models <- c("orl", "variant", "orl_alpha", "variant_alpha", "vse", "vse_split", "eef")  # 43b's names(models)
stan_uses <- function(col) if (startsWith(col, "hier.")) sub("hier.", "", col, fixed = TRUE) else
  if (col == "eight_hier") stan_learning_models else character(0)
stan_compare <- function(sc, a, b, converged_only = FALSE) {
  keep <- rep(TRUE, nrow(sc$per))
  if (converged_only) for (mn in c(stan_uses(a), stan_uses(b))) keep[sc$bad_people[[mn]]] <- FALSE
  x <- sc$per[keep, a] - sc$per[keep, b]
  c(diff = sum(x), se = sd(x) * sqrt(sum(keep)), z = sum(x) / (sd(x) * sqrt(sum(keep))), better_for = sum(x > 0),
    n = sum(keep))
}

# Save so that a fit stopped halfway never leaves a file that looks finished
save_atomic <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(file, ".part")
  saveRDS(object, tmp)
  stopifnot(file.rename(tmp, file))
}
