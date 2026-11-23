# Exploratory: scores the hybrid model (scripts 50 and 50b) against ORL + alpha, the other models and the mixes.
# Uses the hybrid's cross-validation fits (script 50b) and the eight models' held-out predictions (script 46's
# stacking_<data>.rds; for Horstmann, run "Rscript scripts/46_prereg_stacking.R horstmann" first).
# Scores are held-out log-likelihoods summed per participant; a difference is summed over participants, with
# SE = sd of the per-participant differences x sqrt(n), as in scripts 35b and 48.
#   1. What each part adds: each version against ORL + alpha.
#   2. The hybrid against each of the eight models and the shared mix of all eight.
#   3. Mixes (shared weights fit on the other nine folds, as script 46), each against all eight: all eight plus the
#      hybrid; script 52's five_cheapest (ORL + alpha, ORL + split + alpha, VSE-split, EEF, WSLS); and those five with
#      the hybrid in place of ORL + alpha (can the hybrid stand in for it in a smaller set?).
#   4. The fitted starting values and decay rates.
# Usage (from the project folder): Rscript scripts/50c_hybrid_score.R <data>   (horstmann or replication)
library(loo)
data_name <- commandArgs(trailingOnly = TRUE)[1]
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
d <- load_set(data_name)
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
rows <- s$rows
n <- max(rows$subj)
stopifnot(n == nrow(d$choice), nrow(rows) == n * 99)
row_of <- function(i, trials) (i - 1) * 99 + trials - 1

# Held-out predictions of each version, in the rows of script 46
versions <- names(hybrid_models)
lpd_h <- matrix(NA, nrow(rows), length(versions), dimnames = list(NULL, versions))
fits <- list()
check <- 0
for (v in versions) {
  m <- hybrid_models[[v]]
  fits[[v]] <- lapply(1:10, function(f) {
    files <- Sys.glob(sprintf("output/hybrid/cv/%s/%s_fold%02d_part*.rds", data_name, v, f))
    files <- files[!grepl("_quick", files)]
    if (!length(files)) stop("no fits for ", v, " fold ", f)
    fit <- do.call(rbind, lapply(files, function(x) readRDS(x)$fit))
    fit <- fit[order(fit$subj), ]
    stopifnot(identical(as.integer(fit$subj), 1:n))
    fit
  })
  for (f in 1:10) for (i in 1:n) {
    fit <- fits[[v]][[f]]
    lp <- m$lp(unlist(fit[i, m$names]), d$choice[i, ], d$win[i, ], d$loss[i, ])
    k <- row_of(i, test_trials(f))
    lpd_h[k, v] <- lp[rows$trial[k]]
    check <- max(check, abs(sum(lp[rows$trial[k]]) - fit$heldout_ll[i]))
  }
}
check                                                   # refit predictions match the saved totals (should be about 0)
stopifnot(!anyNA(lpd_h))
# Training fits never worse than ORL + alpha's for the net-bending versions (should be about 0 or below)
sapply(setdiff(versions, "hybrid_sep"), function(v) max(sapply(fits[[v]], function(x) max(x$train_nll - x$oa_train_nll))))
# How often the ORL + alpha start gave the best fit (start 21), per version
sapply(versions, function(v) mean(unlist(lapply(fits[[v]], function(x) x$best_start == 21))))

per <- rowsum(cbind(s$lpd, shared_mix = s$cv_group, lpd_h), rows$subj)
compare <- function(a, b) { x <- per[, a] - per[, b]; c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n)), better_for = sum(x > 0)) }

# 1. What each part adds
parts <- round(t(sapply(versions, function(v) compare(v, "orl_alpha"))), 1)
parts
# 2. The hybrid against every model and the shared mix
against <- round(t(sapply(c(model_names, "shared_mix"), function(m) compare("hybrid", m))), 1)
against

# 3. Mixes
mix_score <- function(lp_rows, w) {
  ok <- w > 0
  a <- sweep(lp_rows[, ok, drop = FALSE], 2, log(w[ok]), "+")
  mx <- apply(a, 1, max)
  mx + log(rowSums(exp(a - mx)))
}
shared_cv <- function(L) {
  x <- rep(NA, nrow(L))
  for (f in 1:10) x[rows$fold == f] <- mix_score(L[rows$fold == f, , drop = FALSE], as.numeric(stacking_weights(L[rows$fold != f, , drop = FALSE])))
  x
}
five <- c("orl_alpha", "variant_alpha", "vse_split", "eef", "wsls")
five_h <- c("hybrid", "variant_alpha", "vse_split", "eef", "wsls")
L_all <- cbind(s$lpd, hybrid = lpd_h[, "hybrid"])
mixes <- list(eight_plus_hybrid = shared_cv(L_all), five = shared_cv(L_all[, five]), five_hybrid = shared_cv(L_all[, five_h]))
per <- cbind(per, rowsum(do.call(cbind, mixes), rows$subj))
mix_results <- round(rbind(t(sapply(names(mixes), function(m) compare(m, "shared_mix"))),
                           five_hybrid_vs_five = compare("five_hybrid", "five")), 1)
mix_results
w_nine <- round(setNames(as.numeric(stacking_weights(L_all)), colnames(L_all)), 3)
w_nine                                                  # weights of the eight plus the hybrid, all held-out choices
w_five_h <- round(setNames(as.numeric(stacking_weights(L_all[, five_h])), five_h), 3)
w_five_h

# 4. Fitted starting values and decay rates (each person's average over the ten folds)
par_avg <- function(v, p) rowMeans(sapply(fits[[v]], function(x) x[[p]]))
params <- rbind(Q0_start = summary(par_avg("orl_alpha_start", "Q0")), D_decay = summary(par_avg("orl_alpha_decay", "D")),
                Q0_hybrid = summary(par_avg("hybrid", "Q0")), D_hybrid = summary(par_avg("hybrid", "D")))
round(params, 3)

dir.create("output/hybrid", showWarnings = FALSE)
saveRDS(list(data = data_name, lpd_h = lpd_h, per = per, parts = parts, against = against, mix_results = mix_results,
             w_nine = w_nine, w_five_h = w_five_h, params = params),
        sprintf("output/hybrid/score_%s.rds", data_name))
