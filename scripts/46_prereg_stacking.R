# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), step 3: mixes (as script 36) and the hedged
# single models (H1b, as script 42).
# Uses the cross-validation fits from script 45 (no refitting): each fold's fitted parameters are run again to get the
# probability every model gave every held-out choice (participants x 99 choices x 8 models); WSLS is scored as in
# scripts 35b and 36. Then, with weights never fit on the fold being scored:
#   shared mix:    one set of weights for everyone, fit on the other nine folds of all participants
#   personal mix:  each person's own weights, fit on their other nine folds
#   shrinkage mix: lambda x own weights + (1 - lambda) x shared weights, lambda in 0, 0.1, ..., 1; for each fold, the
#                  lambda that scored best on the other nine folds
#   hedged model:  p = (1 - e) x p_model + e / 4, e fit like the mix's weights (stacking the model with a random guess)
# If stacking weights cannot be fit for one person and fold, that person and fold use the shared weights (registration,
# "If something goes wrong"); the number of times this happens is printed.
# Usage (from the project folder): Rscript scripts/46_prereg_stacking.R <data>
#   data: replication, simulated, replication342 (robustness check 5), replication438 (without the repeated
#         records; exploratory), or horstmann (reproduces script 36, as a check)
# Run it with the project's renv (loo 2.10.1, as for script 36): personal weights are fit to only about 90 choices,
# and other loo versions or machines can stop the optimizer at slightly different weights (a check against script 36
# with loo 2.6.0 on Linux changed the shrinkage total by about 1). The versions used are printed and saved.
library(loo)
data_name <- commandArgs(trailingOnly = TRUE)[1]
versions <- c(R = R.version.string, loo = as.character(packageVersion("loo")), platform = R.version$platform)
versions
source("scripts/43_prereg_data.R")
source("scripts/43b_prereg_models.R")
src <- cv_source(data_name)
d <- load_set(src)
n_all <- nrow(d$choice)
M <- length(model_names)

# Each fold's fits for one model, all participants in order. Script 45's files if there are any (in parts);
# for the Horstmann check, script 35's.
read_fits <- function(m, f) {
  files <- Sys.glob(sprintf("output/prereg/cv/%s/%s_fold%02d_part*.rds", src, m, f))
  files <- files[!grepl("_quick", files)]
  if (length(files)) {
    fit <- do.call(rbind, lapply(files, function(x) readRDS(x)$fit))
  } else if (src == "horstmann") {
    fit <- readRDS(sprintf("output/cv_mle/%s_fold%02d.rds", m, f))
  } else stop("no cross-validation fits for ", m, " fold ", f)
  fit <- fit[order(fit$subj), ]
  stopifnot(identical(as.integer(fit$subj), 1:n_all))
  fit
}

# One row per held-out choice (participant x trial, trials 2-100), one column per model, as in script 36
rows <- expand.grid(trial = 2:100, subj = 1:n_all)
rows$fold <- (rows$trial - 1) %% 10 + 1
row_of <- function(i, trials) (i - 1) * 99 + trials - 1   # row of participant i's trials (rows go participant by participant)
lpd <- matrix(NA, nrow(rows), M, dimnames = list(NULL, model_names))
check <- 0
for (m in names(models)) {
  for (f in 1:10) {
    fit <- read_fits(m, f)
    for (i in 1:n_all) {
      lp <- models[[m]]$lp(unlist(fit[i, models[[m]]$names]), d$choice[i, ], d$win[i, ], d$loss[i, ])
      k <- row_of(i, test_trials(f))
      lpd[k, m] <- lp[rows$trial[k]]
      check <- max(check, abs(sum(lp[rows$trial[k]]) - fit$heldout_ll[i]))
    }
  }
}
check                                                   # rerun fits match the saved held-out totals (should be about 0)
W <- wsls_heldout(d$choice, d$win, d$loss)
lpd[, "wsls"] <- W[cbind(rows$subj, rows$trial - 1)]
stopifnot(!anyNA(lpd), all(rows$trial[row_of(3, test_trials(4))] == test_trials(4)), all(rows$subj[row_of(3, test_trials(4))] == 3))

# The participants this analysis uses (all, except for robustness check 5)
keep <- keep_of(data_name, d$study)
lpd <- lpd[keep[rows$subj], ]
rows <- rows[keep[rows$subj], ]
rows$subj <- match(rows$subj, which(keep))
n <- sum(keep)
study <- d$study[keep]
scheme <- d$scheme[keep]

# Score of a mix: log of the weighted average probability (computed stably), as in script 36
mix_score <- function(lp_rows, w) {
  ok <- w > 0
  a <- sweep(lp_rows[, ok, drop = FALSE], 2, log(w[ok]), "+")
  mx <- apply(a, 1, max)
  mx + log(rowSums(exp(a - mx)))
}

# Shared weights: from all held-out choices (description), and without each fold (for scoring)
w_group <- setNames(as.numeric(stacking_weights(lpd)), model_names)
w_grp_f <- t(sapply(1:10, function(f) as.numeric(stacking_weights(lpd[rows$fold != f, ]))))
cv_group <- rep(NA, nrow(rows))
for (f in 1:10) cv_group[rows$fold == f] <- mix_score(lpd[rows$fold == f, ], w_grp_f[f, ])

# Each person's own weights without each fold; the shared weights when they cannot be fit
fallbacks <- 0
w_own <- array(NA, c(n, 10, M))
own_rows <- split(seq_len(nrow(rows)), rows$subj)       # each participant's rows
for (i in 1:n) for (f in 1:10) {
  k <- own_rows[[i]]
  w <- tryCatch(as.numeric(stacking_weights(lpd[k[rows$fold[k] != f], ])), error = function(e) NULL)
  if (is.null(w) || anyNA(w)) { w <- w_grp_f[f, ]; fallbacks <- fallbacks + 1 }
  w_own[i, f, ] <- w
}
fallbacks                                               # person-folds that used the shared weights

# Personal and shrinkage mixes
lambdas <- seq(0, 1, by = 0.1)
score <- matrix(NA, nrow(rows), length(lambdas))
for (i in 1:n) for (f in 1:10) {
  te <- own_rows[[i]][rows$fold[own_rows[[i]]] == f]
  for (k in seq_along(lambdas)) score[te, k] <- mix_score(lpd[te, , drop = FALSE], lambdas[k] * w_own[i, f, ] + (1 - lambdas[k]) * w_grp_f[f, ])
}
cv_person <- score[, length(lambdas)]
max(abs(score[, 1] - cv_group))                         # lambda = 0 is the shared mix (should print 0)
fold_tot <- rowsum(score, rows$fold)
best_k <- sapply(1:10, function(f) which.max(colSums(fold_tot[-f, ])))
lambda_chosen <- lambdas[best_k]
lambda_chosen                                           # lambda chosen for each fold
cv_shrink <- score[cbind(1:nrow(rows), best_k[rows$fold])]

# Each person's own weights from all their held-out choices (description: personal weight profiles)
w_person <- t(sapply(1:n, function(i) {
  w <- tryCatch(as.numeric(stacking_weights(lpd[own_rows[[i]], ])), error = function(e) rep(NA, M))
  if (anyNA(w)) w_group else w
}))
colnames(w_person) <- model_names

# Hedged single models (H1b and robustness check 3), as in script 42
hedged <- sapply(model_names, function(m) {
  x <- rep(NA, nrow(rows))
  for (f in 1:10) {
    te <- rows$fold == f
    w <- as.numeric(stacking_weights(cbind(lpd[!te, m], log(0.25))))
    x[te] <- log(w[1] * exp(lpd[te, m]) + w[2] * 0.25)
  }
  x
})

# Per-participant totals, for scripts 47 and 48
per <- rowsum(cbind(lpd, shared_mix = cv_group, personal_mix = cv_person, shrink_mix = cv_shrink,
                    setNames(as.data.frame(hedged), paste0(model_names, "_hedged"))), rows$subj)
dir.create("output/prereg", showWarnings = FALSE)
saveRDS(list(data = data_name, versions = versions, lpd = lpd, rows = rows, study = study, scheme = scheme, w_group = w_group,
             w_grp_f = w_grp_f, w_own = w_own, w_person = w_person, fallbacks = fallbacks, lambdas = lambdas,
             lambda_chosen = lambda_chosen, cv_group = cv_group, cv_person = cv_person, cv_shrink = cv_shrink,
             hedged = hedged, per = per),
        sprintf("output/prereg/stacking_%s.rds", data_name))
