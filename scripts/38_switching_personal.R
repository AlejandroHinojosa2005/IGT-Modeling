# Exploratory: the switching model (script 37), redrawing with each person's own weights.
# In script 37, a switch redraws the model in charge using the group weights w. Here it redraws using
#   lambda * the person's own weights + (1 - lambda) * w
# The person's own weights come from stacking on their other folds' choices (as in script 36, section 5).
# No leaks: to score the trials of fold f, the whole session is run with the weights fit without fold f,
# so nothing going into a scored trial has seen it. (Using each trial's own fold weights would leak: the
# belief carries over from earlier trials, whose weights were fit with the current trial included.)
# w, delta and lambda are fit on 9/10 of the participants and scored on the other 1/10, as in script 37.
library(loo)
s <- readRDS("output/stacking.rds")
n <- max(s$rows$subj)
L <- aperm(array(s$lpd, c(99, n, ncol(s$lpd))), c(2, 1, 3))   # participants x trials (2-100) x models
dimnames(L)[[3]] <- colnames(s$lpd)
m <- dim(L)[3]

# Each person's own weights for each fold, fit without that fold (takes a few minutes)
w_own <- array(NA, c(n, 10, m))
for (i in 1:n) {
  for (f in 1:10) {
    w_own[i, f, ] <- as.numeric(stacking_weights(s$lpd[s$rows$subj == i & s$rows$fold != f, ]))
  }
}

# Every session is run 10 times, once per fold ("copies"); copy k uses one person's weights for one fold
# and only that fold's trials are scored
copy_subj <- rep(1:n, each = 10)
copy_fold <- rep(1:10, n)
W_copy <- t(sapply(seq_along(copy_subj), function(k) w_own[copy_subj[k], copy_fold[k], ]))  # copies x models
scored <- outer(copy_fold, (1:99) %% 10 + 1, "==")   # copies x trials: TRUE where the trial is in the copy's fold
# Check: mixing with only the person's own weights gives script 36's personal mixes (should print about 0)
static <- log(apply(exp(L[copy_subj, , ]) * aperm(array(W_copy, c(10 * n, m, 99)), c(1, 3, 2)), c(1, 2), sum))
max(abs(rowsum(rowSums(static * scored), copy_subj)[, 1] - rowsum(s$cv_person, s$rows$subj)[, 1]))

# Forward algorithm, as in script 37, but the redraw weights can differ by copy
hmm_ll <- function(L, W) {                              # W: copies x models, the redraw weights
  out <- matrix(NA, dim(L)[1], dim(L)[2])
  B <- W
  for (t in 1:dim(L)[2]) {
    A <- log(B) + L[, t, ]
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]
    out[, t] <- mx + log(rowSums(exp(A - mx)))
    B <- (1 - attr(W, "delta")) * exp(A - out[, t]) + attr(W, "delta") * W
  }
  out
}
redraw <- function(W_copy, w, delta, lambda) {          # lambda * own weights + (1 - lambda) * w, with delta attached
  W <- lambda * W_copy + (1 - lambda) * matrix(w, nrow(W_copy), m, byrow = TRUE)
  attr(W, "delta") <- delta
  W
}
score_copies <- function(k, w, delta, lambda) {         # each participant's total over their scored trials
  out <- hmm_ll(L[copy_subj[k], , ], redraw(W_copy[k, ], w, delta, lambda))
  rowsum(rowSums(out * scored[k, ]), copy_subj[k])[, 1]
}

fit_hmm <- function(k, lambda = NA) {
  free <- is.na(lambda)
  unpack <- function(par) {
    w <- exp(c(0, par[1:7]))
    list(w = setNames(w / sum(w), dimnames(L)[[3]]), delta = plogis(par[8]),
         lambda = if (free) plogis(par[9]) else lambda)
  }
  nll <- function(par) { u <- unpack(par); -sum(score_copies(k, u$w, u$delta, u$lambda)) }
  best <- NULL
  for (d0 in c(-3, 0)) {                                # delta starts at .05 and .5
    opt <- optim(c(rep(0, 7), d0, if (free) 0), nll, method = "BFGS", control = list(maxit = 500))
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  unpack(best$par)
}

# Check: lambda = 1 with delta = 1 is script 36's personal mix again (should print about 0)
all_k <- seq_along(copy_subj)
max(abs(score_copies(all_k, rep(1 / m, m), 1, 1) - rowsum(s$cv_person, s$rows$subj)[, 1]))

# Cross-validation by participant, with the same groups as script 37
grp <- readRDS("output/switching.rds")$grp
cv <- rep(NA, n)
fits <- NULL
for (g in 1:10) {
  fit <- fit_hmm(which(grp[copy_subj] != g))
  test <- which(grp[copy_subj] == g)
  cv[sort(unique(copy_subj[test]))] <- score_copies(test, fit$w, fit$delta, fit$lambda)
  fits <- rbind(fits, c(delta = fit$delta, lambda = fit$lambda))
}
round(fits, 3)                                          # delta and lambda chosen in each group

# Results (positive = the first predicts better)
cv <- cbind(personal_switching = cv, switching = readRDS("output/switching.rds")$cv[, "switching"],
            shrink_mix = rowsum(s$cv_shrink, s$rows$subj)[, 1], orl_alpha = rowSums(L[, , "orl_alpha"]))
round(colSums(cv), 1)
compare <- function(a, b) {
  d <- cv[, a] - cv[, b]
  round(c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0)), 1)
}
round(rbind(personal_vs_switching = compare("personal_switching", "switching"),
            personal_vs_shrink_mix = compare("personal_switching", "shrink_mix"),
            personal_vs_orl_alpha = compare("personal_switching", "orl_alpha")), 1)

saveRDS(list(cv = cv, fits = fits, w_own = w_own), "output/switching_personal.rds")

