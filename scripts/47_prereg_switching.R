# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), steps 4-5: switching, drifting weights and
# personal switching, as scripts 37, 38, 42 and 42b. Uses the held-out predictions saved by script 46 (no refitting).
# Switching (step 4): one model is in charge of each choice; b_t(m) is the belief that model m is in charge on trial t,
#   starting at b_2 = w. p(c_t) = sum_m b_t(m) p_m(c_t); post_t(m) = b_t(m) p_m(c_t) / p(c_t);
#   b_t+1 = (1 - delta) post_t + delta w. delta = 1 is the fixed blend, delta = 0 one model per person ("types").
#   w and delta are fit by maximum likelihood (BFGS; delta starts at .05, .5 and .95) on 9/10 of the participants and
#   scored on the other 1/10; participants are assigned to the 10 groups at random (seed 37).
# Shuffled (H4): the same model fit to each person's trials in a random order, the same order for all eight models (seed 38).
# Drifting weights (H7): delta = 1, w_t = (1 - g_t) w_early + g_t w_late, g_t = (1 - exp(-k tau_t)) / (1 - exp(-k)),
#   tau_t = (t - 2) / 98; k kept within -30 and 30. The weights start from a straight-line fit (k = 0), which starts
#   from equal weights; k starts at -3, 0 and 3. The straight line (k = 0) is robustness check 4.
# Personal switching (step 5): a switch redraws using lambda x the person's own weights + (1 - lambda) x w; the person's
#   own weights are fit without the fold being scored, and each session is run once per fold (as script 38).
# Exploratory (registration, exploratory analyses): switching whose redraws use drifting weights (straight line, as
#   script 42), and the fits on everyone (w, delta; early and late weights and k).
# If a fit does not converge in a group (optimizer code not 0), that group is refit with up to 2,000 iterations.
# Usage (from the project folder): Rscript scripts/47_prereg_switching.R <data>
#   data: replication, simulated, replication342, or horstmann (reproduces scripts 37, 38, 42 and 42b, as a check)
library(parallel)
data_name <- commandArgs(trailingOnly = TRUE)[1]
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
n <- max(s$rows$subj)
M <- ncol(s$lpd)
L <- aperm(array(s$lpd, c(99, n, M)), c(2, 1, 3))      # participants x trials (2-100) x models
dimnames(L)[[3]] <- colnames(s$lpd)
cores <- max(1, min(10, detectCores()))
softmax <- function(p) { e <- exp(c(0, p)); e / sum(e) }
# optim with BFGS; refit with up to 2,000 iterations if it did not converge
optim2 <- function(par, fn) {
  opt <- optim(par, fn, method = "BFGS", control = list(maxit = 500))
  if (opt$convergence != 0) opt <- optim(opt$par, fn, method = "BFGS", control = list(maxit = 2000))
  opt
}

# Switching: the forward algorithm (script 37)
hmm_ll <- function(L, w, delta, keep = FALSE) {
  Wm <- matrix(w, dim(L)[1], dim(L)[3], byrow = TRUE)
  B <- Wm
  out <- matrix(NA, dim(L)[1], dim(L)[2])
  post <- if (keep) array(NA, dim(L), dimnames(L))
  for (t in 1:dim(L)[2]) {
    A <- log(B) + L[, t, ]
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]
    out[, t] <- mx + log(rowSums(exp(A - mx)))
    P <- exp(A - out[, t])
    if (keep) post[, t, ] <- P
    B <- (1 - delta) * P + delta * Wm
  }
  if (keep) post else out
}
fit_hmm <- function(L, delta = NA) {
  free <- is.na(delta)
  unpack <- function(par) list(w = setNames(softmax(par[1:7]), dimnames(L)[[3]]), delta = if (free) plogis(par[8]) else delta)
  nll <- function(par) { u <- unpack(par); -sum(hmm_ll(L, u$w, u$delta)) }
  best <- NULL
  for (d0 in if (free) c(-3, 0, 3) else NA) {
    opt <- optim2(c(rep(0, 7), if (free) d0), nll)
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  c(unpack(best$par), conv = best$convergence)
}

set.seed(37)
grp <- sample(rep(1:10, length.out = n))
set.seed(38)
L_shuf <- L
for (i in 1:n) L_shuf[i, , ] <- L[i, sample(99), ]

variants <- list(blend = list(L, 1), types = list(L, 0), switching = list(L, NA), shuffled = list(L_shuf, NA))
res <- mclapply(1:10, function(g) lapply(variants, function(v) {
  fit <- fit_hmm(v[[1]][grp != g, , ], v[[2]])
  list(score = rowSums(hmm_ll(v[[1]][grp == g, , , drop = FALSE], fit$w, fit$delta)), delta = fit$delta, conv = fit$conv)
}), mc.cores = cores)
cv <- matrix(NA, n, 4, dimnames = list(NULL, names(variants)))
for (g in 1:10) for (v in names(variants)) cv[grp == g, v] <- res[[g]][[v]]$score
delta_fit <- t(sapply(res, function(r) sapply(r, `[[`, "delta")))
conv_hmm <- t(sapply(res, function(r) sapply(r, `[[`, "conv")))
round(delta_fit[, c("switching", "shuffled")], 3)      # delta chosen in each group
conv_hmm                                                 # optimizer codes (0 = converged)

# Drifting weights (scripts 42 and 42b)
tau <- (0:98) / 98
curve <- function(k) if (abs(k) < 1e-6) tau else (1 - exp(-k * tau)) / (1 - exp(-k))
drift_lp <- function(L, we, wl, k) {                     # no persistence (script 42b)
  g <- curve(k)
  out <- matrix(NA, dim(L)[1], dim(L)[2])
  for (t in 1:dim(L)[2]) {
    w <- (1 - g[t]) * we + g[t] * wl
    A <- sweep(L[, t, ], 2, log(w), "+")
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]
    out[, t] <- mx + log(rowSums(exp(A - mx)))
  }
  out
}
fit_drift <- function(L, k = NA, start = rep(0, 14)) {
  free <- is.na(k)
  unpack <- function(p) list(we = softmax(p[1:7]), wl = softmax(p[8:14]), k = if (free) 30 * tanh(p[15] / 30) else k, par = p)
  nll <- function(p) { u <- unpack(p); -sum(drift_lp(L, u$we, u$wl, u$k)) }
  best <- NULL
  for (k0 in if (free) c(-3, 0, 3) else NA) {
    opt <- optim2(c(start, if (free) k0), nll)
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  c(unpack(best$par), conv = best$convergence)
}
hmm_drift <- function(L, we, wl, delta) {                # switching with drifting redraws, straight line (script 42)
  out <- matrix(NA, dim(L)[1], dim(L)[2])
  B <- matrix(we, dim(L)[1], M, byrow = TRUE)
  for (t in 1:dim(L)[2]) {
    A <- log(B) + L[, t, ]
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]
    out[, t] <- mx + log(rowSums(exp(A - mx)))
    P <- exp(A - out[, t])
    if (t < dim(L)[2]) {
      Wt <- (1 - tau[t + 1]) * we + tau[t + 1] * wl
      B <- (1 - delta) * P + delta * matrix(Wt, dim(L)[1], M, byrow = TRUE)
    }
  }
  out
}
fit_drift_switching <- function(L, start) {
  unpack <- function(p) list(we = softmax(p[1:7]), wl = softmax(p[8:14]), delta = plogis(p[15]))
  nll <- function(p) { u <- unpack(p); -sum(hmm_drift(L, u$we, u$wl, u$delta)) }
  best <- NULL
  for (d0 in c(-3, 0)) {
    opt <- optim2(c(start, d0), nll)
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  c(unpack(best$par), conv = best$convergence)
}
res_d <- mclapply(1:10, function(g) {
  tr <- L[grp != g, , ]
  te <- L[grp == g, , , drop = FALSE]
  f_line <- fit_drift(tr, 0)
  f_curve <- fit_drift(tr, NA, f_line$par[1:14])
  f_ds <- fit_drift_switching(tr, f_line$par[1:14])
  list(line = rowSums(drift_lp(te, f_line$we, f_line$wl, 0)),
       curve = rowSums(drift_lp(te, f_curve$we, f_curve$wl, f_curve$k)),
       drift_switching = rowSums(hmm_drift(te, f_ds$we, f_ds$wl, f_ds$delta)),
       k = f_curve$k, delta = f_ds$delta, conv = c(line = f_line$conv, curve = f_curve$conv, ds = f_ds$conv))
}, mc.cores = cores)
cv <- cbind(cv, drift_line = NA, drift_curve = NA, drift_switching = NA)
for (g in 1:10) {
  cv[grp == g, "drift_line"] <- res_d[[g]]$line
  cv[grp == g, "drift_curve"] <- res_d[[g]]$curve
  cv[grp == g, "drift_switching"] <- res_d[[g]]$drift_switching
}
k_fit <- sapply(res_d, `[[`, "k")
round(k_fit, 2)                                          # k chosen in each group
t(sapply(res_d, `[[`, "conv"))                           # optimizer codes (0 = converged)

# Personal switching (script 38): every session is run once per fold ("copies"), with that person's own weights for
# that fold; only that fold's trials are scored
copy_subj <- rep(1:n, each = 10)
copy_fold <- rep(1:10, n)
W_copy <- t(sapply(seq_along(copy_subj), function(k) s$w_own[copy_subj[k], copy_fold[k], ]))
scored <- outer(copy_fold, (1:99) %% 10 + 1, "==")
hmm_ll_p <- function(Lc, W) {
  out <- matrix(NA, dim(Lc)[1], dim(Lc)[2])
  B <- W
  for (t in 1:dim(Lc)[2]) {
    A <- log(B) + Lc[, t, ]
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]
    out[, t] <- mx + log(rowSums(exp(A - mx)))
    B <- (1 - attr(W, "delta")) * exp(A - out[, t]) + attr(W, "delta") * W
  }
  out
}
redraw <- function(Wc, w, delta, lambda) {
  W <- lambda * Wc + (1 - lambda) * matrix(w, nrow(Wc), M, byrow = TRUE)
  attr(W, "delta") <- delta
  W
}
score_copies <- function(k, w, delta, lambda) {
  out <- hmm_ll_p(L[copy_subj[k], , , drop = FALSE], redraw(W_copy[k, , drop = FALSE], w, delta, lambda))
  rowsum(rowSums(out * scored[k, , drop = FALSE]), copy_subj[k])[, 1]
}
fit_personal <- function(k) {
  unpack <- function(par) list(w = setNames(softmax(par[1:7]), dimnames(L)[[3]]), delta = plogis(par[8]), lambda = plogis(par[9]))
  nll <- function(par) { u <- unpack(par); -sum(score_copies(k, u$w, u$delta, u$lambda)) }
  best <- NULL
  for (d0 in c(-3, 0)) {
    opt <- optim2(c(rep(0, 7), d0, 0), nll)
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  c(unpack(best$par), conv = best$convergence)
}
# Check: lambda = 1 with delta = 1 gives script 46's personal mixes (should print about 0)
max(abs(score_copies(seq_along(copy_subj), rep(1 / M, M), 1, 1) - s$per[, "personal_mix"]))
res_p <- mclapply(1:10, function(g) {
  fit <- fit_personal(which(grp[copy_subj] != g))
  test <- which(grp[copy_subj] == g)
  list(score = score_copies(test, fit$w, fit$delta, fit$lambda), delta = fit$delta, lambda = fit$lambda, conv = fit$conv)
}, mc.cores = cores)
cv <- cbind(cv, personal_switching = NA)
for (g in 1:10) cv[sort(which(grp == g)), "personal_switching"] <- res_p[[g]]$score
personal_fits <- round(t(sapply(res_p, function(r) c(delta = r$delta, lambda = r$lambda, conv = r$conv))), 3)
personal_fits                                           # delta, lambda and optimizer code in each group

# Description: switching and drift fit on everyone
fit_all <- fit_hmm(L)
drift_all <- fit_drift(L, NA, fit_drift(L, 0)$par[1:14])
post <- hmm_ll(L, fit_all$w, fit_all$delta, keep = TRUE)
in_charge <- apply(post, c(2, 3), mean)                # trials x models: average belief about who is in charge

dir.create("output/prereg", showWarnings = FALSE)
saveRDS(list(data = data_name, cv = cv, grp = grp, delta_fit = delta_fit, conv_hmm = conv_hmm, k_fit = k_fit,
             drift_fits = lapply(res_d, function(r) r[c("k", "delta", "conv")]), personal_fits = personal_fits,
             fit_all = fit_all, drift_all = drift_all[c("we", "wl", "k", "conv")], in_charge = in_charge),
        sprintf("output/prereg/switching_%s.rds", data_name))
