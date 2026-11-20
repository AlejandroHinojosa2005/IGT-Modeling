# Exploratory, not predicted (Horstmann only): does a curved drift match switching?
# Script 42 tested weights that drift in a straight line from w_early (trial 2) to w_late (trial 100); switching beat it
# by +103 (z = 5.3). A straight line can't change fast early and then level off, the way learning usually looks.
# Here the drift follows a curve with one extra parameter, k:
#   g(tau) = (1 - exp(-k * tau)) / (1 - exp(-k)),   tau = 0 at trial 2 and 1 at trial 100
#   k > 0: the weights change quickly early and then level off; k < 0: slowly early, then quickly; k = 0: the straight line.
# Weights on trial t: w_t = (1 - g) * w_early + g * w_late. No persistence (delta = 1), as for script 42's drifting blend.
# Cross-validated by participant with script 37's groups. Uses the saved held-out predictions; no refitting of the models.
library(parallel)
s <- readRDS("output/stacking.rds")
sw <- readRDS("output/switching.rds")
n <- max(s$rows$subj)
M <- ncol(s$lpd)
L <- aperm(array(s$lpd, c(99, n, M)), c(2, 1, 3))      # participants x trials (2-100) x models
compare <- function(x, a, b) {
  d <- x[, a] - x[, b]
  round(c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0)), 1)
}

tau <- (0:98) / 98
curve <- function(k) if (abs(k) < 1e-6) tau else (1 - exp(-k * tau)) / (1 - exp(-k))
# Log-probability of each held-out choice under the drifting mix (no persistence)
drift_lp <- function(L, we, wl, k) {
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
softmax <- function(p) { e <- exp(c(0, p)); e / sum(e) }
# k = NA: k is fitted (kept within -30 and 30); otherwise k is fixed (0 = the straight line of script 42)
fit_drift <- function(L, k = NA, start = rep(0, 14)) {
  free <- is.na(k)
  unpack <- function(p) list(we = softmax(p[1:7]), wl = softmax(p[8:14]), k = if (free) 30 * tanh(p[15] / 30) else k, par = p)
  nll <- function(p) { u <- unpack(p); -sum(drift_lp(L, u$we, u$wl, u$k)) }
  best <- NULL
  for (k0 in if (free) c(-3, 0, 3) else NA) {          # the curve starts bending down, straight and bending up
    opt <- optim(c(start, if (free) k0), nll, method = "BFGS", control = list(maxit = 500))
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  u <- unpack(best$par); u$conv <- best$convergence; u
}
# Check: k = 0 gives script 42's straight-line drift (the same function, written without the switching machinery)
stopifnot(max(abs(curve(0) - tau)) == 0)

grp <- sw$grp                                          # script 37's participant groups (seed 37)
res <- mclapply(1:10, function(g) {
  f_line <- fit_drift(L[grp != g, , ], 0)              # straight line (script 42's drifting blend)
  f_curve <- fit_drift(L[grp != g, , ], NA, f_line$par[1:14])
  test <- L[grp == g, , , drop = FALSE]
  list(line = rowSums(drift_lp(test, f_line$we, f_line$wl, 0)),
       curve = rowSums(drift_lp(test, f_curve$we, f_curve$wl, f_curve$k)),
       k = f_curve$k, conv = c(f_line$conv, f_curve$conv))
}, mc.cores = min(10, detectCores()))
cv <- sw$cv
cv <- cbind(cv, drift_line = NA, drift_curve = NA)
for (g in 1:10) {
  cv[grp == g, "drift_line"] <- res[[g]]$line
  cv[grp == g, "drift_curve"] <- res[[g]]$curve
}
k_by_group <- round(sapply(res, `[[`, "k"), 2)
k_by_group                                             # k chosen in each group
sapply(res, `[[`, "conv")                              # optimizer codes (0 = converged)
round(colSums(cv), 1)
curved <- rbind(curve_vs_line = compare(cv, "drift_curve", "drift_line"),
                curve_vs_blend = compare(cv, "drift_curve", "blend"),
                switching_vs_curve = compare(cv, "switching", "drift_curve"))
curved

# Description: the curve and the weights fit on everyone
fit_all <- fit_drift(L, NA, fit_drift(L, 0)$par[1:14])
weights_early_late <- round(rbind(early = fit_all$we, late = fit_all$wl), 3)
colnames(weights_early_late) <- colnames(s$lpd)
k_all <- round(fit_all$k, 2)
k_all
round(curve(fit_all$k)[c(1, 11, 26, 50, 75, 99)], 2)   # how far along the drift is at trials 2, 12, 27, 51, 76, 100
weights_early_late

saveRDS(list(cv = cv, curved = curved, k_by_group = k_by_group, k_all = k_all,
             weights_early_late = weights_early_late), "output/drift_curved.rds")
