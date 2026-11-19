# Exploratory, not predicted (Horstmann only): two checks on scripts 36-37, run before the replication is preregistered.
# Uses the saved held-out predictions (output/stacking.rds) and the switching results (output/switching.rds); no refitting.
#
# 1. Calibration control. Part of the mix's win over ORL + alpha could be hedging: a single model fit by maximum
#    likelihood is overconfident, and averaging softens that. Here each single model is given its own hedge, a small
#    chance of a random guess: p = (1 - e) * p_model + e / 4, with e chosen like the mix's weights (on the other nine
#    folds' held-out choices, scored on the tenth). If the mix still beats the hedged models, its win is not only hedging.
#
# 2. Drifting weights. The main alternative to switching: people don't jump between strategies, but the balance between
#    them slides over the session (for example, from exploring to exploiting). The weights move in a straight line from
#    w_early (trial 2) to w_late (trial 100). With no persistence (delta = 1) this is a "drifting blend"; with delta free,
#    it is switching whose redraws use the drifting weights. Cross-validated by participant with script 37's groups.
library(loo)
library(parallel)
s <- readRDS("output/stacking.rds")
sw <- readRDS("output/switching.rds")
n <- max(s$rows$subj)
M <- ncol(s$lpd)
L <- aperm(array(s$lpd, c(99, n, M)), c(2, 1, 3))      # participants x trials (2-100) x models
dimnames(L)[[3]] <- colnames(s$lpd)
compare <- function(x, a, b) {
  d <- x[, a] - x[, b]
  round(c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0)), 1)
}

# 1. Calibration control -------------------------------------------------------------------------------------
hedged <- sapply(colnames(s$lpd), function(m) {
  x <- rep(NA, nrow(s$rows))
  for (f in 1:10) {
    te <- s$rows$fold == f
    w <- as.numeric(stacking_weights(cbind(s$lpd[!te, m], log(0.25))))   # weight on the model vs. on a random guess
    x[te] <- log(w[1] * exp(s$lpd[te, m]) + w[2] * 0.25)
  }
  x
})
round(colSums(hedged), 1)                              # each hedged model's held-out total
best_hedged <- colnames(hedged)[which.max(colSums(hedged))]
per <- rowsum(cbind(orl_alpha = s$lpd[, "orl_alpha"], orl_alpha_hedged = hedged[, "orl_alpha"],
                    best_hedged = hedged[, best_hedged], shared = s$cv_group, personal = s$cv_person), s$rows$subj)
calibration <- rbind(hedged_vs_orl_alpha = compare(per, "orl_alpha_hedged", "orl_alpha"),
                     shared_vs_hedged_orl_alpha = compare(per, "shared", "orl_alpha_hedged"),
                     shared_vs_best_hedged = compare(per, "shared", "best_hedged"),
                     personal_vs_shared = compare(per, "personal", "shared"))
best_hedged
calibration

# 2. Drifting weights ----------------------------------------------------------------------------------------
tau <- (0:98) / 98                                      # 0 at trial 2, 1 at trial 100
hmm_drift <- function(L, we, wl, delta) {
  out <- matrix(NA, dim(L)[1], dim(L)[2])
  B <- matrix(we, dim(L)[1], M, byrow = TRUE)          # belief about who is in charge, before trial 2
  for (t in 1:dim(L)[2]) {
    A <- log(B) + L[, t, ]
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]
    out[, t] <- mx + log(rowSums(exp(A - mx)))         # log-probability of this choice
    P <- exp(A - out[, t])                             # belief after seeing the choice
    if (t < dim(L)[2]) {
      W <- (1 - tau[t + 1]) * we + tau[t + 1] * wl     # the weights for the next trial
      B <- (1 - delta) * P + delta * matrix(W, dim(L)[1], M, byrow = TRUE)
    }
  }
  out
}
softmax <- function(p) { e <- exp(c(0, p)); e / sum(e) }
fit_drift <- function(L, delta = NA, start = rep(0, 14)) {
  free <- is.na(delta)
  unpack <- function(p) list(we = softmax(p[1:7]), wl = softmax(p[8:14]), delta = if (free) plogis(p[15]) else delta, par = p)
  nll <- function(p) { u <- unpack(p); -sum(hmm_drift(L, u$we, u$wl, u$delta)) }
  best <- NULL
  for (d0 in if (free) c(-3, 0) else NA) {             # delta starts at .05 and .5
    opt <- optim(c(start, if (free) d0), nll, method = "BFGS", control = list(maxit = 500))
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  unpack(best$par)
}
# Check: with w_early = w_late and delta = 1 this is the plain shared mix (should print about 0)
max(abs(hmm_drift(L, s$w_group, s$w_group, 1) - log(apply(exp(L), c(1, 2), function(p) sum(p * s$w_group)))))

grp <- sw$grp                                          # script 37's participant groups (seed 37)
res <- mclapply(1:10, function(g) {
  f1 <- fit_drift(L[grp != g, , ], 1)                  # drifting blend
  f2 <- fit_drift(L[grp != g, , ], NA, f1$par[1:14])   # switching with drifting redraws (starts from the blend fit)
  list(blend = rowSums(hmm_drift(L[grp == g, , , drop = FALSE], f1$we, f1$wl, 1)),
       switching = rowSums(hmm_drift(L[grp == g, , , drop = FALSE], f2$we, f2$wl, f2$delta)),
       delta = f2$delta)
}, mc.cores = min(10, detectCores()))
cv <- sw$cv
cv <- cbind(cv, drift_blend = NA, drift_switching = NA)
for (g in 1:10) {
  cv[grp == g, "drift_blend"] <- res[[g]]$blend
  cv[grp == g, "drift_switching"] <- res[[g]]$switching
}
round(sapply(res, `[[`, "delta"), 3)                   # delta chosen in each group
round(colSums(cv), 1)
drift <- rbind(drift_vs_blend = compare(cv, "drift_blend", "blend"),
               switching_vs_drift = compare(cv, "switching", "drift_blend"),
               drift_switching_vs_drift = compare(cv, "drift_switching", "drift_blend"),
               drift_switching_vs_switching = compare(cv, "drift_switching", "switching"))
drift

# Description: early and late weights fit on everyone (drifting blend)
fit_all <- fit_drift(L, 1)
weights_early_late <- round(rbind(early = fit_all$we, late = fit_all$wl), 3)
colnames(weights_early_late) <- colnames(s$lpd)
weights_early_late

saveRDS(list(hedged = hedged, calibration = calibration, cv = cv, drift = drift,
             weights_early_late = weights_early_late), "output/drift_calibration.rds")
