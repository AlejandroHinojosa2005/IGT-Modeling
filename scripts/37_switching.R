# Exploratory: do people switch between strategies during the task?
# Uses the 8 models' held-out predictions saved by script 36 (no refitting). A hidden Markov model: on each
# trial one model is "in charge" of the choice. After each trial, the same model stays in charge with
# probability 1 - delta; otherwise a model is drawn again using the group weights w (which can redraw the same one).
#   delta = 1: a fresh draw every trial, so no memory. This is the shared stacking mix ("blend").
#   delta = 0: each person keeps one model for the whole session ("types").
#   in between: strategies that last a while, then switch ("switching").
# Who is in charge is never observed; it is inferred from the choices so far. So each choice is predicted
# only from earlier choices (plus the model fits, as in scripts 35-36).
# w and delta are fit on 9/10 of the participants and scored on the other 1/10 (10-fold, by participant).
s <- readRDS("output/stacking.rds")
n <- max(s$rows$subj)
L <- aperm(array(s$lpd, c(99, n, ncol(s$lpd))), c(2, 1, 3))   # participants x trials (2-100) x models
dimnames(L)[[3]] <- colnames(s$lpd)
# Check: L holds the same numbers as lpd (should print 0)
max(abs(L[5, 10, ] - s$lpd[s$rows$subj == 5 & s$rows$trial == 11, ]))

# The "forward algorithm": log-probability of each choice given only the earlier choices
hmm_ll <- function(L, w, delta, keep = FALSE) {
  W <- matrix(w, dim(L)[1], dim(L)[3], byrow = TRUE)
  B <- W                                               # belief about who is in charge, before trial 2
  out <- matrix(NA, dim(L)[1], dim(L)[2])
  post <- if (keep) array(NA, dim(L), dimnames(L))
  for (t in 1:dim(L)[2]) {
    A <- log(B) + L[, t, ]                             # log of (belief x probability of the actual choice)
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]     # largest value in each row (for stable sums)
    out[, t] <- mx + log(rowSums(exp(A - mx)))         # log-probability of this choice
    P <- exp(A - out[, t])                             # belief after seeing the choice
    if (keep) post[, t, ] <- P
    B <- (1 - delta) * P + delta * W                   # before the next trial: maybe a switch
  }
  if (keep) post else out
}

# Fit w (and delta, unless fixed) by maximum likelihood
fit_hmm <- function(L, delta = NA) {
  free <- is.na(delta)
  unpack <- function(par) {
    w <- exp(c(0, par[1:7]))
    list(w = setNames(w / sum(w), dimnames(L)[[3]]), delta = if (free) plogis(par[8]) else delta)
  }
  nll <- function(par) { u <- unpack(par); -sum(hmm_ll(L, u$w, u$delta)) }
  best <- NULL
  for (d0 in if (free) c(-3, 0, 3) else NA) {          # delta starts at .05, .5, .95 (plogis turns -3 into .05)
    opt <- optim(c(rep(0, 7), if (free) d0), nll, method = "BFGS", control = list(maxit = 500))
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  unpack(best$par)
}

# Check: with delta = 1 the model is the plain stacking mix (should print about 0)
w_test <- s$w_group
max(abs(hmm_ll(L, w_test, 1) - log(apply(exp(L), c(1, 2), function(p) sum(p * w_test)))))

# Cross-validation by participant
set.seed(37)
grp <- sample(rep(1:10, length.out = n))               # each participant is put in one of 10 groups
variants <- c(blend = 1, types = 0, switching = NA)
cv <- matrix(NA, n, 3, dimnames = list(NULL, names(variants)))
delta_fit <- matrix(NA, 10, 3, dimnames = list(NULL, names(variants)))
for (g in 1:10) {
  for (v in names(variants)) {
    fit <- fit_hmm(L[grp != g, , ], variants[[v]])
    cv[grp == g, v] <- rowSums(hmm_ll(L[grp == g, , , drop = FALSE], fit$w, fit$delta))
    delta_fit[g, v] <- fit$delta
  }
}
round(delta_fit[, "switching"], 3)                     # delta chosen in each group

# Results: totals, and pairwise comparisons (positive = the first predicts better)
cv <- cbind(cv, orl_alpha = rowSums(L[, , "orl_alpha"]))
round(colSums(cv), 1)
compare <- function(a, b) {
  d <- cv[, a] - cv[, b]
  round(c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0)), 1)
}
round(rbind(switching_vs_blend = compare("switching", "blend"),
            switching_vs_types = compare("switching", "types"),
            types_vs_blend = compare("types", "blend"),
            switching_vs_orl_alpha = compare("switching", "orl_alpha")), 1)

# Check: does the order of trials matter? Shuffle each person's trials (the same shuffle for all 8 models).
# This keeps how well each model fits each person but breaks up any runs in time. If the switching model
# gains from runs in time, it should do worse on the shuffled data.
set.seed(38)
L_shuf <- L
for (i in 1:n) L_shuf[i, , ] <- L[i, sample(99), ]
cv <- cbind(cv, shuffled = NA)
for (g in 1:10) {
  fit <- fit_hmm(L_shuf[grp != g, , ])
  cv[grp == g, "shuffled"] <- rowSums(hmm_ll(L_shuf[grp == g, , , drop = FALSE], fit$w, fit$delta))
}
round(rbind(switching_vs_shuffled = compare("switching", "shuffled"),
            shuffled_vs_blend = compare("shuffled", "blend")), 1)

# Description: fit on everyone, then the average belief about who is in charge on each trial
fit_all <- fit_hmm(L)
round(fit_all$w, 3)
round(fit_all$delta, 3)                                # on average a redraw every 1/delta trials
post <- hmm_ll(L, fit_all$w, fit_all$delta, keep = TRUE)
in_charge <- apply(post, c(2, 3), mean)                # trials x models
png("output/figures/switching.png", width = 1600, height = 1000, res = 200)
matplot(2:100, in_charge, type = "l", lty = 1, lwd = 2, col = 1:8, xlab = "Trial",
        ylab = "Average belief: model in charge", ylim = c(0, max(in_charge)))
legend("topright", colnames(in_charge), col = 1:8, lty = 1, lwd = 2, cex = .7, bty = "n")
dev.off()

saveRDS(list(cv = cv, grp = grp, delta_fit = delta_fit, fit_all = fit_all, in_charge = in_charge),
        "output/switching.rds")