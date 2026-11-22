# Exploratory (after the preregistered replication, https://doi.org/10.17605/OSF.IO/U86SY): what the switching
# model is picking up. Run on the replication sample and on the simulated switchers of script 44b, whose switches are
# known and happen at random times, so they show what the method reports when switches have nothing to do with
# outcomes. Uses the held-out predictions of script 46 and the switching fits of script 47 (no refitting of models).
#
# 1. When does a switch happen? With the switching model fit on everyone (script 47's fit_all: w and delta), the
#    forward-backward algorithm gives, for each pair of neighboring trials, the probability that the model in charge
#    changed between them, given all of the person's choices. This is compared after a net gain, after a net loss,
#    and after a large net loss (-500 or worse): each person's average after losses minus after gains, then averaged
#    over people (SE across people). For the simulated switchers, the same is done with the true switches.
# 2. Is part of the switching gain runs of random choice? A ninth "model" that picks a deck at random (p = 1/4 on
#    every trial) is added. The blend (delta = 1) and the switching model are cross-validated with it, with the
#    groups and fitting of script 47, and compared with script 47's eight-model versions. If the switching gain
#    shrinks once a random state is available, part of it was stretches of near-random choice (lapses of attention)
#    rather than switching between strategies. Also reported: the random state's weight and how often it is in charge.
# Usage (from the project folder): Rscript scripts/51_switch_triggers.R <data>   (replication or sim_switch)
library(parallel)
data_name <- commandArgs(trailingOnly = TRUE)[1]
source("scripts/43_prereg_data.R")
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
sw <- readRDS(sprintf("output/prereg/switching_%s.rds", data_name))
d <- load_set(data_name)
n <- max(s$rows$subj)
M <- ncol(s$lpd)
L <- aperm(array(s$lpd, c(99, n, M)), c(2, 1, 3))      # participants x trials (2-100) x models, as script 47
dimnames(L)[[3]] <- colnames(s$lpd)
stopifnot(nrow(d$choice) == n)
cores <- max(1, min(10, detectCores()))

# Forward-backward for the switching model: P(model in charge changes between trial t and t + 1 | all choices)
switch_prob <- function(L, w, delta) {
  n <- dim(L)[1]; Tn <- dim(L)[2]; M <- dim(L)[3]
  W <- matrix(w, n, M, byrow = TRUE)
  E <- array(NA, dim(L)); for (t in 1:Tn) E[, t, ] <- exp(L[, t, ] - apply(L[, t, ], 1, max))
  post <- array(NA, dim(L)); B <- W                       # forward: filtered belief after each choice
  for (t in 1:Tn) { A <- B * E[, t, ]; post[, t, ] <- A / rowSums(A); B <- (1 - delta) * post[, t, ] + delta * W }
  beta <- matrix(1, n, M); p_stay <- matrix(NA, n, Tn - 1)
  for (t in (Tn - 1):1) {                                 # backward
    v <- E[, t + 1, ] * beta
    b_t <- (1 - delta) * v + delta * rowSums(v * W)
    p_stay[, t] <- rowSums(post[, t, ] * ((1 - delta) + delta * W) * v) / rowSums(post[, t, ] * b_t)
    beta <- b_t / rowSums(b_t)
  }
  1 - p_stay                                             # column t: between trials t + 1 and t + 2 (trials 2-100)
}
P_sw <- switch_prob(L, sw$fit_all$w, sw$fit_all$delta)
net <- (d$win + d$loss)[, 2:99]                          # the outcome just before each possible switch
after <- function(P) {
  g <- rowMeans(ifelse(net >= 0, P, NA), na.rm = TRUE)
  l <- rowMeans(ifelse(net < 0, P, NA), na.rm = TRUE)
  bl <- rowMeans(ifelse(net <= -500, P, NA), na.rm = TRUE)
  dl <- l - g; db <- bl - g
  c(after_gain = mean(g), after_loss = mean(l, na.rm = TRUE), after_big_loss = mean(bl, na.rm = TRUE),
    loss_minus_gain = mean(dl, na.rm = TRUE), se1 = sd(dl, na.rm = TRUE) / sqrt(sum(!is.na(dl))),
    big_minus_gain = mean(db, na.rm = TRUE), se2 = sd(db, na.rm = TRUE) / sqrt(sum(!is.na(db))),
    people_with_big_loss = sum(!is.na(db)))
}
triggers <- rbind(inferred = after(P_sw))
c(switches_per_person = mean(rowSums(P_sw)), delta = sw$fit_all$delta)
if (!is.null(d$in_charge)) {                             # simulated switchers: the true switches
  true_sw <- 1 * (d$in_charge[, 2:99] != d$in_charge[, 3:100])
  triggers <- rbind(triggers, true = after(true_sw))
  c(true_switches_per_person = mean(rowSums(true_sw)),
    r_inferred_true = cor(as.vector(P_sw), as.vector(true_sw)))
}
round(triggers, 4)

# 2. A random state: the blend and the switching model with a ninth model that guesses (log 1/4 on every trial)
softmax <- function(p) { e <- exp(c(0, p)); e / sum(e) }
optim2 <- function(par, fn) {
  opt <- optim(par, fn, method = "BFGS", control = list(maxit = 500))
  if (opt$convergence != 0) opt <- optim(opt$par, fn, method = "BFGS", control = list(maxit = 2000))
  opt
}
hmm_ll <- function(L, w, delta) {                        # script 47's forward algorithm
  Wm <- matrix(w, dim(L)[1], dim(L)[3], byrow = TRUE)
  B <- Wm
  out <- matrix(NA, dim(L)[1], dim(L)[2])
  for (t in 1:dim(L)[2]) {
    A <- log(B) + L[, t, ]
    mx <- A[cbind(1:nrow(A), max.col(A, "first"))]
    out[, t] <- mx + log(rowSums(exp(A - mx)))
    B <- (1 - delta) * exp(A - out[, t]) + delta * Wm
  }
  out
}
fit_hmm <- function(L, delta = NA) {                     # script 47's fit, for any number of models
  K <- dim(L)[3]; free <- is.na(delta)
  unpack <- function(par) list(w = softmax(par[1:(K - 1)]), delta = if (free) plogis(par[K]) else delta)
  nll <- function(par) { u <- unpack(par); -sum(hmm_ll(L, u$w, u$delta)) }
  # Weights start equal; with the random state (last model), also from equal weights with the random state near 0
  # (weight about 0.0003), because when the random state is not needed BFGS can stop while still pushing its weight
  # toward 0 before the other weights are fit
  w_starts <- list(rep(0, K - 1))
  if (dimnames(L)[[3]][K] == "random") w_starts[[2]] <- c(rep(0, K - 2), -8)
  best <- NULL
  for (w0 in w_starts) for (d0 in if (free) c(-3, 0, 3) else NA) {
    opt <- optim2(c(w0, if (free) d0), nll)
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  c(unpack(best$par), conv = best$convergence)
}
L9 <- array(c(L, rep(log(0.25), n * 99)), c(n, 99, M + 1), dimnames = list(NULL, NULL, c(colnames(s$lpd), "random")))
grp <- sw$grp                                            # script 47's groups
res <- mclapply(1:10, function(g) lapply(list(blend9 = 1, switching9 = NA), function(dl) {
  fit <- fit_hmm(L9[grp != g, , ], dl)
  list(score = rowSums(hmm_ll(L9[grp == g, , , drop = FALSE], fit$w, fit$delta)), w_random = fit$w[M + 1],
       delta = fit$delta, conv = fit$conv)
}), mc.cores = cores)
cv <- cbind(sw$cv[, c("blend", "switching")], blend9 = NA, switching9 = NA)
for (g in 1:10) for (v in c("blend9", "switching9")) cv[grp == g, v] <- res[[g]][[v]]$score
fits9 <- round(t(sapply(res, function(r) c(w_random_blend = r$blend9$w_random, w_random_switching = r$switching9$w_random,
                                         delta = r$switching9$delta, conv = r$blend9$conv + r$switching9$conv))), 3)
fits9
compare <- function(a, b) { x <- cv[, a] - cv[, b]; c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n))) }
random_state <- round(rbind(switching_gain_8 = compare("switching", "blend"),
                            switching_gain_9 = compare("switching9", "blend9"),
                            random_in_blend = compare("blend9", "blend"),
                            random_in_switching = compare("switching9", "switching")), 1)
random_state
all9 <- fit_hmm(L9)                                      # fit on everyone: how often the random state is in charge
B <- matrix(all9$w, n, M + 1, byrow = TRUE); share <- matrix(NA, n, 99)
for (t in 1:99) { A <- B * exp(L9[, t, ] - apply(L9[, t, ], 1, max)); P <- A / rowSums(A); share[, t] <- P[, M + 1]
                  B <- (1 - all9$delta) * P + all9$delta * matrix(all9$w, n, M + 1, byrow = TRUE) }
random_all <- c(w_random = all9$w[M + 1], delta = all9$delta, mean_belief_random = mean(share),
                people_with_random_over_half_on_10_trials = sum(rowSums(share > 0.5) >= 10))
round(random_all, 3)

saveRDS(list(data = data_name, P_switch = P_sw, triggers = triggers, cv = cv, fits9 = fits9,
             random_state = random_state, random_all = random_all, random_share = share),
        sprintf("output/prereg/triggers_%s.rds", data_name))
