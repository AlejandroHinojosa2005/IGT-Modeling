# Exploratory, not predicted: stacking the cross-validated predictions (Yao, Vehtari, Simpson & Gelman, 2018).
# Instead of picking one model, stacking finds the mix of models whose combined predictions score best on
# held-out choices: each held-out choice gets probability w1 * p_model1 + w2 * p_model2 + ..., with the
# weights summing to 1. A model that is badly wrong on one choice is then rescued by the others.
# Uses the fits saved by script 35 (no refitting): each fold's fitted parameters are run again to get the
# probability of every held-out choice; these must add up to the totals script 35 saved (checked below).
# Weights are judged honestly: for each fold, they are estimated on the other nine folds' held-out choices
# and scored on that fold's (cross-validated stacking). Done once for everyone (group weights) and once
# per participant (each person's own mix: a learner-type profile).
library(loo)
source("scripts/27_alpha_model.R")   # ORL family (also loads 01, 04, 10)
source("scripts/31_vse_models.R")    # VSE family and w_ini

# Per-trial log-probability of each choice, copied from script 35 (same equations)
orl_family_lp <- function(par7, choice, win, loss) {
  Arew <- par7[1]; Apun <- par7[2]; K <- 3^par7[3] - 1
  betaF <- par7[4]; betaPw <- par7[5]; betaPl <- par7[6]; alpha <- par7[7]
  Ev <- EF <- PSw <- PSl <- c(0, 0, 0, 0)
  lp <- rep(NA, length(choice))
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    x <- (win[t] + loss[t]) / 100
    sgn <- sign(x)
    u <- sgn * abs(x)^alpha
    if (x >= 0) { a_chosen <- Arew; a_unchosen <- Apun } else { a_chosen <- Apun; a_unchosen <- Arew }
    Ev[j] <- Ev[j] + a_chosen * (u - Ev[j])
    EF[j] <- EF[j] + a_chosen * (sgn - EF[j])
    EF[-j] <- EF[-j] + a_unchosen * (-sgn / 3 - EF[-j])
    PSw <- PSw / (1 + K)
    PSl <- PSl / (1 + K)
    if (x >= 0) { PSw[j] <- 1 / (1 + K); PSl[j] <- 0 } else { PSl[j] <- 1 / (1 + K); PSw[j] <- 0 }
  }
  lp
}
vse_family_lp <- function(par6, choice, win, loss) {
  theta <- par6[1]; Delta <- par6[2]; alpha <- par6[3]; phi <- par6[4]; C <- 3^par6[5] - 1; psi <- par6[6]
  Exploit <- c(0, 0, 0, 0)
  Explore <- rep(phi, 4)
  lp <- rep(NA, length(choice))
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- (Exploit + Explore) * C
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- Exploit * Delta
    Exploit[j] <- Exploit[j] + v
    Explore[-j] <- Explore[-j] + alpha * (phi - Explore[-j])
    Explore[j] <- if (g - l < 0) psi else 0
  }
  lp
}
eef_lp <- function(par, choice, win, loss) {
  theta <- par[1]; lambda <- par[2]; phi <- par[3]; C <- 3^par[4] - 1
  Exploit <- w_ini
  Explore <- c(0, 0, 0, 0)
  lp <- rep(NA, length(choice))
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- (Exploit + Explore) * C
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- (1 - lambda) * Exploit
    Exploit[j] <- Exploit[j] + v
    Explore[-j] <- lambda * Explore[-j] + (1 - lambda) * phi
    Explore[j] <- 0
  }
  lp
}

models <- c("orl", "variant", "orl_alpha", "variant_alpha", "vse", "vse_split", "eef")
par_names <- list(orl = c("Arew", "Apun", "Kprime", "betaF", "betaP"),
                  variant = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl"),
                  orl_alpha = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha"),
                  variant_alpha = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl", "alpha"),
                  vse = c("theta", "Delta", "alpha", "phi", "beta"),
                  vse_split = c("theta", "Delta", "alpha", "phi", "beta", "psi"),
                  eef = c("theta", "lambda", "phi", "beta"))
lp_fun <- list(orl = function(p, ...) orl_family_lp(c(p[1:5], p[5], 1), ...),
               variant = function(p, ...) orl_family_lp(c(p, 1), ...),
               orl_alpha = function(p, ...) orl_family_lp(c(p[1:5], p[5], p[6]), ...),
               variant_alpha = function(p, ...) orl_family_lp(p, ...),
               vse = function(p, ...) vse_family_lp(c(p, 0), ...),
               vse_split = function(p, ...) vse_family_lp(p, ...),
               eef = function(p, ...) eef_lp(p, ...))

# One row per held-out choice (participant x trial; trial 1 is never held out), one column per model
n <- nrow(choice_h)
rows <- expand.grid(trial = 2:100, subj = 1:n)
rows$fold <- (rows$trial - 1) %% 10 + 1        # fold f holds out trials f, f + 10, ..., f + 90
lpd <- matrix(NA, nrow(rows), length(models) + 1, dimnames = list(NULL, c(models, "wsls")))
for (m in models) {
  for (f in 1:10) {
    fit <- readRDS(sprintf("output/cv_mle/%s_fold%02d.rds", m, f))
    for (i in 1:n) {
      lp <- lp_fun[[m]](unlist(fit[i, par_names[[m]]]), choice_h[i, ], wi_h[i, ], lo_h[i, ])
      k <- which(rows$subj == i & rows$fold == f)
      lpd[k, m] <- lp[rows$trial[k]]
    }
  }
}
# WSLS, exactly as in script 35b
for (f in 1:10) {
  test <- seq(f, 100, by = 10)
  test <- test[test != 1]
  train <- setdiff(2:100, test)
  for (i in 1:n) {
    ch <- choice_h[i, ]
    stay <- c(NA, ch[-1] == ch[-100])
    gain <- c(NA, net_h[i, -100] >= 0)
    p_sw <- min(max(if (any(gain[train])) mean(stay[train][gain[train]]) else 0.5, .01), .99)
    p_sl <- min(max(if (any(!gain[train])) mean(!stay[train][!gain[train]]) else 0.5, .01), .99)
    p <- ifelse(gain, ifelse(stay, p_sw, (1 - p_sw) / 3), ifelse(stay, 1 - p_sl, p_sl / 3))
    k <- which(rows$subj == i & rows$fold == f)
    lpd[k, "wsls"] <- log(p[rows$trial[k]])
  }
}

# Check: per-participant sums match what script 35b computed (should print about 0)
saved <- readRDS("output/cv_mle_summary.rds")$ll
max(abs(rowsum(lpd, rows$subj) - saved))

# Score of a mix: log of the weighted average probability (computed stably)
mix_score <- function(lp_rows, w) {
  keep <- w > 0
  a <- sweep(lp_rows[, keep, drop = FALSE], 2, log(w[keep]), "+")
  mx <- apply(a, 1, max)
  mx + log(rowSums(exp(a - mx)))
}

# 1. Group weights (one mix for everyone), from all held-out choices
w_group <- as.numeric(stacking_weights(lpd))
names(w_group) <- colnames(lpd)
round(w_group, 3)

# 2. Cross-validated scores: group mix and each person's own mix, weights always fit without the scored fold
cv_group <- cv_person <- rep(NA, nrow(rows))
for (f in 1:10) {
  te <- rows$fold == f
  cv_group[te] <- mix_score(lpd[te, ], as.numeric(stacking_weights(lpd[!te, ])))
}
for (i in 1:n) {
  for (f in 1:10) {
    tr <- rows$subj == i & rows$fold != f
    te <- rows$subj == i & rows$fold == f
    cv_person[te] <- mix_score(lpd[te, , drop = FALSE], as.numeric(stacking_weights(lpd[tr, ])))
  }
}

# 3. Compare with the best single model (ORL + alpha); positive = the mix predicts better
per <- rowsum(cbind(lpd, group_mix = cv_group, person_mix = cv_person), rows$subj)
vs_best <- function(x) {
  d <- per[, x] - per[, "orl_alpha"]
  c(total = sum(per[, x]), diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0))
}
round(rbind(orl_alpha = vs_best("orl_alpha"), group_mix = vs_best("group_mix"), person_mix = vs_best("person_mix")), 1)

# Personal mixes vs. the shared mix (positive = personal mixes predict better)
d <- per[, "person_mix"] - per[, "group_mix"]
round(c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0)), 1)

# The shared mix vs. ORL + alpha without the 5 participants left out in script 35b's sensitivity check
bad <- readRDS("output/cv_mle_sensitivity.rds")$left_out
d <- per[-bad, "group_mix"] - per[-bad, "orl_alpha"]
round(c(diff = sum(d), se = sd(d) * sqrt(n - length(bad)), z = sum(d) / (sd(d) * sqrt(n - length(bad))),
        better_for = sum(d > 0)), 1)

# 4. Each person's own weights, from all their held-out choices: a learner-type profile
w_person <- t(sapply(1:n, function(i) as.numeric(stacking_weights(lpd[rows$subj == i, ]))))
colnames(w_person) <- colnames(lpd)
round(colMeans(w_person), 3)                                   # average weight per model
table(factor(colnames(w_person)[apply(w_person, 1, which.max)], levels = colnames(lpd)))  # largest weight
sum(apply(w_person, 1, max) > .5)                              # people with one model above half

saveRDS(list(lpd = lpd, rows = rows, w_group = w_group, w_person = w_person,
             cv_group = cv_group, cv_person = cv_person),
        "output/stacking.rds")

# 5. Shrinkage (exploratory): pull each person's mix toward the shared mix
# mix = lambda * own weights + (1 - lambda) * group weights. lambda = 0 is the shared mix, 1 is the personal mix.
# As above, all weights are fit without the fold being scored.
lambdas <- seq(0, 1, by = 0.1)
w_grp_f <- t(sapply(1:10, function(f) as.numeric(stacking_weights(lpd[rows$fold != f, ]))))  # group weights without fold f
score <- matrix(NA, nrow(rows), length(lambdas))   # every held-out choice's score under every lambda
for (i in 1:n) {
  for (f in 1:10) {
    tr <- rows$subj == i & rows$fold != f
    te <- rows$subj == i & rows$fold == f
    w_own <- as.numeric(stacking_weights(lpd[tr, ]))
    for (k in seq_along(lambdas)) {
      w <- lambdas[k] * w_own + (1 - lambdas[k]) * w_grp_f[f, ]
      score[te, k] <- mix_score(lpd[te, , drop = FALSE], w)
    }
  }
}
# Checks: lambda = 0 must match the shared mix and lambda = 1 the personal mixes (both should print 0)
max(abs(score[, 1] - cv_group))
max(abs(score[, length(lambdas)] - cv_person))

# Total score for each lambda (higher = better)
fold_tot <- rowsum(score, rows$fold)               # 10 folds x 11 lambdas
round(setNames(colSums(fold_tot), lambdas), 1)
# Pick lambda honestly: each fold uses the lambda that did best on the other nine folds
best_k <- sapply(1:10, function(f) which.max(colSums(fold_tot[-f, ])))
lambdas[best_k]
cv_shrink <- score[cbind(1:nrow(rows), best_k[rows$fold])]
per <- cbind(per, shrink_mix = rowsum(cv_shrink, rows$subj)[, 1])
round(vs_best("shrink_mix"), 1)                    # vs. ORL + alpha
d <- per[, "shrink_mix"] - per[, "group_mix"]      # vs. the shared mix
round(c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0)), 1)