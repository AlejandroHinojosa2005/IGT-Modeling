# Trial-wise cross-validation with maximum-likelihood fits, for the main model comparison.
# Same fold design as the Stan cross-validation (script 23): fold f holds out trials f, f+10, ..., f+90.
# Each participant's model is fit to the other 90% of their choices; held-out outcomes still update
# the model's values (the model sees what happened, it just isn't scored on those choices while fitting).
# Score = log-probability of the held-out choices. Trial 1 is never held out (most models give
# every deck 1/4 there, and EEF's starting weights were set from these participants' first choices).
# Limits are the same as in the main fits. WSLS is done in script 35b (no optimizer needed).
# Usage (from the project folder):  Rscript scripts/35_cv_mle.R <model> <fold 1-10> [quick]
#   model: orl, variant, orl_alpha, variant_alpha, vse, vse_split, eef
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
fold <- as.integer(args[2])
quick <- length(args) >= 3 && args[3] == "quick"

source("scripts/27_alpha_model.R")   # ORL family (also loads 01, 04, 10)
source("scripts/31_vse_models.R")    # VSE family

# Per-trial log-probability of each choice. All models take (par, choice, win, loss).
# ORL family, one function: par7 = c(Arew, Apun, Kprime, betaF, betaPw, betaPl, alpha);
# ORL is betaPw = betaPl and alpha = 1, the variant is alpha = 1 (equations as in scripts 04, 10, 27)
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
# VSE family (equations as in script 31): par6 = c(theta, Delta, alpha, phi, beta, psi)
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

# Each model: its per-trial function, parameter names, and limits (same as the main fits)
models <- list(
  orl = list(lp = function(p, ...) orl_family_lp(c(p[1:5], p[5], 1), ...),
             names = c("Arew", "Apun", "Kprime", "betaF", "betaP"),
             lower = c(.01, .01, 0, -10, -10), upper = c(1, 1, 5, 10, 10)),
  variant = list(lp = function(p, ...) orl_family_lp(c(p, 1), ...),
                 names = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl"),
                 lower = c(.01, .01, 0, -10, -10, -10), upper = c(1, 1, 5, 10, 10, 10)),
  orl_alpha = list(lp = function(p, ...) orl_family_lp(c(p[1:5], p[5], p[6]), ...),
                   names = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha"),
                   lower = lower_orl_alpha, upper = upper_orl_alpha),
  variant_alpha = list(lp = function(p, ...) orl_family_lp(p, ...),
                       names = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl", "alpha"),
                       lower = lower_wl_alpha, upper = upper_wl_alpha),
  vse = list(lp = function(p, ...) vse_family_lp(c(p, 0), ...),
             names = c("theta", "Delta", "alpha", "phi", "beta"),
             lower = lower_vse, upper = upper_vse),
  vse_split = list(lp = function(p, ...) vse_family_lp(p, ...),
                   names = c("theta", "Delta", "alpha", "phi", "beta", "psi"),
                   lower = lower_vse_split, upper = upper_vse_split),
  eef = list(lp = function(p, ...) eef_lp(p, ...),
             names = c("theta", "lambda", "phi", "beta"),
             lower = lower_eef, upper = upper_eef))

# Check: with every trial included, each per-trial function sums to the original model's
# log-likelihood (should print 0)
set.seed(2026)
nll_orig <- list(orl = function(p, i) orl_nll(p, choice_h[i, ], net_h[i, ]),
                 variant = function(p, i) orl_wl_nll(p, choice_h[i, ], net_h[i, ]),
                 orl_alpha = function(p, i) orl_alpha_nll(p, choice_h[i, ], net_h[i, ]),
                 variant_alpha = function(p, i) orl_wl_alpha_nll(p, choice_h[i, ], net_h[i, ]),
                 vse = function(p, i) vse_nll(p, choice_h[i, ], wi_h[i, ], lo_h[i, ]),
                 vse_split = function(p, i) vse_split_nll(p, choice_h[i, ], wi_h[i, ], lo_h[i, ]),
                 eef = function(p, i) eef_nll(p, choice_h[i, ], wi_h[i, ], lo_h[i, ]))
m <- models[[model]]
chk <- replicate(10, {
  i <- sample(nrow(choice_h), 1)
  p <- runif(length(m$lower), m$lower, m$upper)
  sum(m$lp(p, choice_h[i, ], wi_h[i, ], lo_h[i, ])) + nll_orig[[model]](p, i)
})
print(max(abs(chk)))

# Fit on the training trials, score the held-out trials
test <- seq(fold, 100, by = 10)
test <- test[test != 1]
train <- !(1:100 %in% test)
n <- if (quick) 3 else nrow(choice_h)
res <- data.frame(subj = 1:n)
res[c(m$names, "train_nll", "heldout_ll", "conv")] <- NA

set.seed(2026 + fold)
for (i in 1:n) {
  nll_train <- function(p) -sum(m$lp(p, choice_h[i, ], wi_h[i, ], lo_h[i, ])[train])
  best <- list(value = Inf)
  for (s in 1:20) {
    f <- optim(runif(length(m$lower), m$lower, m$upper), nll_train, method = "L-BFGS-B",
               lower = m$lower, upper = m$upper, control = list(maxit = 1000))
    if (f$value < best$value) best <- f
  }
  res[i, m$names] <- best$par
  res$train_nll[i] <- best$value
  res$heldout_ll[i] <- sum(m$lp(best$par, choice_h[i, ], wi_h[i, ], lo_h[i, ])[!train])
  res$conv[i] <- best$convergence
}

dir.create("output/cv_mle", showWarnings = FALSE)
saveRDS(res, sprintf("output/cv_mle/%s_fold%02d%s.rds", model, fold, if (quick) "_quick" else ""))
print(table(res$conv))
print(sum(res$heldout_ll))
