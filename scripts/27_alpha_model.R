# Diminishing sensitivity: the ORL variant with PVL-Delta's outcome curvature (alpha).
# Each outcome x (in hundreds) is transformed before the expected-value update:
#   u = sign(x) * |x|^alpha        (same form as PVL-Delta, script 02; alpha 0.01-2)
# alpha = 1 gives back the variant; alpha < 1 means large outcomes count less than in proportion.
# Win frequency (EF) uses only the sign of x and perseveration ignores size, so both are unchanged.
source("scripts/10_orl_variant.R")

# par = c(Arew, Apun, Kprime, betaF, betaPw, betaPl, alpha)
orl_wl_alpha_nll <- function(par, choice, outcome) {
  Arew <- par[1]
  Apun <- par[2]
  Kprime <- par[3]
  betaF <- par[4]
  betaPw <- par[5]
  betaPl <- par[6]
  alpha <- par[7]

  K <- 3^Kprime - 1
  Ev <- c(0, 0, 0, 0)
  EF <- c(0, 0, 0, 0)
  PSw <- c(0, 0, 0, 0)
  PSl <- c(0, 0, 0, 0)
  ll <- 0

  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl
    logp <- V - max(V) - log(sum(exp(V - max(V))))
    ll <- ll + logp[j]
    x <- outcome[t] / 100
    sgn <- sign(x)
    u <- sgn * abs(x)^alpha                     # the change: curved outcome

    if (x >= 0) {
      a_chosen <- Arew
      a_unchosen <- Apun
    } else {
      a_chosen <- Apun
      a_unchosen <- Arew
    }

    Ev[j] <- Ev[j] + a_chosen * (u - Ev[j])     # was (x - Ev[j])
    EF[j] <- EF[j] + a_chosen * (sgn - EF[j])
    EF[-j] <- EF[-j] + a_unchosen * (-sgn / 3 - EF[-j])

    PSw <- PSw / (1 + K)
    PSl <- PSl / (1 + K)
    if (x >= 0) {
      PSw[j] <- 1 / (1 + K); PSl[j] <- 0
    } else {
      PSl[j] <- 1 / (1 + K); PSw[j] <- 0
    }
  }
  -ll
}

# ORL with alpha = the variant with betaPw = betaPl; par = c(Arew, Apun, Kprime, betaF, betaP, alpha)
orl_alpha_nll <- function(par, choice, outcome) {
  orl_wl_alpha_nll(c(par[1:5], par[5], par[6]), choice, outcome)
}

# Simulator for the variant with alpha (copied from sim_orl_wl, script 10)
sim_orl_wl_alpha <- function(par, deck_outcomes, n_trials = 100) {
  Arew <- par[1]
  Apun <- par[2]
  Kprime <- par[3]
  betaF <- par[4]
  betaPw <- par[5]
  betaPl <- par[6]
  alpha <- par[7]

  K <- 3^Kprime - 1
  Ev <- c(0, 0, 0, 0)
  EF <- c(0, 0, 0, 0)
  PSw <- c(0, 0, 0, 0)
  PSl <- c(0, 0, 0, 0)
  choice <- rep(NA, n_trials)
  outcome <- rep(NA, n_trials)

  for (t in 1:n_trials) {
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl
    p <- exp(V - max(V)) / sum(exp(V - max(V)))
    j <- sample(1:4, size = 1, prob = p)
    choice[t] <- j
    outcome[t] <- sample(deck_outcomes[[j]], size = 1)

    x <- outcome[t] / 100
    sgn <- sign(x)
    u <- sgn * abs(x)^alpha
    if (x >= 0) {
      a_chosen <- Arew
      a_unchosen <- Apun
    } else {
      a_chosen <- Apun
      a_unchosen <- Arew
    }
    Ev[j] <- Ev[j] + a_chosen * (u - Ev[j])
    EF[j] <- EF[j] + a_chosen * (sgn - EF[j])
    EF[-j] <- EF[-j] + a_unchosen * (-sgn / 3 - EF[-j])
    PSw <- PSw / (1 + K)
    PSl <- PSl / (1 + K)
    if (x >= 0) {
      PSw[j] <- 1 / (1 + K); PSl[j] <- 0
    } else {
      PSl[j] <- 1 / (1 + K); PSw[j] <- 0
    }
  }
  list(choice = choice, outcome = outcome)
}

# Multi-start fit; the first start is given (the fit without alpha, with alpha = 1),
# so adding alpha can never make a fit worse
fit_alpha <- function(nll_fn, choice, outcome, lower, upper, start_par, n_starts = 20) {
  best_value <- Inf
  for (s in 1:n_starts) {
    if (s == 1) {
      start <- start_par
    } else {
      start <- runif(length(lower), min = lower, max = upper)
    }
    fit <- optim(par = start, fn = nll_fn, choice = choice, outcome = outcome,
                 method = "L-BFGS-B", lower = lower, upper = upper,
                 control = list(maxit = 1000))
    if (fit$value < best_value) {
      best_value <- fit$value
      best_fit <- fit
    }
  }
  best_fit
}
lower_wl_alpha <- c(.01, .01, 0, -10, -10, -10, .01)
upper_wl_alpha <- c(1, 1, 5, 10, 10, 10, 2)
lower_orl_alpha <- c(.01, .01, 0, -10, -10, .01)
upper_orl_alpha <- c(1, 1, 5, 10, 10, 2)

# Check: with alpha = 1 the new model must equal the variant (difference should be ~0)
results_wl_check <- readRDS("output/orl_wl_fit_horstmann.rds")
d_alpha <- rep(NA, 10)
for (i in 1:10) {
  p <- unlist(results_wl_check[i, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")])
  d_alpha[i] <- orl_wl_alpha_nll(c(p, 1), choice_h[i, ], net_h[i, ]) - orl_wl_nll(p, choice_h[i, ], net_h[i, ])
}
max(abs(d_alpha))
