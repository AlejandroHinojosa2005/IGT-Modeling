# Lapse-term check: does the win/loss perseveration asymmetry survive a random-choice term?
# Each model gets epsilon: on any trial, with probability epsilon a deck is picked at random,
#   P(deck) = (1 - epsilon) * softmax + epsilon / 4     (Shteingart, Neiman & Loewenstein, 2013)
# Usage (from the project folder):  Rscript scripts/25_lapse.R <orl|wl> [quick]
#   quick: first 5 participants only, for testing
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
quick <- length(args) >= 2 && args[2] == "quick"

source("scripts/10_orl_variant.R")   # data, ORL and variant functions

# Variant with a lapse term; par = c(Arew, Apun, Kprime, betaF, betaPw, betaPl, eps)
# Copied from orl_wl_nll (script 10); only the choice probability changes
wl_lapse_nll <- function(par, choice, outcome) {
  Arew <- par[1]
  Apun <- par[2]
  Kprime <- par[3]
  betaF <- par[4]
  betaPw <- par[5]
  betaPl <- par[6]
  eps <- par[7]

  K <- 3^Kprime - 1
  Ev <- c(0, 0, 0, 0)
  EF <- c(0, 0, 0, 0)
  PSw <- c(0, 0, 0, 0)
  PSl <- c(0, 0, 0, 0)
  ll <- 0

  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl
    p <- exp(V - max(V)) / sum(exp(V - max(V)))
    ll <- ll + log((1 - eps) * p[j] + eps / 4)   # the change: mix in a random choice
    x <- outcome[t] / 100
    sgn <- sign(x)

    if (x >= 0) {
      a_chosen <- Arew
      a_unchosen <- Apun
    } else {
      a_chosen <- Apun
      a_unchosen <- Arew
    }

    Ev[j] <- Ev[j] + a_chosen * (x - Ev[j])
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

# ORL with a lapse term = the variant with betaPw = betaPl (script 10 showed the two are identical)
# par = c(Arew, Apun, Kprime, betaF, betaP, eps)
orl_lapse_nll <- function(par, choice, outcome) {
  wl_lapse_nll(c(par[1:5], par[5], par[6]), choice, outcome)
}

# Multi-start fit; the first start is the participant's fit without a lapse term (eps = 0),
# so the lapse fit can never be worse than the fit without it
fit_lapse <- function(nll_fn, choice, outcome, lower, upper, start_par, n_starts = 20) {
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

if (model == "orl") {
  base <- readRDS("output/orl_fit_horstmann.rds")
  par_names <- c("Arew", "Apun", "Kprime", "betaF", "betaP")
  nll_fn <- orl_lapse_nll
} else {
  base <- readRDS("output/orl_wl_fit_horstmann.rds")
  par_names <- c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")
  nll_fn <- wl_lapse_nll
}
n_beta <- length(par_names) - 3
lower <- c(0.01, 0.01, 0, rep(-10, n_beta), 0)   # eps from 0 (no lapses) to 1 (all random)
upper <- c(1, 1, 5, rep(10, n_beta), 1)

n <- if (quick) 5 else nrow(choice_h)
res <- data.frame(subj = 1:n)
res[c(par_names, "eps", "nll", "conv")] <- NA

set.seed(2026)
for (i in 1:n) {
  start <- c(unlist(base[i, par_names]), 0)
  f <- fit_lapse(nll_fn, choice_h[i, ], net_h[i, ], lower, upper, start)
  res[i, c(par_names, "eps")] <- f$par
  res$nll[i] <- f$value
  res$conv[i] <- f$convergence
}

# Checks: never worse than the same model without a lapse term; optimizer messages
all(res$nll <= base$nll[1:n] + 1e-6)
table(res$conv)

dir.create("output/lapse", showWarnings = FALSE)
saveRDS(res, sprintf("output/lapse/lapse_%s%s.rds", model, if (quick) "_quick" else ""))
