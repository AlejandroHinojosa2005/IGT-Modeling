# PVL-Delta model: log-likelihood for one participant
# Model: Ahn et al. (2008); Equations as given in Haines et al (2018).
# Fitting plan: 
# 1. Maximum likelihood estimation: each participant fit separately, using R's optim() function.
# This script builds the likelihood.
# 2. Hierarchical Bayesian estimation with hBayesDM (Ahn, Haines & Zhang, 2017),
# as the stronger estimate and as a cross-check on step 1.

source("scripts/01_load_data.R")

# PVL-Delta negative log-likelihood for one participant
# par = c(A, alpha, cons, lambda), same order as hBayesDM
# ranges: A 0-1, alpha 0-2, cons 0-5, lambda 0-10
# choice = one participant's deck choices (1-4)
# outcome = that participant's raw net outcomes (divided by 100 inside)
pvl_delta_nll <- function(par, choice, outcome) {
  A <- par[1]
  alpha <- par[2]
  cons <- par[3]
  lambda <- par[4]
  
  theta <- 3^cons - 1
  Ev <- c(0, 0, 0, 0)
  ll <- 0
  
  for (t in 1:length(choice)) {
    j <- choice[t]
    
    # log-softmax, computed stably (subtracting the max avoids rounding to zero)
    v <- theta*Ev
    logp <- v - max(v) - log(sum(exp(v-max(v))))
    ll <- ll + logp[j]
    x <- outcome[t] / 100
    if (x >= 0) {
      u <- x^alpha
    } else {
      u <- -lambda * abs(x)^alpha
    }
    Ev[j] <- Ev[j] + A * (u-Ev[j])
  }
  
  -ll
}

# Simulate one PVL-Delta participant (used for parameter recovery).
# Same equations as pvl_delta_nll, but the model makes the choices
# instead of scoring them.
# par = c(A, alpha, cons, lambda)
# deck_outcomes = list of 4 vectors: observed net outcomes for each deck.
# # outcomes are drawn from these at random, an approximation
# of the real task's fixed card sequences.
sim_pvl_delta <- function(par, deck_outcomes, n_trials = 100) {
  A <- par[1]
  alpha <- par[2]
  cons <- par[3]
  lambda <- par[4]
  
  theta <- 3^cons-1
  Ev <- c(0, 0, 0,0)
  choice <- rep(NA, n_trials)
  outcome <- rep(NA, n_trials)
  
  for (t in 1:n_trials) {
    v <- theta * Ev
    p <- exp(v-max(v))/sum(exp(v-max(v)))
    j <- sample(1:4, size = 1, prob = p)
    choice[t] <- j
    outcome[t] <- sample(deck_outcomes[[j]], size=1)
    
    x <- outcome[t]/100
    if (x>=0) {
      u <- x^alpha
    } else {
      u <- -lambda * abs(x)^alpha
    }
    Ev[j] <- Ev[j] + A * (u-Ev[j])
  }
  list(choice=choice,outcome=outcome)
}

# Multi-Start Maximum likelihood fit of PVL-Delta for one participant.
# optim() can get stuck, so run it from n_starts random starting points
# and keep best. Lower limits are 0.01, not 0:
# at exactly A=0 & cons =0 model is flat in every direction & optim() stops there.
fit_pvl_delta <- function(choice, outcome, n_starts = 20) {
  lower <- c(.01, .01, .01, .01) # A, alpha, cons, lambda
  upper <- c(1, 2, 5, 10)
  best_value <- Inf

for (s in 1:n_starts) {
  start <- runif(4, min=lower, max=upper)
  fit <- optim(
    par = start,
    fn = pvl_delta_nll,
    choice = choice,
    outcome = outcome,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )
  if (fit$value < best_value) {
    best_value <- fit$value
    best_fit <- fit
  }
}

best_fit
}

# TEST: participant 1 should reproduce 128.0598
set.seed(2026)
fit_pvl_delta(choice_h[1, ], net_h[1, ])