source("scripts/01_load_data.R")

# VPP (Value-Plus-Perseverance): Worthy, Pang & Byrne (2013);
# equations as given in Haines et al. (2018)
# VPP negative log-likelihood for one participant
# par = c(A, alpha, cons, lambda, epP, epN, K, w), same order as hBayesDM
#       range: A 0-1, alpha 0-2, cons 0-5, lambda 0-10,
#        epP & epN -5 to 5 (my choice; unbounded in hBayesDM), k 0-1, w 0-1
# choice = one participant's deck choice (1-4)
# outcome = that participant's raw net outcome (divided by 100 inside)

vpp_nll <- function(par, choice, outcome) {
  A <- par[1]
  alpha <- par[2]
  cons <- par[3]
  lambda <- par[4]
  epP <- par[5]
  epN <- par[6]
  K <- par[7]
  w <- par[8]
  
  theta <- 3^cons - 1
  Ev <- c(0, 0, 0, 0)
  Pers <- c(0, 0, 0, 0)
  ll <- 0
  
  for (t in 1:length(choice)) {
    j <- choice[t]
    # log-softmax, computed stably (subtracting the max avoids rounding to zero)
    V <- w * Ev + (1-w) * Pers
    v <- theta * V
    logp <- v - max(v) - log(sum(exp(v-max(v))))
    ll <- ll + logp[j]
    x <- outcome[t]/100
    Pers <- K * Pers
    if (x >=0) {
      u <- x^alpha
      Pers[j] <- Pers[j] + epP
    } else {
      u <- -lambda * abs(x)^alpha
      Pers[j] <- Pers[j] + epN
    }
    Ev[j] <- Ev[j] + A * (u-Ev[j])
  }
  
  -ll
}

# Simulate one VPP participant (used for parameter recovery).

sim_vpp <- function(par, deck_outcomes, n_trials = 100) {
  A <- par[1]
  alpha <- par[2]
  cons <- par[3]
  lambda <- par[4]
  epP <- par[5]
  epN <- par[6]
  K <- par[7]
  w <- par[8]
  
  theta <- 3^cons-1
  Ev <- c(0, 0, 0, 0)
  Pers <- c(0, 0, 0, 0)
  choice <- rep(NA, n_trials)
  outcome <- rep(NA, n_trials)
  
  for (t in 1:n_trials) {
    V <- w * Ev + (1-w) * Pers
    v <- theta * V
    p <- exp(v-max(v))/sum(exp(v-max(v)))
    j <- sample(1:4, size = 1, prob = p)
    choice[t] <- j
    outcome[t] <- sample(deck_outcomes[[j]], size=1)
    
    x <- outcome[t]/100
    Pers <- K * Pers
    if (x>=0) {
      u <- x^alpha
      Pers[j] <- Pers[j] + epP
    } else {
      u <- -lambda * abs(x)^alpha
      Pers[j] <- Pers[j] + epN
    }
    Ev[j] <- Ev[j] + A * (u-Ev[j])
  }
  list(choice=choice,outcome=outcome)
}

# Multi-Start Maximum likelihood fit of VPP for one participant.
# optim() can get stuck, so run it from n_starts random starting points
# and keep best. Lower limits are 0.01, not 0:
# at exactly A=0 & cons =0 model is flat in every direction & optim() stops there.
fit_vpp <- function(choice, outcome, n_starts = 20) {
  lower <- c(.01, .01, .01, .01, -5, -5, .01, .01) # A, alpha, cons, lambda, epP, epN, K, w
  upper <- c(1, 2, 5, 10, 5, 5, 1, 1)
  best_value <- Inf
  
  for (s in 1:n_starts) {
    start <- runif(8, min=lower, max=upper)
    fit <- optim(
      par = start,
      fn = vpp_nll,
      choice = choice,
      outcome = outcome,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 1000)
    )
    if (fit$value < best_value) {
      best_value <- fit$value
      best_fit <- fit
    }
  }
  
  best_fit
}

# TEST: participant 1 should give 128.0598 or lower (VPP contains PVL-Delta as w=1)
set.seed(2026)
fit_vpp(choice_h[1, ], net_h[1, ])