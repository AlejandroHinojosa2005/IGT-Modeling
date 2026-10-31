# ORL (Outcome-Representation Learning): Haines, Vassileva & Ahn (2018)
# ORL negative log-likelihood for one participant
# par   = c(Arew, Apun, K, betaF, betaP), same order as hBayesDM
#       ranges: Arew 0-1, Apun 0-1, K 0-5 (transformed inside to 0-242)
#       betaF & betaP unbounded in the paper (limits set when fitting)
# choice = one participant's deck choice (1-4)
# outcome = that participant's raw net outcomes (divided by 100 inside)
source("scripts/01_load_data.R")

orl_nll <- function(par, choice, outcome) {
  
  Arew <- par[1]
  Apun <- par[2]
  Kprime <- par[3]
  betaF <- par[4]
  betaP <- par[5]
  
  # decay: K = 3^Kprime - 1 (Eq. 11)
  K <- 3^Kprime - 1
  
  # starting values: EV, EF, PS (four zeros each); ll
  Ev <- c(0, 0, 0, 0)
  EF <- c(0, 0, 0, 0)
  PS <- c(0, 0, 0, 0)
  ll <- 0
  
  for (t in 1:length(choice)) {
  
    # 1. which deck was chose,
    j <- choice[t]
    
    # 2. combined value V (Eq. 12); log-softmax with theta = 1
    V <- Ev + EF * betaF + PS * betaP
    logp <- V - max(V) - log(sum(exp(V - max(V))))
    
    # 3. add log probability of the chosen deck to ll
    ll <- ll + logp[j]
    
    # 4. outcome / 100, & its sign
    x <- outcome[t] / 100
    sgn <- sign(x)
    
    # 5. learning rates:chosen deck uses Arew if x >= 0, elseApun;
    #    unchosen decks use the other one
    if (x >= 0) {
      a_chosen <- Arew
      a_unchosen <- Apun
    } else {
      a_chosen <- Apun
      a_unchosen <- Arew
    }
    
    # 6. update the chosen deck's EV (Eq. 7) & EF (Eq. 8)
    Ev[j] <- Ev[j] + a_chosen * (x - Ev[j])
    EF[j] <- EF[j] + a_chosen * (sgn - EF[j])
    
    # 7. update the unchosen deck's EF (Eq. 9), with C = 3
    EF[-j] <- EF[-j] + a_unchosen * (-sgn / 3 - EF[-j])
    
    # 8. update PS (Eq. 10)
    PS <- PS / (1 + K)
    PS[j] <- 1 / (1+K)
  }
  -ll
}

# Simulate one ORL participant (used for parameter recovery).
sim_orl <- function(par, deck_outcomes, n_trials = 100) {
  Arew <- par[1]
  Apun <- par[2]
  Kprime <- par[3]
  betaF <- par[4]
  betaP <- par[5]
  
  K <- 3^Kprime - 1
  Ev <- c(0, 0, 0, 0)
  EF <- c(0, 0, 0, 0)
  PS <- c(0, 0, 0, 0)
  choice <- rep(NA, n_trials)
  outcome <- rep(NA, n_trials)
  
  for (t in 1:n_trials) {
    # combined value (Eq. 12) -> choice probabilities (softmax, theta=1)
    V <- Ev + EF * betaF + PS * betaP
    p <- exp(V - max(V)) / sum(exp(V - max(V)))
    j <- sample(1:4, size = 1, prob=p)
    choice[t] <- j
    outcome[t] <- sample(deck_outcomes[[j]], size = 1)
    
    # learning: identical to orl_nll steps 4-8
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
    PS <- PS / (1 + K)
    PS[j] <- 1 / (1+K)
  }
  list(choice=choice,outcome=outcome)
}

# Multi-start maximum likelihood fit or ORL for one participant.
# Learning-rate lower limits are 0.01, not 0: both at 0 & betaP=0
# nothing in the model can change and optim() stops at chance.
# betaF & betaP limits (+/-10) are my choice; unbounded in the paper.
fit_orl <- function(choice, outcome, n_starts = 20, beta_lim = 10) {
  lower <- c(0.01, 0.01, 0, -beta_lim, -beta_lim)
  upper <- c(1, 1, 5, beta_lim, beta_lim)
  best_value <- Inf
  
  for  (s in 1:n_starts) {
    start <- runif(5, min=lower, max=upper)
    fit <- optim(
      par = start,
      fn = orl_nll,
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