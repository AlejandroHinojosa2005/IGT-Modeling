# ORL Variant
# ORL with outcome-dependent perseveration: separate trace after gains & after losses
# par = c(Arew, Apun, Kprime, betaF, betaPw, betaPl)
# Reduces exactly to ORL when betaPw == betaPl

# I copy pasted orl_nll function from 04_orl.R and modified it from their

source("scripts/01_load_data.R")
source("scripts/04_orl.R") # source 04 for when betaPw == betaPl so original ORL can be used as a fallback
orl_wl_nll <- function(par, choice, outcome) {
  
  Arew <- par[1]
  Apun <- par[2]
  Kprime <- par[3]
  betaF <- par[4]
  betaPw <- par[5] # change 1: replace betaP with betaPw & betaPl 
  betaPl <- par[6]
  
  K <- 3^Kprime - 1
  
  Ev <- c(0, 0, 0, 0)
  EF <- c(0, 0, 0, 0)
  PSw <- c(0, 0, 0, 0) # replace PS <- c(0, 0, 0, 0) with PSw <- c(0, 0, 0, 0) & PSl <- c(0, 0, 0, 0)
  PSl <- c(0, 0, 0, 0)
  ll <- 0
  
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl # change 3: both traces in V
    logp <- V - max(V) - log(sum(exp(V - max(V))))
    ll <- ll + logp[j]
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
    
    PSw <- PSw / (1 + K) # Change 4: both fade,
    PSl <- PSl / (1 + K) # chosen deck reset in one trace, zeroed in the other
    if (x >= 0) {
      PSw[j] <- 1 / (1 + K); PSl[j] <- 0
    } else {
      PSl[j] <- 1 / (1 + K); PSw[j] <- 0
    }
  }
  -ll
}

orl_fit <- readRDS("output/orl_fit_horstmann.rds")
orl_p <- as.matrix(orl_fit[, c("Arew", "Apun", "Kprime", "betaF", "betaP")])

d <- sapply(1:162, function(i)
  orl_nll(orl_p[i, ], choice_h[i, ], net_h[i, ]) -
    orl_wl_nll(c(orl_p[i, ], orl_p[i, 5]), choice_h[i, ], net_h[i, ]))
max(abs(d))


# Sim_ORL copyed from script 04_orl.R
# Same changes made as in orl_nll above

sim_orl_wl <- function(par, deck_outcomes, n_trials = 100) {
  Arew <- par[1]
  Apun <- par[2]
  Kprime <- par[3]
  betaF <- par[4]
  betaPw <- par[5]
  betaPl <- par[6]
  
  K <- 3^Kprime - 1
  Ev <- c(0, 0, 0, 0)
  EF <- c(0, 0, 0, 0)
  PSw <- c(0, 0, 0, 0) 
  PSl <- c(0, 0, 0, 0)
  choice <- rep(NA, n_trials)
  outcome <- rep(NA, n_trials)
  
  for (t in 1:n_trials) {
    # combined value (Eq. 12) -> choice probabilities (softmax, theta=1)
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl
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
    PSw <- PSw / (1 + K) 
    PSl <- PSl / (1 + K) 
    if (x >= 0) {
      PSw[j] <- 1 / (1 + K); PSl[j] <- 0
    } else {
      PSl[j] <- 1 / (1 + K); PSw[j] <- 0
    }
  }
  list(choice=choice,outcome=outcome)
}

# Copied fit ORL and made necessary adjustments as before

fit_orl_wl <- function(choice, outcome, n_starts = 20, beta_lim = 10, start_par = NULL) {
  lower <- c(0.01, 0.01, 0, -beta_lim, -beta_lim, -beta_lim)
  upper <- c(1, 1, 5, beta_lim, beta_lim, beta_lim)
  best_value <- Inf
  
  for  (s in 1:n_starts) {
    if (s == 1 && !is.null(start_par)) {
      start <- start_par
    } else {
      start <- runif(6, min=lower, max=upper)
    }
    fit <- optim(
      par = start,
      fn = orl_wl_nll,
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











