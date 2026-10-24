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
    
    PSw <- PSw / (1 + K) # Chage 4: both fade,
    PSl <- PSl / (1 + K) # chosen deck reset in one trace, zeroed in the other
    if (x >= 0) {
      PSw[j] <- 1 / (1 + K); PSl[j] <- 0
    } else {
      PSl[j] <- 1 / (1 + K); PSw[j] <- 0
    }
  }
  -ll
}









