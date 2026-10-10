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
    # 2. combined value V (Eq. 12); log-sotmax with theta = 1
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