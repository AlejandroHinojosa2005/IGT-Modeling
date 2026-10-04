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
    x <- outcome[t]/100
    if (x >=0) {
      u <- x^alpha
    } else {
      u <- -lambda * abs(x)^alpha
    }
    Ev[j] <- Ev[j] + A * (u-Ev[j])
  }
  
  -ll
}