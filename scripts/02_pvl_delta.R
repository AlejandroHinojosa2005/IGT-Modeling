# PVL-Delta model: log-likelihood for one participant
# Model: Ahn et al. (2008); Equations as given in Haines et al.
# Fitting plan: 
# 1. Maximum likelihood estimation: each participant fit separately, using R's optim() function.
# This script builds the likelihood.
# 2. Hierarchical Bayesian estimation with hBayesDM (Ahn, Haines & Zhang, 2017),
# as the stronger estimate and as a cross-check on step 1.

source("scripts/01_load_data.R")

# TEST SETUP: arbitrary values to check code runs
# These are not fitted estimates. Later, i will loop over all 162 participants
# and the parameters will be estimated by optim().
A <- 0.2 # learning rate (0-1)
alpha <- 0.5 # outcome sensitivity (0-1)
lambda <- 1.5 # loss aversion (0-10)
cons <- 1 # consistency (0-5)
i <- 1 # participant index

# Choice sensitivity (computed once;doesn't change across trials)
theta <- 3^cons - 1

# Starting values
Ev <- c(0, 0, 0, 0) # value of decks A, B, C, D
ll <- 0 # running log-likelihood total

for (t in 1:100) {
  # 1. which deck was chosen on this trial
  j <- choice_h[i,t]
  # 2. choice probabilities from the current deck values
  p <- exp(theta*Ev)/sum(exp(theta*Ev))
  # 3. add the log of the chosen deck's probability to the running total
  ll <- ll+log(p[j])
  # 4. this trial's outcome, divided by 100
  x <- net_h[i,t]/100
  # 5. utility of that outcome
  if (x >=0) {
    u <- x^alpha
  } else {
    u <- -lambda * abs(x)^alpha
  }
  # 6. update the chosen deck's value
  Ev[j] <- Ev[j] + A * (u-Ev[j])
}

ll