# Newer IGT models, written from the published equations:
#   VSE: Value plus Sequential Exploration (Ligneul, 2019), Eqs. 1-4
#   VSE-split: VSE, but after a net loss the chosen deck's exploration weight resets to psi
#              instead of 0 (my extension; psi = 0 is VSE exactly)
#   EEF: Exploitation and Exploration with Forgetting (Yang, Xie & Wang, 2025), Eqs. 1-7
#   WSLS: win-stay/lose-shift (Worthy, Hawthorne & Otto, 2013)
# Conventions checked against Ligneul's igt-toolbox (GPL; no code copied): outcomes / 100,
# gains and losses bent separately, exploration weights start at phi (the paper text says 1),
# phi bounded to [-10, 10], consistency C = 3^beta - 1 with beta in [0, 5].
# win = that trial's gain (>= 0); loss = that trial's loss (stored as a negative number)
source("scripts/01_load_data.R")

# VSE negative log-likelihood for one participant
# par = c(theta, Delta, alpha, phi, beta)
#   theta: value sensitivity (0-1); Delta: decay of exploitation weights (0-1)
#   alpha: exploration learning rate (0-1); phi: exploration bonus (-10 to 10)
#   beta: consistency (0-5), C = 3^beta - 1
vse_nll <- function(par, choice, win, loss) {
  theta <- par[1]
  Delta <- par[2]
  alpha <- par[3]
  phi <- par[4]
  beta <- par[5]
  C <- 3^beta - 1

  Exploit <- c(0, 0, 0, 0)
  Explore <- rep(phi, 4)   # toolbox convention: every deck starts at the bonus
  ll <- 0

  for (t in 1:length(choice)) {
    # 1. which deck was chosen; log-softmax of the summed weights (Eq. 4)
    j <- choice[t]
    V <- (Exploit + Explore) * C
    logp <- V - max(V) - log(sum(exp(V - max(V))))
    ll <- ll + logp[j]

    # 2. utility: gain and loss bent separately (Eq. 1); 0^theta treated as 0
    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)

    # 3. exploitation: every deck decays, the chosen deck adds v (Eqs. 2.A, 2.B)
    Exploit <- Exploit * Delta
    Exploit[j] <- Exploit[j] + v

    # 4. exploration: unchosen decks move toward phi (Eq. 3.B), chosen deck resets to 0 (Eq. 3.A)
    Explore[-j] <- Explore[-j] + alpha * (phi - Explore[-j])
    Explore[j] <- 0
  }
  -ll
}

# VSE-split: copied from vse_nll; one change, marked below
# par = c(theta, Delta, alpha, phi, beta, psi); psi: reset value after a net loss (-10 to 10)
vse_split_nll <- function(par, choice, win, loss) {
  theta <- par[1]
  Delta <- par[2]
  alpha <- par[3]
  phi <- par[4]
  beta <- par[5]
  psi <- par[6]
  C <- 3^beta - 1

  Exploit <- c(0, 0, 0, 0)
  Explore <- rep(phi, 4)
  ll <- 0

  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- (Exploit + Explore) * C
    logp <- V - max(V) - log(sum(exp(V - max(V))))
    ll <- ll + logp[j]

    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)

    Exploit <- Exploit * Delta
    Exploit[j] <- Exploit[j] + v

    Explore[-j] <- Explore[-j] + alpha * (phi - Explore[-j])
    # change: after a net loss the chosen deck resets to psi; after a net gain (or 0) to 0, as in VSE
    Explore[j] <- if (g - l < 0) psi else 0
  }
  -ll
}

# EEF negative log-likelihood for one participant
# par = c(theta, lambda, phi, beta)
#   theta: value sensitivity (0-1); lambda: forgetting (0-1); phi: exploration bonus (-5 to 5)
#   beta: consistency (0-5), C = 3^beta - 1
# Same as VSE with Delta = alpha = 1 - lambda, except the starting weights:
# exploitation starts at the published first-choice weights (Eq. 7), exploration at 0.
w_ini <- c(0.0184, 0.0020, -0.0050, -0.0155)
eef_nll <- function(par, choice, win, loss) {
  theta <- par[1]
  lambda <- par[2]
  phi <- par[3]
  beta <- par[4]
  C <- 3^beta - 1

  Exploit <- w_ini
  Explore <- c(0, 0, 0, 0)
  ll <- 0

  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- (Exploit + Explore) * C                     # Eq. 6
    logp <- V - max(V) - log(sum(exp(V - max(V))))
    ll <- ll + logp[j]

    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)   # Eq. 2

    Exploit <- (1 - lambda) * Exploit                # Eq. 1 (unchosen decks)
    Exploit[j] <- Exploit[j] + v                     # Eq. 3 (chosen deck)
    Explore[-j] <- lambda * Explore[-j] + (1 - lambda) * phi   # Eq. 5
    Explore[j] <- 0                                  # Eq. 4
  }
  -ll
}

# WSLS: P(stay | previous net outcome >= 0) and P(shift | previous net outcome < 0);
# a shift is spread evenly over the other 3 decks; trial 1 is a random pick.
# The best-fitting values are just the observed proportions, so no optimizer is needed.
fit_wsls <- function(choice, win, loss) {
  net <- win + loss
  m <- length(choice)
  stay <- choice[-1] == choice[-m]
  gain <- net[-m] >= 0
  p_stay_win <- if (any(gain)) mean(stay[gain]) else 0.5
  p_shift_loss <- if (any(!gain)) mean(!stay[!gain]) else 0.5
  p <- ifelse(gain, ifelse(stay, p_stay_win, (1 - p_stay_win) / 3),
                    ifelse(stay, 1 - p_shift_loss, p_shift_loss / 3))
  list(par = c(p_stay_win = p_stay_win, p_shift_loss = p_shift_loss),
       value = -(log(1 / 4) + sum(log(p))), convergence = 0)
}

# Simulate one participant. deck_draws: for each deck, a 2-column matrix of observed
# (win, loss) pairs from that deck; each draw takes one whole pair.
sim_vse_split <- function(par, deck_draws, n_trials = 100) {
  theta <- par[1]
  Delta <- par[2]
  alpha <- par[3]
  phi <- par[4]
  beta <- par[5]
  psi <- par[6]
  C <- 3^beta - 1

  Exploit <- c(0, 0, 0, 0)
  Explore <- rep(phi, 4)
  choice <- win <- loss <- rep(NA, n_trials)

  for (t in 1:n_trials) {
    V <- (Exploit + Explore) * C
    p <- exp(V - max(V)) / sum(exp(V - max(V)))
    j <- sample(1:4, size = 1, prob = p)
    k <- sample(nrow(deck_draws[[j]]), size = 1)
    choice[t] <- j
    win[t] <- deck_draws[[j]][k, 1]
    loss[t] <- deck_draws[[j]][k, 2]

    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- Exploit * Delta
    Exploit[j] <- Exploit[j] + v
    Explore[-j] <- Explore[-j] + alpha * (phi - Explore[-j])
    Explore[j] <- if (g - l < 0) psi else 0
  }
  list(choice = choice, win = win, loss = loss, outcome = win + loss)
}

sim_vse <- function(par, deck_draws, n_trials = 100) sim_vse_split(c(par, 0), deck_draws, n_trials)

sim_eef <- function(par, deck_draws, n_trials = 100) {
  theta <- par[1]
  lambda <- par[2]
  phi <- par[3]
  beta <- par[4]
  C <- 3^beta - 1

  Exploit <- w_ini
  Explore <- c(0, 0, 0, 0)
  choice <- win <- loss <- rep(NA, n_trials)

  for (t in 1:n_trials) {
    V <- (Exploit + Explore) * C
    p <- exp(V - max(V)) / sum(exp(V - max(V)))
    j <- sample(1:4, size = 1, prob = p)
    k <- sample(nrow(deck_draws[[j]]), size = 1)
    choice[t] <- j
    win[t] <- deck_draws[[j]][k, 1]
    loss[t] <- deck_draws[[j]][k, 2]

    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- (1 - lambda) * Exploit
    Exploit[j] <- Exploit[j] + v
    Explore[-j] <- lambda * Explore[-j] + (1 - lambda) * phi
    Explore[j] <- 0
  }
  list(choice = choice, win = win, loss = loss, outcome = win + loss)
}

# Multi-start maximum likelihood (L-BFGS-B), same procedure as the other models;
# if start_par is given it is the first start (used so VSE-split starts at the VSE fit, psi = 0)
fit_multi <- function(nll_fn, choice, win, loss, lower, upper, n_starts = 20, start_par = NULL) {
  best_value <- Inf
  for (s in 1:n_starts) {
    if (s == 1 && !is.null(start_par)) {
      start <- start_par
    } else {
      start <- runif(length(lower), min = lower, max = upper)
    }
    fit <- optim(par = start, fn = nll_fn, choice = choice, win = win, loss = loss,
                 method = "L-BFGS-B", lower = lower, upper = upper,
                 control = list(maxit = 1000))
    if (fit$value < best_value) {
      best_value <- fit$value
      best_fit <- fit
    }
  }
  best_fit
}

lower_vse <- c(0, 0, 0, -10, 0)
upper_vse <- c(1, 1, 1, 10, 5)
lower_vse_split <- c(lower_vse, -10)
upper_vse_split <- c(upper_vse, 10)
lower_eef <- c(0, 0, -5, 0)
upper_eef <- c(1, 1, 5, 5)

# Outcome pools for simulation: every observed (win, loss) pair, by deck
deck_draws <- lapply(1:4, function(d) cbind(wi_h[choice_h == d], lo_h[choice_h == d]))

# Checks (should both print 0):
# 1. VSE-split with psi = 0 is VSE exactly
# 2. EEF equals VSE with Delta = alpha = 1 - lambda when both start from the same weights
set.seed(2026)
chk <- replicate(20, {
  i <- sample(nrow(choice_h), 1)
  p <- runif(5, lower_vse, upper_vse)
  vse_nll(p, choice_h[i, ], wi_h[i, ], lo_h[i, ]) -
    vse_split_nll(c(p, 0), choice_h[i, ], wi_h[i, ], lo_h[i, ])
})
max(abs(chk))
w_ini_saved <- w_ini
w_ini <- c(0, 0, 0, 0)                    # temporarily give EEF VSE's starting exploitation
chk2 <- replicate(20, {
  i <- sample(nrow(choice_h), 1)
  p <- runif(4, lower_eef, upper_eef)
  vse_like <- c(p[1], 1 - p[2], 1 - p[2], p[3], p[4])
  # VSE starts exploration at phi, EEF at 0; compare from trial 2 on by giving both phi = 0
  vse_like[4] <- 0
  eef_nll(c(p[1], p[2], 0, p[4]), choice_h[i, ], wi_h[i, ], lo_h[i, ]) -
    vse_nll(vse_like, choice_h[i, ], wi_h[i, ], lo_h[i, ])
})
w_ini <- w_ini_saved
max(abs(chk2))
