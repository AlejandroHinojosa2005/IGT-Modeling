# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY): the eight models, shared by scripts 44-48.
# Equations and parameter limits are those of the main fits (scripts 27 and 31); the per-trial log-probability
# functions are copied from script 35 unchanged. For the no-switching simulation (script 44), each model also has a
# step-by-step version (start, deck probabilities, update), checked below against the per-trial functions.
source("scripts/27_alpha_model.R")   # ORL family and its limits (also loads 01, 04, 10)
source("scripts/31_vse_models.R")    # VSE family, its limits, w_ini and deck_draws

# Per-trial log-probability of each choice (copied from script 35)
orl_family_lp <- function(par7, choice, win, loss) {
  Arew <- par7[1]; Apun <- par7[2]; K <- 3^par7[3] - 1
  betaF <- par7[4]; betaPw <- par7[5]; betaPl <- par7[6]; alpha <- par7[7]
  Ev <- EF <- PSw <- PSl <- c(0, 0, 0, 0)
  lp <- rep(NA, length(choice))
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    x <- (win[t] + loss[t]) / 100
    sgn <- sign(x)
    u <- sgn * abs(x)^alpha
    if (x >= 0) { a_chosen <- Arew; a_unchosen <- Apun } else { a_chosen <- Apun; a_unchosen <- Arew }
    Ev[j] <- Ev[j] + a_chosen * (u - Ev[j])
    EF[j] <- EF[j] + a_chosen * (sgn - EF[j])
    EF[-j] <- EF[-j] + a_unchosen * (-sgn / 3 - EF[-j])
    PSw <- PSw / (1 + K)
    PSl <- PSl / (1 + K)
    if (x >= 0) { PSw[j] <- 1 / (1 + K); PSl[j] <- 0 } else { PSl[j] <- 1 / (1 + K); PSw[j] <- 0 }
  }
  lp
}
vse_family_lp <- function(par6, choice, win, loss) {
  theta <- par6[1]; Delta <- par6[2]; alpha <- par6[3]; phi <- par6[4]; C <- 3^par6[5] - 1; psi <- par6[6]
  Exploit <- c(0, 0, 0, 0)
  Explore <- rep(phi, 4)
  lp <- rep(NA, length(choice))
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- (Exploit + Explore) * C
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- Exploit * Delta
    Exploit[j] <- Exploit[j] + v
    Explore[-j] <- Explore[-j] + alpha * (phi - Explore[-j])
    Explore[j] <- if (g - l < 0) psi else 0
  }
  lp
}
eef_lp <- function(par, choice, win, loss) {
  theta <- par[1]; lambda <- par[2]; phi <- par[3]; C <- 3^par[4] - 1
  Exploit <- w_ini
  Explore <- c(0, 0, 0, 0)
  lp <- rep(NA, length(choice))
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- (Exploit + Explore) * C
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- (1 - lambda) * Exploit
    Exploit[j] <- Exploit[j] + v
    Explore[-j] <- lambda * Explore[-j] + (1 - lambda) * phi
    Explore[j] <- 0
  }
  lp
}

# The seven fitted models: per-trial function, parameter names and limits (as in script 35)
models <- list(
  orl = list(lp = function(p, ...) orl_family_lp(c(p[1:5], p[5], 1), ...),
             names = c("Arew", "Apun", "Kprime", "betaF", "betaP"),
             lower = c(.01, .01, 0, -10, -10), upper = c(1, 1, 5, 10, 10)),
  variant = list(lp = function(p, ...) orl_family_lp(c(p, 1), ...),
                 names = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl"),
                 lower = c(.01, .01, 0, -10, -10, -10), upper = c(1, 1, 5, 10, 10, 10)),
  orl_alpha = list(lp = function(p, ...) orl_family_lp(c(p[1:5], p[5], p[6]), ...),
                   names = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha"),
                   lower = lower_orl_alpha, upper = upper_orl_alpha),
  variant_alpha = list(lp = function(p, ...) orl_family_lp(p, ...),
                       names = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl", "alpha"),
                       lower = lower_wl_alpha, upper = upper_wl_alpha),
  vse = list(lp = function(p, ...) vse_family_lp(c(p, 0), ...),
             names = c("theta", "Delta", "alpha", "phi", "beta"),
             lower = lower_vse, upper = upper_vse),
  vse_split = list(lp = function(p, ...) vse_family_lp(p, ...),
                   names = c("theta", "Delta", "alpha", "phi", "beta", "psi"),
                   lower = lower_vse_split, upper = upper_vse_split),
  eef = list(lp = function(p, ...) eef_lp(p, ...),
             names = c("theta", "lambda", "phi", "beta"),
             lower = lower_eef, upper = upper_eef))
model_names <- c(names(models), "wsls")

# Held-out trials of fold f (trial 1 is never held out)
test_trials <- function(f) { tt <- seq(f, 100, by = 10); tt[tt != 1] }

# WSLS held-out log-probabilities, exactly as in scripts 35b and 36: stay/shift proportions from the training
# trials, kept within .01-.99. Returns participants x trials 2-100.
wsls_heldout <- function(choice, win, loss) {
  n <- nrow(choice)
  out <- matrix(NA, n, 99)
  for (f in 1:10) {
    test <- test_trials(f)
    train <- setdiff(2:100, test)
    for (i in 1:n) {
      ch <- choice[i, ]
      stay <- c(NA, ch[-1] == ch[-100])
      gain <- c(NA, (win[i, -100] + loss[i, -100]) >= 0)
      p_sw <- min(max(if (any(gain[train])) mean(stay[train][gain[train]]) else 0.5, .01), .99)
      p_sl <- min(max(if (any(!gain[train])) mean(!stay[train][!gain[train]]) else 0.5, .01), .99)
      p <- ifelse(gain, ifelse(stay, p_sw, (1 - p_sw) / 3), ifelse(stay, 1 - p_sl, p_sl / 3))
      out[i, test - 1] <- log(p[test])
    }
  }
  out
}

# Step-by-step versions, for simulation: start(p) gives the starting state, probs(s, p) each deck's probability,
# update(s, p, j, win, loss) the state after choosing deck j and receiving (win, loss). Same equations as above.
softmax4 <- function(V) { e <- exp(V - max(V)); e / sum(e) }
orl_steps <- function(par7) list(
  start = function() list(Ev = rep(0, 4), EF = rep(0, 4), PSw = rep(0, 4), PSl = rep(0, 4)),
  probs = function(s) softmax4(s$Ev + s$EF * par7[4] + s$PSw * par7[5] + s$PSl * par7[6]),
  update = function(s, j, win, loss) {
    Arew <- par7[1]; Apun <- par7[2]; K <- 3^par7[3] - 1; alpha <- par7[7]
    x <- (win + loss) / 100
    sgn <- sign(x)
    u <- sgn * abs(x)^alpha
    if (x >= 0) { a_chosen <- Arew; a_unchosen <- Apun } else { a_chosen <- Apun; a_unchosen <- Arew }
    s$Ev[j] <- s$Ev[j] + a_chosen * (u - s$Ev[j])
    s$EF[j] <- s$EF[j] + a_chosen * (sgn - s$EF[j])
    s$EF[-j] <- s$EF[-j] + a_unchosen * (-sgn / 3 - s$EF[-j])
    s$PSw <- s$PSw / (1 + K)
    s$PSl <- s$PSl / (1 + K)
    if (x >= 0) { s$PSw[j] <- 1 / (1 + K); s$PSl[j] <- 0 } else { s$PSl[j] <- 1 / (1 + K); s$PSw[j] <- 0 }
    s
  })
vse_steps <- function(par6) list(
  start = function() list(Exploit = rep(0, 4), Explore = rep(par6[4], 4)),
  probs = function(s) softmax4((s$Exploit + s$Explore) * (3^par6[5] - 1)),
  update = function(s, j, win, loss) {
    theta <- par6[1]; Delta <- par6[2]; alpha <- par6[3]; phi <- par6[4]; psi <- par6[6]
    g <- win / 100
    l <- abs(loss) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    s$Exploit <- s$Exploit * Delta
    s$Exploit[j] <- s$Exploit[j] + v
    s$Explore[-j] <- s$Explore[-j] + alpha * (phi - s$Explore[-j])
    s$Explore[j] <- if (g - l < 0) psi else 0
    s
  })
eef_steps <- function(par) list(
  start = function() list(Exploit = w_ini, Explore = rep(0, 4)),
  probs = function(s) softmax4((s$Exploit + s$Explore) * (3^par[4] - 1)),
  update = function(s, j, win, loss) {
    theta <- par[1]; lambda <- par[2]; phi <- par[3]
    g <- win / 100
    l <- abs(loss) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    s$Exploit <- (1 - lambda) * s$Exploit
    s$Exploit[j] <- s$Exploit[j] + v
    s$Explore[-j] <- lambda * s$Explore[-j] + (1 - lambda) * phi
    s$Explore[j] <- 0
    s
  })
# WSLS (script 31's fit_wsls): trial 1 a random pick; then stay with P(stay | last net outcome >= 0) after a gain or 0,
# shift with P(shift | loss) after a loss, a shift spread evenly over the other three decks
wsls_steps <- function(par) list(
  start = function() list(last = NA, gain = NA),
  probs = function(s) {
    if (is.na(s$last)) return(rep(.25, 4))
    p_stay <- if (s$gain) par[1] else 1 - par[2]
    p <- rep((1 - p_stay) / 3, 4)
    p[s$last] <- p_stay
    p
  },
  update = function(s, j, win, loss) list(last = j, gain = win + loss >= 0))
steps_of <- list(orl = function(p) orl_steps(c(p[1:5], p[5], 1)),
                 variant = function(p) orl_steps(c(p, 1)),
                 orl_alpha = function(p) orl_steps(c(p[1:5], p[5], p[6])),
                 variant_alpha = function(p) orl_steps(p),
                 vse = function(p) vse_steps(c(p, 0)),
                 vse_split = function(p) vse_steps(p),
                 eef = function(p) eef_steps(p),
                 wsls = function(p) wsls_steps(p))
# Log-probability of each actual choice, computed step by step (used to check the step versions)
steps_lp <- function(st, choice, win, loss) {
  s <- st$start()
  lp <- rep(NA, length(choice))
  for (t in seq_along(choice)) {
    lp[t] <- log(st$probs(s)[choice[t]])
    s <- st$update(s, choice[t], win[t], loss[t])
  }
  lp
}
