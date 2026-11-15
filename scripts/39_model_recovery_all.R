# Model recovery for the eight models in the cross-validated comparison: can AIC and BIC tell them apart?
# For one generating model: simulate one synthetic participant per real participant, using that
# participant's fitted parameters (from the main fits), then fit all eight models to each synthetic
# dataset with the same procedure as the main fits and save each model's NLL. Script 39b counts the winners.
# All simulations draw real (win, loss) cards from each deck (deck_draws, script 31), so every model,
# including the exploration models that need gains and losses separately, can be fit to every dataset.
# Fitting: 20 random starts within the main fits' limits; the variant, the alpha models and VSE-split
# also start once from the best fit of the model they extend (as in the main fits). WSLS is closed form.
# Usage (from the project folder):  Rscript scripts/39_model_recovery_all.R <generator> [quick]
#   generator: orl, variant, orl_alpha, variant_alpha, vse, vse_split, eef, wsls
args <- commandArgs(trailingOnly = TRUE)
gen <- args[1]
quick <- length(args) >= 2 && args[2] == "quick"

source("scripts/27_alpha_model.R")   # ORL family (also loads 01, 04, 10)
source("scripts/31_vse_models.R")    # VSE family, w_ini, deck_draws, fit_wsls, sim_vse, sim_vse_split, sim_eef

# Per-trial log-probability of each choice (same equations as scripts 35 and 36)
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

# Simulators for the ORL family and WSLS, drawing (win, loss) cards like sim_vse_split
sim_orl_family <- function(par7, deck_draws, n_trials = 100) {
  Arew <- par7[1]; Apun <- par7[2]; K <- 3^par7[3] - 1
  betaF <- par7[4]; betaPw <- par7[5]; betaPl <- par7[6]; alpha <- par7[7]
  Ev <- EF <- PSw <- PSl <- c(0, 0, 0, 0)
  choice <- win <- loss <- rep(NA, n_trials)
  for (t in 1:n_trials) {
    V <- Ev + EF * betaF + PSw * betaPw + PSl * betaPl
    j <- sample(1:4, size = 1, prob = exp(V - max(V)) / sum(exp(V - max(V))))
    k <- sample(nrow(deck_draws[[j]]), size = 1)
    choice[t] <- j; win[t] <- deck_draws[[j]][k, 1]; loss[t] <- deck_draws[[j]][k, 2]
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
  list(choice = choice, win = win, loss = loss)
}
sim_wsls <- function(par, deck_draws, n_trials = 100) {
  choice <- win <- loss <- rep(NA, n_trials)
  for (t in 1:n_trials) {
    if (t == 1) {
      j <- sample(1:4, size = 1)
    } else {
      gain <- win[t - 1] + loss[t - 1] >= 0
      stay <- runif(1) < (if (gain) par[1] else 1 - par[2])
      j <- if (stay) choice[t - 1] else sample(setdiff(1:4, choice[t - 1]), size = 1)
    }
    k <- sample(nrow(deck_draws[[j]]), size = 1)
    choice[t] <- j; win[t] <- deck_draws[[j]][k, 1]; loss[t] <- deck_draws[[j]][k, 2]
  }
  list(choice = choice, win = win, loss = loss)
}

# The generating models: where the fitted parameters are, and how to simulate from them
gens <- list(
  orl = list(file = "output/orl_fit_horstmann.rds", cols = c("Arew", "Apun", "Kprime", "betaF", "betaP"),
             sim = function(p) sim_orl_family(c(p[1:5], p[5], 1), deck_draws)),
  variant = list(file = "output/orl_wl_fit_horstmann.rds",
                 cols = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl"),
                 sim = function(p) sim_orl_family(c(p, 1), deck_draws)),
  orl_alpha = list(file = "output/alpha/alpha_orl.rds",
                   cols = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha"),
                   sim = function(p) sim_orl_family(c(p[1:5], p[5], p[6]), deck_draws)),
  variant_alpha = list(file = "output/alpha/alpha_wl.rds",
                       cols = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl", "alpha"),
                       sim = function(p) sim_orl_family(p, deck_draws)),
  vse = list(file = "output/vse/fit_vse.rds", cols = c("theta", "Delta", "alpha", "phi", "beta"),
             sim = function(p) sim_vse(p, deck_draws)),
  vse_split = list(file = "output/vse/fit_vse_split.rds", cols = c("theta", "Delta", "alpha", "phi", "beta", "psi"),
                   sim = function(p) sim_vse_split(p, deck_draws)),
  eef = list(file = "output/vse/fit_eef.rds", cols = c("theta", "lambda", "phi", "beta"),
             sim = function(p) sim_eef(p, deck_draws)),
  wsls = list(file = "output/vse/fit_wsls.rds", cols = c("p_stay_win", "p_shift_loss"),
              sim = function(p) sim_wsls(p, deck_draws)))

# Best fit from 20 random starts within the limits, plus an optional start (the nested model's best fit)
fit_best <- function(lp, lower, upper, d, start = NULL) {
  nll <- function(p) -sum(lp(p, d$choice, d$win, d$loss))
  inits <- c(if (!is.null(start)) list(start), lapply(1:20, function(s) runif(length(lower), lower, upper)))
  best <- list(value = Inf, convergence = NA)
  for (st in inits) {
    f <- tryCatch(optim(st, nll, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = 1000)),
                  error = function(e) list(value = Inf))
    if (f$value < best$value) best <- f
  }
  best
}

# Fit all eight models to one dataset; returns each model's NLL and optim's convergence code
fit_all <- function(d) {
  f_orl <- fit_best(function(p, ...) orl_family_lp(c(p[1:5], p[5], 1), ...),
                    c(.01, .01, 0, -10, -10), c(1, 1, 5, 10, 10), d)
  f_var <- fit_best(function(p, ...) orl_family_lp(c(p, 1), ...),
                    c(.01, .01, 0, -10, -10, -10), c(1, 1, 5, 10, 10, 10), d, c(f_orl$par, f_orl$par[5]))
  f_orla <- fit_best(function(p, ...) orl_family_lp(c(p[1:5], p[5], p[6]), ...),
                     lower_orl_alpha, upper_orl_alpha, d, c(f_orl$par, 1))
  f_vara <- fit_best(function(p, ...) orl_family_lp(p, ...), lower_wl_alpha, upper_wl_alpha, d, c(f_var$par, 1))
  f_vse <- fit_best(function(p, ...) vse_family_lp(c(p, 0), ...), lower_vse, upper_vse, d)
  f_split <- fit_best(function(p, ...) vse_family_lp(p, ...), lower_vse_split, upper_vse_split, d, c(f_vse$par, 0))
  f_eef <- fit_best(function(p, ...) eef_lp(p, ...), lower_eef, upper_eef, d)
  f_wsls <- fit_wsls(d$choice, d$win, d$loss)
  fits <- list(orl = f_orl, variant = f_var, orl_alpha = f_orla, variant_alpha = f_vara,
               vse = f_vse, vse_split = f_split, eef = f_eef, wsls = f_wsls)
  c(nll = sapply(fits, function(f) f$value), conv = sapply(fits, function(f) f$convergence))
}

# Check: each per-trial function sums to the original model's log-likelihood (should print about 0)
set.seed(1)
i <- 5
chk <- c(sum(orl_family_lp(c(.2, .05, .5, 2, 1, 1, 1), choice_h[i, ], wi_h[i, ], lo_h[i, ])) +
           orl_nll(c(.2, .05, .5, 2, 1), choice_h[i, ], net_h[i, ]),
         sum(orl_family_lp(c(.2, .05, .5, 2, 1, -1, .6), choice_h[i, ], wi_h[i, ], lo_h[i, ])) +
           orl_wl_alpha_nll(c(.2, .05, .5, 2, 1, -1, .6), choice_h[i, ], net_h[i, ]),
         sum(vse_family_lp(c(.5, .9, .3, 2, 1, -1), choice_h[i, ], wi_h[i, ], lo_h[i, ])) +
           vse_split_nll(c(.5, .9, .3, 2, 1, -1), choice_h[i, ], wi_h[i, ], lo_h[i, ]),
         sum(eef_lp(c(.5, .1, 1, 1), choice_h[i, ], wi_h[i, ], lo_h[i, ])) +
           eef_nll(c(.5, .1, 1, 1), choice_h[i, ], wi_h[i, ], lo_h[i, ]))
print(max(abs(chk)))

# Run
g <- gens[[gen]]
true_par <- as.matrix(readRDS(g$file)[, g$cols])
n <- if (quick) 3 else nrow(true_par)
set.seed(2026 + match(gen, names(gens)))
res <- t(sapply(1:n, function(i) fit_all(g$sim(true_par[i, ]))))
res <- data.frame(subj = 1:n, res)

dir.create("output/model_recovery_all", showWarnings = FALSE)
saveRDS(list(gen = gen, true_par = true_par[1:n, , drop = FALSE], res = res),
        sprintf("output/model_recovery_all/%s%s.rds", gen, if (quick) "_quick" else ""))
k <- c(orl = 5, variant = 6, orl_alpha = 6, variant_alpha = 7, vse = 5, vse_split = 6, eef = 4, wsls = 2)
nll <- as.matrix(res[, paste0("nll.", names(k))])
print(table(factor(names(k)[apply(sweep(2 * nll, 2, 2 * k, "+"), 1, which.min)], levels = names(k))))  # AIC winners
