# The hybrid model (exploratory): ORL + alpha with the three parts of VSE and EEF that it lacks
# (Hinojosa, 2026, Future Directions). VSE's exploration term needs no addition: it is ORL's perseverance with the
# sign reversed (the correspondence in the Discussion of Hinojosa, 2026).
#
# ORL + alpha, for reference. On each trial the deck values are
#   V = Ev + betaF * EF + betaP * PS,   p(deck) = softmax(V)
# and after choosing deck j with net outcome x (win + loss, in hundreds):
#   u = sign(x) |x|^alpha                          (outcome, bent)
#   Ev[j] <- Ev[j] + A (u - Ev[j])                 (A = Arew after a gain, Apun after a loss; other decks unchanged)
#   EF: win frequency;  PS: perseverance, decaying by 1 / (1 + K) each trial
# and every Ev starts at 0.
#
# The hybrid changes three things, each of which gives back ORL + alpha at one setting:
#   1. Starting value: every deck's Ev starts at Q0. Q0 > 0 is optimism about decks not yet tried, the novelty
#      mechanism of Cockburn et al. (2022). Q0 = 0 is ORL + alpha.
#   2. Decay: after each choice, the three unchosen decks' values shrink toward 0, Ev[k] <- D * Ev[k], as every deck's
#      value fades in VSE and EEF. D = 1 is ORL + alpha. (Decay is toward 0, as in VSE and EEF, so a deck's
#      starting optimism also fades if it is left alone.)
#   3. Bending: ORL + alpha bends the net outcome; VSE bends the win and the loss separately, u = g^alpha - l^alpha
#      (g = win, l = |loss|, in hundreds). The two differ only on trials with both a win and a loss. They are not
#      special cases of each other, so the full hybrid is fit both ways ("net" and "separate").
# Everything else (learning rates chosen by the sign of the net outcome, win frequency, perseverance) is ORL + alpha's.
#
# Four versions, so cross-validation can show what each part adds:
#   orl_alpha_start: + Q0              (7 parameters)
#   orl_alpha_decay: + D               (7)
#   hybrid:          + Q0 + D          (8)
#   hybrid_sep:      + Q0 + D, bending the win and the loss separately (8)
# Limits: ORL + alpha's (script 27); Q0 in [-5, 5]; D in [0, 1].
source("scripts/43b_prereg_models.R")

hybrid_lp <- function(par, choice, win, loss, bend = "net") {
  Arew <- par[1]; Apun <- par[2]; K <- 3^par[3] - 1
  betaF <- par[4]; betaP <- par[5]; alpha <- par[6]; Q0 <- par[7]; D <- par[8]
  Ev <- rep(Q0, 4)                                        # 1. starting value
  EF <- PS <- c(0, 0, 0, 0)
  lp <- rep(NA, length(choice))
  for (t in seq_along(choice)) {
    j <- choice[t]
    V <- Ev + EF * betaF + PS * betaP
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    x <- (win[t] + loss[t]) / 100
    sgn <- sign(x)
    u <- if (bend == "net") sgn * abs(x)^alpha else (win[t] / 100)^alpha - (abs(loss[t]) / 100)^alpha   # 3. bending
    if (x >= 0) { a_chosen <- Arew; a_unchosen <- Apun } else { a_chosen <- Apun; a_unchosen <- Arew }
    Ev[-j] <- D * Ev[-j]                                  # 2. decay of the unchosen decks
    Ev[j] <- Ev[j] + a_chosen * (u - Ev[j])
    EF[j] <- EF[j] + a_chosen * (sgn - EF[j])
    EF[-j] <- EF[-j] + a_unchosen * (-sgn / 3 - EF[-j])
    PS <- PS / (1 + K)
    PS[j] <- 1 / (1 + K)
  }
  lp
}

hybrid_models <- list(
  orl_alpha_start = list(lp = function(p, ...) hybrid_lp(c(p, 1), ...),
                         names = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha", "Q0"),
                         lower = c(lower_orl_alpha, -5), upper = c(upper_orl_alpha, 5)),
  orl_alpha_decay = list(lp = function(p, ...) hybrid_lp(c(p[1:6], 0, p[7]), ...),
                         names = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha", "D"),
                         lower = c(lower_orl_alpha, 0), upper = c(upper_orl_alpha, 1)),
  hybrid = list(lp = function(p, ...) hybrid_lp(p, ...),
                names = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha", "Q0", "D"),
                lower = c(lower_orl_alpha, -5, 0), upper = c(upper_orl_alpha, 5, 1)),
  hybrid_sep = list(lp = function(p, ...) hybrid_lp(p, ..., bend = "separate"),
                    names = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha", "Q0", "D"),
                    lower = c(lower_orl_alpha, -5, 0), upper = c(upper_orl_alpha, 5, 1)))

# Check: at Q0 = 0, D = 1 and net bending, the hybrid is ORL + alpha exactly (script 35's function), for random
# parameters within the limits and every Horstmann participant
check_hybrid <- function() {
  d <- load_set("horstmann")
  set.seed(50)
  worst <- 0
  for (i in 1:nrow(d$choice)) {
    p <- runif(6, lower_orl_alpha, upper_orl_alpha)
    a <- hybrid_lp(c(p, 0, 1), d$choice[i, ], d$win[i, ], d$loss[i, ])
    b <- models$orl_alpha$lp(p, d$choice[i, ], d$win[i, ], d$loss[i, ])
    worst <- max(worst, abs(a - b))
  }
  worst
}
