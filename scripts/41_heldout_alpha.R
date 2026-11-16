# Held-out test of outcome curvature (ORL + alpha) on the SteingroverInPrep sample (N = 70).
# ORL + alpha was found on the Horstmann sample, the sample that motivated it; this tests it on the
# held-out sample. Predictions are in the README (written before this script was run); 41b scores them.
# Same procedures and limits as the Horstmann fits (scripts 28 and 32) and cross-validation (script 35).
# Usage (from the project folder):
#   Rscript scripts/41_heldout_alpha.R fit [quick]
#     ORL + alpha and the variant + alpha (first start: that participant's held-out ORL or variant
#     fit with alpha = 1, plus 20 random starts) and VSE (20 random starts)
#   Rscript scripts/41_heldout_alpha.R cv <model> <fold> [quick]
#     ten-fold cross-validation; model: orl, variant, orl_alpha, variant_alpha; fold: 1-10
#   quick: first 3 participants only, for testing
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
quick <- "quick" %in% args
tag <- if (quick) "_quick" else ""

source("scripts/27_alpha_model.R")   # ORL family (also loads 01, 04, 10)
source("scripts/31_vse_models.R")    # VSE family
source("scripts/17_load_heldout.R")  # choice_ho, wi_ho, lo_ho, net_ho
n <- if (quick) 3 else nrow(choice_ho)
dir.create("output/heldout_alpha", showWarnings = FALSE)

if (mode == "fit") {
  orl_ho <- readRDS("output/orl_fit_heldout.rds")
  wl_ho <- readRDS("output/orl_wl_fit_heldout.rds")

  # Adds alpha to a fitted model, every participant; the first start is the fit without alpha
  fit_alpha_all <- function(base, par_names, nll_fn, lower, upper) {
    res <- data.frame(subj = 1:n)
    res[c(par_names, "alpha", "nll", "conv")] <- NA
    for (i in 1:n) {
      start <- c(unlist(base[i, par_names]), 1)
      f <- fit_alpha(nll_fn, choice_ho[i, ], net_ho[i, ], lower, upper, start)
      res[i, c(par_names, "alpha")] <- f$par
      res$nll[i] <- f$value
      res$conv[i] <- f$convergence
    }
    res
  }
  set.seed(2026)
  orl_alpha <- fit_alpha_all(orl_ho, c("Arew", "Apun", "Kprime", "betaF", "betaP"),
                             orl_alpha_nll, lower_orl_alpha, upper_orl_alpha)
  set.seed(2026)
  variant_alpha <- fit_alpha_all(wl_ho, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl"),
                                 orl_wl_alpha_nll, lower_wl_alpha, upper_wl_alpha)

  set.seed(2026)
  vse <- data.frame(subj = 1:n)
  vse[c("theta", "Delta", "alpha", "phi", "beta", "nll", "conv")] <- NA
  for (i in 1:n) {
    f <- fit_multi(vse_nll, choice_ho[i, ], wi_ho[i, ], lo_ho[i, ], lower_vse, upper_vse)
    vse[i, c("theta", "Delta", "alpha", "phi", "beta")] <- f$par
    vse$nll[i] <- f$value
    vse$conv[i] <- f$convergence
  }

  # Checks: adding alpha never makes a fit worse (both should print TRUE); optimizer messages
  print(all(orl_alpha$nll <= orl_ho$nll[1:n] + 1e-6))
  print(all(variant_alpha$nll <= wl_ho$nll[1:n] + 1e-6))
  print(table(orl_alpha$conv))
  print(table(variant_alpha$conv))
  print(table(vse$conv))
  saveRDS(list(orl_alpha = orl_alpha, variant_alpha = variant_alpha, vse = vse),
          sprintf("output/heldout_alpha/fits%s.rds", tag))
}

if (mode == "cv") {
  model <- args[2]
  fold <- as.integer(args[3])

  # Per-trial log-probability for the ORL family (copied from script 35):
  # par7 = c(Arew, Apun, Kprime, betaF, betaPw, betaPl, alpha);
  # ORL is betaPw = betaPl and alpha = 1, the variant is alpha = 1
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
  models <- list(
    orl = list(lp = function(p, ...) orl_family_lp(c(p[1:5], p[5], 1), ...),
               nll = function(p, i) orl_nll(p, choice_ho[i, ], net_ho[i, ]),
               names = c("Arew", "Apun", "Kprime", "betaF", "betaP"),
               lower = c(.01, .01, 0, -10, -10), upper = c(1, 1, 5, 10, 10)),
    variant = list(lp = function(p, ...) orl_family_lp(c(p, 1), ...),
                   nll = function(p, i) orl_wl_nll(p, choice_ho[i, ], net_ho[i, ]),
                   names = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl"),
                   lower = c(.01, .01, 0, -10, -10, -10), upper = c(1, 1, 5, 10, 10, 10)),
    orl_alpha = list(lp = function(p, ...) orl_family_lp(c(p[1:5], p[5], p[6]), ...),
                     nll = function(p, i) orl_alpha_nll(p, choice_ho[i, ], net_ho[i, ]),
                     names = c("Arew", "Apun", "Kprime", "betaF", "betaP", "alpha"),
                     lower = lower_orl_alpha, upper = upper_orl_alpha),
    variant_alpha = list(lp = function(p, ...) orl_family_lp(p, ...),
                         nll = function(p, i) orl_wl_alpha_nll(p, choice_ho[i, ], net_ho[i, ]),
                         names = c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl", "alpha"),
                         lower = lower_wl_alpha, upper = upper_wl_alpha))
  m <- models[[model]]

  # Check: with every trial included, the per-trial function sums to the model's log-likelihood
  # (should print 0, up to rounding)
  set.seed(2026)
  chk <- replicate(10, {
    i <- sample(nrow(choice_ho), 1)
    p <- runif(length(m$lower), m$lower, m$upper)
    sum(m$lp(p, choice_ho[i, ], wi_ho[i, ], lo_ho[i, ])) + m$nll(p, i)
  })
  print(max(abs(chk)))

  # Fit on the training trials, score the held-out trials (trial 1 is never held out)
  test <- seq(fold, 100, by = 10)
  test <- test[test != 1]
  train <- !(1:100 %in% test)
  res <- data.frame(subj = 1:n)
  res[c(m$names, "train_nll", "heldout_ll", "conv")] <- NA

  set.seed(2026 + fold)
  for (i in 1:n) {
    nll_train <- function(p) -sum(m$lp(p, choice_ho[i, ], wi_ho[i, ], lo_ho[i, ])[train])
    best <- list(value = Inf)
    for (s in 1:20) {
      f <- optim(runif(length(m$lower), m$lower, m$upper), nll_train, method = "L-BFGS-B",
                 lower = m$lower, upper = m$upper, control = list(maxit = 1000))
      if (f$value < best$value) best <- f
    }
    res[i, m$names] <- best$par
    res$train_nll[i] <- best$value
    res$heldout_ll[i] <- sum(m$lp(best$par, choice_ho[i, ], wi_ho[i, ], lo_ho[i, ])[!train])
    res$conv[i] <- best$convergence
  }
  saveRDS(res, sprintf("output/heldout_alpha/cv_%s_fold%02d%s.rds", model, fold, tag))
  print(table(res$conv))
  print(sum(res$heldout_ll))
}
