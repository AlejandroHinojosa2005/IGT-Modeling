# Limit check for VSE, VSE-split and EEF (same idea as script 11b): many phi (and some psi)
# estimates sit at their limits (+/-10; +/-5 for EEF), mostly in people with low consistency.
# Refit only those people with the limits widened 5x, starting from their current fit,
# and see whether the model comparison changes. Sensitivity analysis; not a prediction.
quick <- FALSE   # TRUE: only the first 2 people at a limit per model, for testing
source("scripts/31_vse_models.R")
n <- nrow(choice_h)

vse <- readRDS("output/vse/fit_vse.rds")
split <- readRDS("output/vse/fit_vse_split.rds")
eef <- readRDS("output/vse/fit_eef.rds")
names_vse <- c("theta", "Delta", "alpha", "phi", "beta")
names_split <- c(names_vse, "psi")
names_eef <- c("theta", "lambda", "phi", "beta")

# Wider limits: phi and psi to +/-50 (EEF phi to +/-25); everything else unchanged
lower_vse_w <- replace(lower_vse, 4, -50)
upper_vse_w <- replace(upper_vse, 4, 50)
lower_split_w <- replace(lower_vse_split, c(4, 6), -50)
upper_split_w <- replace(upper_vse_split, c(4, 6), 50)
lower_eef_w <- replace(lower_eef, 3, -25)
upper_eef_w <- replace(upper_eef, 3, 25)

# Refit the people at a limit; keep the new fit only if it is better
refit_at_limit <- function(fits, par_names, nll_fn, lower, upper, at_limit) {
  out <- fits
  who <- which(at_limit)
  if (quick) who <- head(who, 2)
  for (i in who) {
    f <- fit_multi(nll_fn, choice_h[i, ], wi_h[i, ], lo_h[i, ], lower, upper,
                   start_par = unlist(fits[i, par_names]))
    if (f$value < out$nll[i]) {
      out[i, par_names] <- f$par
      out$nll[i] <- f$value
      out$conv[i] <- f$convergence
    }
  }
  out
}

at_vse <- abs(vse$phi) >= 9.99
at_split <- abs(split$phi) >= 9.99 | abs(split$psi) >= 9.99
at_eef <- abs(eef$phi) >= 4.99
c(vse = sum(at_vse), vse_split = sum(at_split), eef = sum(at_eef))   # people refit

set.seed(2026)
vse_w <- refit_at_limit(vse, names_vse, vse_nll, lower_vse_w, upper_vse_w, at_vse)
set.seed(2026)
split_w <- refit_at_limit(split, names_split, vse_split_nll, lower_split_w, upper_split_w, at_split)
set.seed(2026)
eef_w <- refit_at_limit(eef, names_eef, eef_nll, lower_eef_w, upper_eef_w, at_eef)

# Summed AIC and BIC with the original and the widened limits
ic <- function(fits, k) c(AIC = sum(2 * fits$nll + 2 * k), BIC = sum(2 * fits$nll + k * log(100)))
round(rbind(vse = ic(vse, 5), vse_wide = ic(vse_w, 5),
            vse_split = ic(split, 6), vse_split_wide = ic(split_w, 6),
            eef = ic(eef, 4), eef_wide = ic(eef_w, 4)))

# Still at the new limits?
c(vse_phi = sum(abs(vse_w$phi) >= 49.9), split_phi = sum(abs(split_w$phi) >= 49.9),
  split_psi = sum(abs(split_w$psi) >= 49.9), eef_phi = sum(abs(eef_w$phi) >= 24.9))

# Does prediction 2's sign result change?
c(psi_below_0 = sum(split_w$psi < 0), was = sum(split$psi < 0))

saveRDS(list(vse = vse_w, vse_split = split_w, eef = eef_w),
        sprintf("output/vse/limit_check%s.rds", if (quick) "_quick" else ""))
