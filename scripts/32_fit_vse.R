# Fits the newer models to the Horstmann sample (predictions in the README, written before fitting).
# Usage (from the project folder):  Rscript scripts/32_fit_vse.R <vse|eef|wsls> [quick]
#   vse:  fits VSE, then VSE-split with its first start at the VSE fit and psi = 0
#         (so the split can never fit worse than VSE)
#   quick: first 3 participants only, for testing
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
quick <- length(args) >= 2 && args[2] == "quick"
tag <- if (quick) "_quick" else ""

source("scripts/31_vse_models.R")
n <- if (quick) 3 else nrow(choice_h)
dir.create("output/vse", showWarnings = FALSE)

# Fits every participant with one model; returns a data frame of parameters, nll, conv
fit_all <- function(par_names, fit_one) {
  res <- data.frame(subj = 1:n)
  res[c(par_names, "nll", "conv")] <- NA
  for (i in 1:n) {
    f <- fit_one(i)
    res[i, par_names] <- f$par
    res$nll[i] <- f$value
    res$conv[i] <- f$convergence
  }
  res
}

if (model == "vse") {
  names_vse <- c("theta", "Delta", "alpha", "phi", "beta")
  set.seed(2026)
  vse <- fit_all(names_vse, function(i)
    fit_multi(vse_nll, choice_h[i, ], wi_h[i, ], lo_h[i, ], lower_vse, upper_vse))
  saveRDS(vse, sprintf("output/vse/fit_vse%s.rds", tag))

  set.seed(2026)
  split <- fit_all(c(names_vse, "psi"), function(i)
    fit_multi(vse_split_nll, choice_h[i, ], wi_h[i, ], lo_h[i, ], lower_vse_split, upper_vse_split,
              start_par = c(unlist(vse[i, names_vse]), 0)))
  saveRDS(split, sprintf("output/vse/fit_vse_split%s.rds", tag))

  # Checks: the split never fits worse than VSE; optimizer messages
  print(all(split$nll <= vse$nll + 1e-6))
  print(table(vse$conv))
  print(table(split$conv))
}

if (model == "eef") {
  set.seed(2026)
  eef <- fit_all(c("theta", "lambda", "phi", "beta"), function(i)
    fit_multi(eef_nll, choice_h[i, ], wi_h[i, ], lo_h[i, ], lower_eef, upper_eef))
  saveRDS(eef, sprintf("output/vse/fit_eef%s.rds", tag))
  print(table(eef$conv))
}

if (model == "wsls") {
  wsls <- fit_all(c("p_stay_win", "p_shift_loss"), function(i)
    fit_wsls(choice_h[i, ], wi_h[i, ], lo_h[i, ]))
  saveRDS(wsls, sprintf("output/vse/fit_wsls%s.rds", tag))
}

# Parameters at a limit (for the limit check, as with the earlier models)
res <- readRDS(sprintf("output/vse/fit_%s%s.rds", model, tag))
print(summary(res))
