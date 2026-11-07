# Fits ORL + alpha and the variant + alpha to the Horstmann sample (exploratory; see README).
# Usage (from the project folder):  Rscript scripts/28_fit_alpha.R <orl|wl> [quick]
#   quick: first 5 participants only, for testing
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
quick <- length(args) >= 2 && args[2] == "quick"

source("scripts/27_alpha_model.R")

if (model == "orl") {
  base <- readRDS("output/orl_fit_horstmann.rds")
  par_names <- c("Arew", "Apun", "Kprime", "betaF", "betaP")
  nll_fn <- orl_alpha_nll
  lower <- lower_orl_alpha
  upper <- upper_orl_alpha
} else {
  base <- readRDS("output/orl_wl_fit_horstmann.rds")
  par_names <- c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")
  nll_fn <- orl_wl_alpha_nll
  lower <- lower_wl_alpha
  upper <- upper_wl_alpha
}

n <- if (quick) 5 else nrow(choice_h)
res <- data.frame(subj = 1:n)
res[c(par_names, "alpha", "nll", "conv")] <- NA

set.seed(2026)
for (i in 1:n) {
  start <- c(unlist(base[i, par_names]), 1)   # first start: the fit without alpha, alpha = 1
  f <- fit_alpha(nll_fn, choice_h[i, ], net_h[i, ], lower, upper, start)
  res[i, c(par_names, "alpha")] <- f$par
  res$nll[i] <- f$value
  res$conv[i] <- f$convergence
}

# Checks: never worse than the same model without alpha; optimizer messages
all(res$nll <= base$nll[1:n] + 1e-6)
table(res$conv)

dir.create("output/alpha", showWarnings = FALSE)
saveRDS(res, sprintf("output/alpha/alpha_%s%s.rds", model, if (quick) "_quick" else ""))
