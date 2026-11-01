# Score the held-out test (SteingroverInPrep, N = 70) prediction by prediction,
# exactly as worded in the held-out test plan (October 30, 2026).
source("scripts/10_orl_variant.R")
source("scripts/17_load_heldout.R")
orl_ho <- readRDS("output/orl_fit_heldout.rds")
wl_ho <- readRDS("output/orl_wl_fit_heldout.rds")
n_ho <- nrow(orl_ho)

# Checks
all(wl_ho$nll <= orl_ho$nll + 1e-8)   # should be TRUE: never worse than ORL
table(orl_ho$conv)
table(wl_ho$conv)

# Prediction 1: summed AIC favors the variant (BIC reported either way)
c(orl_aic = sum(2 * orl_ho$nll + 10), wl_aic = sum(2 * wl_ho$nll + 12),
  orl_bic = sum(2 * orl_ho$nll + 5 * log(100)), wl_bic = sum(2 * wl_ho$nll + 6 * log(100)))

# Prediction 2: variant wins on AIC for at least 21 of 70 (one-sided binomial test against 20%)
wins_aic <- sum(orl_ho$nll - wl_ho$nll > 1)
wins_aic
binom.test(wins_aic, n_ho, p = 0.2, alternative = "greater")$p.value
sum(orl_ho$nll - wl_ho$nll > log(100) / 2)   # BIC wins (reported, not a prediction)

# Prediction 3: WSLS. Paired difference's 95% CI includes 0, and r with observed >= ORL's
wsls_of <- function(choice, outcome) {
  n <- length(choice)
  stay <- choice[-1] == choice[-n]
  gain <- outcome[-n] >= 0
  mean(stay[gain]) - mean(stay[!gain])
}
deck_ho <- list(net_ho[choice_ho == 1], net_ho[choice_ho == 2],
                net_ho[choice_ho == 3], net_ho[choice_ho == 4])   # this sample's own outcomes

wsls_obs_ho <- rep(NA, n_ho)
for (i in 1:n_ho) wsls_obs_ho[i] <- wsls_of(choice_ho[i, ], net_ho[i, ])

est_orl <- as.matrix(orl_ho[, c("Arew", "Apun", "Kprime", "betaF", "betaP")])
est_wl <- as.matrix(wl_ho[, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")])
wsls_orl_ho <- rep(NA, n_ho)
wsls_wl_ho <- rep(NA, n_ho)

set.seed(2026)
for (i in 1:n_ho) {
  vals <- rep(NA, 20)
  for (r in 1:20) {
    s <- sim_orl(est_orl[i, ], deck_ho)
    vals[r] <- wsls_of(s$choice, s$outcome)
  }
  wsls_orl_ho[i] <- mean(vals, na.rm = TRUE)
}
set.seed(2026)
for (i in 1:n_ho) {
  vals <- rep(NA, 20)
  for (r in 1:20) {
    s <- sim_orl_wl(est_wl[i, ], deck_ho)
    vals[r] <- wsls_of(s$choice, s$outcome)
  }
  wsls_wl_ho[i] <- mean(vals, na.rm = TRUE)
}

c(observed = mean(wsls_obs_ho, na.rm = TRUE), orl = mean(wsls_orl_ho), variant = mean(wsls_wl_ho))
t.test(wsls_obs_ho, wsls_wl_ho, paired = TRUE)    # variant: does the 95% CI include 0?
t.test(wsls_obs_ho, wsls_orl_ho, paired = TRUE)   # ORL, for comparison
c(r_variant = cor(wsls_obs_ho, wsls_wl_ho, use = "complete.obs"),
  r_orl = cor(wsls_obs_ho, wsls_orl_ho, use = "complete.obs"))

# Prediction 4: A_pun's same-person r between ORL and the variant falls below .80
p <- c("Arew", "Apun", "Kprime", "betaF")
round(diag(cor(orl_ho[, p], wl_ho[, p])), 2)

# Prediction 5: betaP_win > betaP_loss for at least 43 of 70 (one-sided sign test against 50%)
n_pw <- sum(wl_ho$betaPw > wl_ho$betaPl)
n_pw
binom.test(n_pw, n_ho, p = 0.5, alternative = "greater")$p.value

saveRDS(list(wsls_obs = wsls_obs_ho, wsls_orl = wsls_orl_ho, wsls_wl = wsls_wl_ho),
        "output/heldout_wsls.rds")







