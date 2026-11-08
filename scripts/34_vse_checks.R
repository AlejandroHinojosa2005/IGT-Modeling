# Scores the six predictions for the newer models (written before fitting; see README),
# plus the full model comparison. Needs scripts 32 (vse, eef, wsls) and 33 (vse_split) to have run,
# and the alpha fits (script 28) for the comparison table.
quick <- FALSE
tag <- if (quick) "_quick" else ""
source("scripts/31_vse_models.R")
source("scripts/10_orl_variant.R")   # sim_orl and sim_orl_wl, for the simulations
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2],
                      net_h[choice_h == 3], net_h[choice_h == 4])

vse <- readRDS(sprintf("output/vse/fit_vse%s.rds", tag))
split <- readRDS(sprintf("output/vse/fit_vse_split%s.rds", tag))
eef <- readRDS(sprintf("output/vse/fit_eef%s.rds", tag))
wsls <- readRDS(sprintf("output/vse/fit_wsls%s.rds", tag))
n <- nrow(vse)
pvl <- readRDS("output/pvl_delta_fit_horstmann.rds")[1:n, ]
vpp <- readRDS("output/vpp_fit_horstmann.rds")[1:n, ]
orl <- readRDS("output/orl_fit_horstmann.rds")[1:n, ]
wl <- readRDS("output/orl_wl_fit_horstmann.rds")[1:n, ]
oa <- readRDS("output/alpha/alpha_orl.rds")[1:n, ]
wa <- readRDS("output/alpha/alpha_wl.rds")[1:n, ]

# Model comparison: summed AIC and BIC (n = 100 trials per person), best model per person
models <- list(pvl_delta = list(pvl, 4), vpp = list(vpp, 8), orl = list(orl, 5),
               variant = list(wl, 6), orl_alpha = list(oa, 6), variant_alpha = list(wa, 7),
               vse = list(vse, 5), vse_split = list(split, 6), eef = list(eef, 4), wsls = list(wsls, 2))
aic_i <- sapply(models, function(m) 2 * m[[1]]$nll + 2 * m[[2]])
bic_i <- sapply(models, function(m) 2 * m[[1]]$nll + m[[2]] * log(100))
tab <- data.frame(k = sapply(models, function(m) m[[2]]),
                  AIC = round(colSums(aic_i)), BIC = round(colSums(bic_i)),
                  best_AIC = tabulate(apply(aic_i, 1, which.min), length(models)),
                  best_BIC = tabulate(apply(bic_i, 1, which.min), length(models)))
tab[order(tab$AIC), ]
aic <- colSums(aic_i)
bic <- colSums(bic_i)

# 1. VSE beats ORL on summed AIC
pred1 <- aic["vse"] < aic["orl"]

# 2. The split improves VSE: summed AIC, and psi < 0 for most participants (one-sided sign test;
#    psi exactly 0 counts against)
summary(split$psi)
c(below_0 = sum(split$psi < 0), exactly_0 = sum(split$psi == 0), above_0 = sum(split$psi > 0))
sum(vse$nll - split$nll > 1)               # people for whom the split wins on AIC
sum(vse$nll - split$nll > log(100) / 2)    # ... on BIC
p2 <- binom.test(sum(split$psi < 0), n, alternative = "greater")
p2
pred2 <- aic["vse_split"] < aic["vse"] & p2$p.value < .05

# Simulations: 20 datasets per person from their fitted parameters; from each dataset,
# win-stay/lose-shift (script 14), the SE index (Ligneul, 2019) and the deck B curve (script 22)
lags <- 20
wsls_of <- function(choice, outcome) {
  m <- length(choice)
  stay <- choice[-1] == choice[-m]
  gain <- outcome[-m] >= 0
  mean(stay[gain]) - mean(stay[!gain])
}
# SE index: share of 4-trial windows (1-4, 2-5, ...) in which all 4 decks were chosen
se_index <- function(choice) {
  m <- length(choice)
  mean(sapply(1:(m - 3), function(t) length(unique(choice[t:(t + 3)])) == 4))
}
after_B <- function(choice, outcome, target) {
  events <- which(choice == 2 & outcome == target)
  events <- events[events <= length(choice) - lags]
  if (length(events) == 0) return(rep(NA, lags))
  back <- matrix(NA, length(events), lags)
  for (e in 1:length(events)) back[e, ] <- choice[events[e] + 1:lags] == 2
  colMeans(back)
}
sim_stats <- function(sim_fun, est, pool) {
  out <- list(wsls = rep(NA, n), se = rep(NA, n), bloss = matrix(NA, n, lags))
  for (i in 1:n) {
    w <- s_e <- rep(NA, 20)
    reps <- matrix(NA, 20, lags)
    for (r in 1:20) {
      s <- sim_fun(est[i, ], pool)
      w[r] <- wsls_of(s$choice, s$outcome)
      s_e[r] <- se_index(s$choice)
      reps[r, ] <- after_B(s$choice, s$outcome, -1150)
    }
    out$wsls[i] <- mean(w, na.rm = TRUE)
    out$se[i] <- mean(s_e)
    out$bloss[i, ] <- colMeans(reps, na.rm = TRUE)
  }
  out
}
est <- list(
  orl = as.matrix(orl[, c("Arew", "Apun", "Kprime", "betaF", "betaP")]),
  variant = as.matrix(wl[, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")]),
  vse = as.matrix(vse[, c("theta", "Delta", "alpha", "phi", "beta")]),
  vse_split = as.matrix(split[, c("theta", "Delta", "alpha", "phi", "beta", "psi")]),
  eef = as.matrix(eef[, c("theta", "lambda", "phi", "beta")]))
sims <- list()
set.seed(2026); sims$orl <- sim_stats(sim_orl, est$orl, deck_outcomes)
set.seed(2026); sims$variant <- sim_stats(sim_orl_wl, est$variant, deck_outcomes)
set.seed(2026); sims$vse <- sim_stats(sim_vse, est$vse, deck_draws)
set.seed(2026); sims$vse_split <- sim_stats(sim_vse_split, est$vse_split, deck_draws)
set.seed(2026); sims$eef <- sim_stats(sim_eef, est$eef, deck_draws)

# 3. VSE under-produces WSLS (simulated < observed, paired t, p < .05);
#    VSE-split reproduces it (95% CI of the paired difference includes 0)
wsls_obs <- sapply(1:n, function(i) wsls_of(choice_h[i, ], net_h[i, ]))
round(c(observed = mean(wsls_obs), sapply(sims, function(s) mean(s$wsls))), 3)
p3_vse <- t.test(wsls_obs, sims$vse$wsls, paired = TRUE)
p3_split <- t.test(wsls_obs, sims$vse_split$wsls, paired = TRUE)
p3_vse
p3_split
pred3 <- p3_vse$p.value < .05 & mean(wsls_obs - sims$vse$wsls) > 0 &
  p3_split$conf.int[1] <= 0 & p3_split$conf.int[2] >= 0

# 4. Deck B: VSE's gap (observed minus simulated P(choose B), trials 2-10 after the big loss),
#    in absolute value, at most half the variant's gap from script 22
b <- readRDS("output/deckB_after_loss.rds")
obs <- rowMeans(b$obs_loss[1:n, 2:10])
gap_variant <- mean(obs - rowMeans(b$wl_loss[1:n, 2:10]), na.rm = TRUE)
gaps <- sapply(sims, function(s) mean(obs - rowMeans(s$bloss[, 2:10]), na.rm = TRUE))
round(c(variant_script22 = gap_variant, gaps), 3)
round(rbind(observed = colMeans(b$obs_loss[1:n, ], na.rm = TRUE),
            t(sapply(sims, function(s) colMeans(s$bloss, na.rm = TRUE))))[, c(1, 2, 3, 5, 10, 20)], 2)
t.test(obs, rowMeans(sims$vse$bloss[, 2:10]), paired = TRUE)
pred4 <- abs(gaps["vse"]) <= gap_variant / 2

# 5. SE index: VSE and VSE-split closer to the observed than the variant (mean |observed - simulated|)
se_obs <- sapply(1:n, function(i) se_index(choice_h[i, ]))
round(c(observed = mean(se_obs), sapply(sims, function(s) mean(s$se))), 3)
se_err <- sapply(sims, function(s) mean(abs(se_obs - s$se)))
round(se_err, 3)
pred5 <- se_err["vse"] < se_err["variant"] & se_err["vse_split"] < se_err["variant"]

# 6. psi recovers at r >= .70 (script 33, data simulated from VSE-split)
rec_file <- sprintf("output/vse/recovery_vse_split%s.rds", tag)
pred6 <- NA
if (file.exists(rec_file)) {
  rec <- readRDS(rec_file)
  r_split <- diag(cor(rec$true, rec$rec_split))
  print(round(r_split, 2))
  pred6 <- r_split["psi"] >= .70
}

# Not predicted, reported either way
round(c(variant_minus_vse_AIC = aic["variant"] - aic["vse"], variant_minus_vse_BIC = bic["variant"] - bic["vse"],
        eef_minus_vse_AIC = aic["eef"] - aic["vse"], eef_minus_vse_BIC = bic["eef"] - bic["vse"]))
sum(aic_i[, "variant"] < aic_i[, "vse"])   # people better fit by the variant than VSE (AIC)
rec0_file <- sprintf("output/vse/recovery_vse%s.rds", tag)
if (file.exists(rec0_file)) {           # how often the split is wrongly picked when VSE is true
  rec0 <- readRDS(rec0_file)
  print(c(aic = mean(rec0$nll_vse - rec0$nll_split > 1), bic = mean(rec0$nll_vse - rec0$nll_split > log(100) / 2)))
}

# Scorecard
data.frame(prediction = 1:6, met = c(pred1, pred2, pred3, pred4, pred5, pred6), row.names = NULL)

saveRDS(list(table = tab, aic = aic, bic = bic, sims = sims,
             observed = list(wsls = wsls_obs, se = se_obs, bloss = b$obs_loss[1:n, ]),
             gaps = c(variant_script22 = gap_variant, gaps), se_err = se_err,
             met = c(pred1, pred2, pred3, pred4, pred5, pred6)),
        sprintf("output/vse_checks%s.rds", tag))
