# Scores the five predictions for the variant + alpha (written before fitting; see README).
# Needs scripts 28 (both models) and 29 to have run.
quick <- FALSE
tag <- if (quick) "_quick" else ""
source("scripts/27_alpha_model.R")
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2],
                      net_h[choice_h == 3], net_h[choice_h == 4])

wl <- readRDS("output/orl_wl_fit_horstmann.rds")              # variant, no alpha
orl <- readRDS("output/orl_fit_horstmann.rds")                # ORL, no alpha
wa <- readRDS(sprintf("output/alpha/alpha_wl%s.rds", tag))    # variant + alpha
oa <- readRDS(sprintf("output/alpha/alpha_orl%s.rds", tag))   # ORL + alpha
rec <- readRDS(sprintf("output/alpha/alpha_recovery%s.rds", tag))
n <- nrow(wa)
wl <- wl[1:n, ]
orl <- orl[1:n, ]

# 1. alpha < 1 for most participants (sign test) and median alpha < 1
summary(wa$alpha)
sum(wa$alpha < 1)
p1 <- binom.test(sum(wa$alpha < 1), n, alternative = "greater")
p1
pred1 <- p1$p.value < .05 & median(wa$alpha) < 1

# 2. The return to deck B: observed-minus-simulated gap over trials 2-10 after a big B loss,
#    at most half the variant's gap (same curves and method as script 22)
lags <- 20
after_B <- function(choice, outcome, target) {
  events <- which(choice == 2 & outcome == target)
  events <- events[events <= length(choice) - lags]
  if (length(events) == 0) return(rep(NA, lags))
  back <- matrix(NA, length(events), lags)
  for (e in 1:length(events)) back[e, ] <- choice[events[e] + 1:lags] == 2
  colMeans(back)
}
b <- readRDS("output/deckB_after_loss.rds")   # observed and variant curves from script 22
est_wa <- as.matrix(wa[, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl", "alpha")])
wa_loss <- matrix(NA, n, lags)
set.seed(2026)
for (i in 1:n) {
  reps <- matrix(NA, 20, lags)
  for (r in 1:20) {
    s <- sim_orl_wl_alpha(est_wa[i, ], deck_outcomes)
    reps[r, ] <- after_B(s$choice, s$outcome, -1150)
  }
  wa_loss[i, ] <- colMeans(reps, na.rm = TRUE)
}
obs <- rowMeans(b$obs_loss[1:n, 2:10])
gap_variant <- mean(obs - rowMeans(b$wl_loss[1:n, 2:10]), na.rm = TRUE)
gap_alpha <- mean(obs - rowMeans(wa_loss[, 2:10]), na.rm = TRUE)
c(gap_variant = gap_variant, gap_alpha = gap_alpha)
t.test(obs, rowMeans(wa_loss[, 2:10]), paired = TRUE)
round(rbind(observed = colMeans(b$obs_loss[1:n, ], na.rm = TRUE),
            variant = colMeans(b$wl_loss[1:n, ], na.rm = TRUE),
            variant_alpha = colMeans(wa_loss, na.rm = TRUE))[, c(1, 2, 3, 5, 10, 20)], 2)
pred2 <- gap_alpha <= gap_variant / 2

# 3. Summed AIC favors the variant + alpha over the variant (BIC reported either way)
aic <- c(variant = sum(2 * wl$nll + 2 * 6), variant_alpha = sum(2 * wa$nll + 2 * 7),
         orl = sum(2 * orl$nll + 2 * 5), orl_alpha = sum(2 * oa$nll + 2 * 6))
bic <- c(variant = sum(2 * wl$nll + 6 * log(100)), variant_alpha = sum(2 * wa$nll + 7 * log(100)),
         orl = sum(2 * orl$nll + 5 * log(100)), orl_alpha = sum(2 * oa$nll + 6 * log(100)))
round(aic)
round(bic)
sum(wl$nll - wa$nll > 1)               # people for whom alpha wins on AIC
sum(wl$nll - wa$nll > log(100) / 2)    # ... on BIC
pred3 <- aic["variant_alpha"] < aic["variant"]

# 4. alpha recovers (r >= .70); other parameters within .10 of the variant's earlier recovery
r_new <- diag(cor(rec$true, rec$recovered))
old <- readRDS("output/orl_wl_recovery.rds")
r_old <- diag(cor(old$true, old$recovered))
round(rbind(variant = c(r_old, alpha = NA), variant_alpha = r_new), 2)
pred4 <- r_new["alpha"] >= .70 & all(r_new[1:6] >= r_old - .10)

# 5. The variant's gains hold: betaPw > betaPl for at least 2/3 of people,
#    and simulated win-stay/lose-shift still matches the observed (paired t-test, p > .05)
sum(wa$betaPw > wa$betaPl)
wsls_of <- function(choice, outcome) {   # same as script 14
  m <- length(choice)
  stay <- choice[-1] == choice[-m]
  gain <- outcome[-m] >= 0
  mean(stay[gain]) - mean(stay[!gain])
}
wsls_obs <- rep(NA, n)
wsls_sim <- rep(NA, n)
set.seed(2026)
for (i in 1:n) {
  wsls_obs[i] <- wsls_of(choice_h[i, ], net_h[i, ])
  vals <- rep(NA, 20)
  for (r in 1:20) {
    s <- sim_orl_wl_alpha(est_wa[i, ], deck_outcomes)
    vals[r] <- wsls_of(s$choice, s$outcome)
  }
  wsls_sim[i] <- mean(vals, na.rm = TRUE)
}
p5 <- t.test(wsls_obs, wsls_sim, paired = TRUE)
p5
pred5 <- sum(wa$betaPw > wa$betaPl) >= ceiling(2 * n / 3) & p5$p.value > .05

# Scorecard
data.frame(prediction = 1:5, met = c(pred1, pred2, pred3, pred4, pred5), row.names = NULL)

saveRDS(list(alpha = summary(wa$alpha), gaps = c(variant = gap_variant, alpha = gap_alpha),
             wa_loss = wa_loss, aic = aic, bic = bic, recovery = rbind(variant = c(r_old, alpha = NA), variant_alpha = r_new),
             wsls = list(obs = wsls_obs, sim = wsls_sim),
             met = c(pred1, pred2, pred3, pred4, pred5)),
        sprintf("output/alpha_checks%s.rds", tag))
