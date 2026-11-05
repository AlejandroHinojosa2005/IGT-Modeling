# Exploratory (not predicted): what happens after deck B's big loss (-1,150)?
# For each participant, P(choose deck B) on each of the 20 trials after a big B loss,
# compared with after an ordinary B win (+100). Observed vs. simulated from the fitted ORL
# and variant (same logic as the WSLS check). Deeper avoidance and a sharper return than the
# models produce would hint at latent-state inference ("deck B has changed").
source("scripts/10_orl_variant.R")
orl <- readRDS("output/orl_fit_horstmann.rds")
wl <- readRDS("output/orl_wl_fit_horstmann.rds")
deck_outcomes <- list(net_h[choice_h == 1], net_h[choice_h == 2], net_h[choice_h == 3], net_h[choice_h == 4])
lags <- 20
n_subj <- nrow(choice_h)

# P(choose B) on each of the next 20 trials after a deck B outcome equal to `target`
after_B <- function(choice, outcome, target) {
  events <- which(choice == 2 & outcome == target)
  events <- events[events <= length(choice) - lags]   # keep events with 20 trials after them
  if (length(events) == 0) return(rep(NA, lags))
  back <- matrix(NA, length(events), lags)
  for (e in 1:length(events)) back[e, ] <- choice[events[e] + 1:lags] == 2
  colMeans(back)
}

# Observed: one curve per person
obs_loss <- matrix(NA, n_subj, lags)
obs_win <- matrix(NA, n_subj, lags)
for (i in 1:n_subj) {
  obs_loss[i, ] <- after_B(choice_h[i, ], net_h[i, ], -1150)
  obs_win[i, ] <- after_B(choice_h[i, ], net_h[i, ], 100)
}

# Simulated: 20 datasets per person from their fitted parameters
sim_curves <- function(sim_fun, est, target) {
  curves <- matrix(NA, n_subj, lags)
  for (i in 1:n_subj) {
    reps <- matrix(NA, 20, lags)
    for (r in 1:20) {
      s <- sim_fun(est[i, ], deck_outcomes)
      reps[r, ] <- after_B(s$choice, s$outcome, target)
    }
    curves[i, ] <- colMeans(reps, na.rm = TRUE)
  }
  curves
}
est_orl <- as.matrix(orl[, c("Arew", "Apun", "Kprime", "betaF", "betaP")])
est_wl <- as.matrix(wl[, c("Arew", "Apun", "Kprime", "betaF", "betaPw", "betaPl")])
set.seed(2026)
orl_loss <- sim_curves(sim_orl, est_orl, -1150)
set.seed(2026)
wl_loss <- sim_curves(sim_orl_wl, est_wl, -1150)

# Average across people
curves <- rbind(observed = colMeans(obs_loss, na.rm = TRUE),
                ORL = colMeans(orl_loss, na.rm = TRUE),
                variant = colMeans(wl_loss, na.rm = TRUE),
                observed_after_win = colMeans(obs_win, na.rm = TRUE))
round(curves[, c(1, 2, 3, 5, 10, 20)], 2)   # selected lags
sum(!is.na(obs_loss[, 1]))                  # people with at least one usable big loss
saveRDS(list(obs_loss = obs_loss, obs_win = obs_win, orl_loss = orl_loss, wl_loss = wl_loss),
        "output/deckB_after_loss.rds")

# Figure
png("output/figures/deckB_after_bigloss.png", width = 1600, height = 1000, res = 200)
matplot(1:lags, t(curves), type = "l", lty = 1, lwd = 2,
        col = c("black", "red", "blue", "gray50"),
        xlab = "Trials after the deck B outcome", ylab = "P(choose deck B)",
        main = "Returning to deck B after its big loss")
legend("topright", lty = 1, lwd = 2, bty = "n", col = c("black", "red", "blue", "gray50"),
       legend = c("Observed, after big loss", "ORL simulated", "Variant simulated", "Observed, after a B win"))
dev.off()


# Paired tests: observed vs. simulated P(choose B), averaged over trials 2-10 after the big loss
obs <- rowMeans(obs_loss[, 2:10])
t.test(obs, rowMeans(orl_loss[, 2:10]), paired = TRUE)
t.test(obs, rowMeans(wl_loss[, 2:10]), paired = TRUE)
mean(obs > rowMeans(orl_loss[, 2:10]))   # share of people above ORL
mean(obs > rowMeans(wl_loss[, 2:10]))    # share above the variant







