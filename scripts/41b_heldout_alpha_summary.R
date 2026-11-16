# Scores the held-out test of outcome curvature (script 41) prediction by prediction,
# exactly as worded in the README (written November 16, 2026, before script 41 was run).
quick <- FALSE
tag <- if (quick) "_quick" else ""
orl_ho <- readRDS("output/orl_fit_heldout.rds")
wl_ho <- readRDS("output/orl_wl_fit_heldout.rds")
fits <- readRDS(sprintf("output/heldout_alpha/fits%s.rds", tag))
n <- nrow(fits$orl_alpha)

# Summed AIC and BIC (k = number of free parameters, 100 trials per participant)
k <- c(orl = 5, variant = 6, orl_alpha = 6, variant_alpha = 7, vse = 5)
nll <- cbind(orl = orl_ho$nll[1:n], variant = wl_ho$nll[1:n], orl_alpha = fits$orl_alpha$nll,
             variant_alpha = fits$variant_alpha$nll, vse = fits$vse$nll)
aic <- colSums(sweep(2 * nll, 2, 2 * k, "+"))
bic <- colSums(sweep(2 * nll, 2, k * log(100), "+"))
round(rbind(AIC = aic, BIC = bic), 1)
# Participants each model fits best (reported, not predicted)
models_fit <- colnames(nll)
rbind(AIC = table(factor(models_fit[apply(sweep(2 * nll, 2, 2 * k, "+"), 1, which.min)], levels = models_fit)),
      BIC = table(factor(models_fit[apply(sweep(2 * nll, 2, k * log(100), "+"), 1, which.min)], levels = models_fit)))

# Prediction 1: summed AIC favors ORL + alpha over ORL
pred1 <- aic[["orl_alpha"]] < aic[["orl"]]

# Prediction 2: alpha < 1 for more than 55% of participants (one-sided binomial test against .55,
# the share below 1 when data were simulated with alpha = 1, script 29b; p < .05), and median alpha < 1
a <- fits$orl_alpha$alpha
below <- sum(a < 1)
p2 <- binom.test(below, n, p = .55, alternative = "greater")$p.value
c(below = below, share = below / n, p = p2, median = median(a), at_lower_limit = sum(a <= .0101))
pred2 <- p2 < .05 & median(a) < 1

# Cross-validation: held-out log-likelihood per participant, summed over the 10 folds
models <- c("orl", "variant", "orl_alpha", "variant_alpha")
ll <- sapply(models, function(m) {
  Reduce(`+`, lapply(1:10, function(f)
    readRDS(sprintf("output/heldout_alpha/cv_%s_fold%02d%s.rds", m, f, tag))$heldout_ll))
})
round(colSums(ll), 1)
pair <- function(a, b) {   # positive = first model predicts better
  d <- ll[, a] - ll[, b]
  c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), first_better_for = sum(d > 0))
}
pairs <- rbind(`orl_alpha - orl` = pair("orl_alpha", "orl"),
               `variant_alpha - orl_alpha` = pair("variant_alpha", "orl_alpha"),
               `variant - orl` = pair("variant", "orl"),
               `orl_alpha - variant` = pair("orl_alpha", "variant"))
round(pairs, 2)

# Prediction 3: ORL + alpha predicts held-out choices better than ORL (one-sided, z > 1.645)
pred3 <- pairs["orl_alpha - orl", "z"] > qnorm(.95)
# Prediction 4: once alpha is included, the split adds nothing (variant + alpha and ORL + alpha within 2 SE)
pred4 <- abs(pairs["variant_alpha - orl_alpha", "z"]) < 2
# Prediction 5: summed AIC favors ORL + alpha over VSE
pred5 <- aic[["orl_alpha"]] < aic[["vse"]]

met <- c(pred1, pred2, pred3, pred4, pred5)
data.frame(prediction = 1:5, met = met, row.names = NULL)

saveRDS(list(nll = nll, aic = aic, bic = bic, alpha = a, ll = ll, pairs = pairs, met = met),
        sprintf("output/heldout_alpha_summary%s.rds", tag))
