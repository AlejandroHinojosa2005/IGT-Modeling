# Summarizes the maximum-likelihood cross-validation (script 35) and scores the dated prediction.
# Higher held-out log-likelihood = better prediction of choices the model was not fit to.
# Differences are summed per participant first; SE = sd of the per-participant differences * sqrt(n),
# as in the Stan cross-validation (script 23b).
quick <- FALSE
tag <- if (quick) "_quick" else ""
source("scripts/01_load_data.R")
models <- c("orl", "variant", "orl_alpha", "variant_alpha", "vse", "vse_split", "eef")
n <- if (quick) 3 else nrow(choice_h)

# Held-out log-likelihood per participant, summed over the 10 folds
ll <- matrix(NA, n, length(models) + 1, dimnames = list(NULL, c(models, "wsls")))
missing <- c()
for (m in models) {
  tot <- rep(0, n)
  for (f in 1:10) {
    file <- sprintf("output/cv_mle/%s_fold%02d%s.rds", m, f, tag)
    if (!file.exists(file)) {
      missing <- c(missing, file)
      next
    }
    tot <- tot + readRDS(file)$heldout_ll
  }
  ll[, m] <- tot
}
missing   # should be empty (NULL)

# WSLS: stay/shift proportions from the training trials, held-out trials scored.
# Probabilities kept within .01-.99 so one event never seen in training isn't scored as impossible.
ll[, "wsls"] <- 0
for (f in 1:10) {
  test <- seq(f, 100, by = 10)
  test <- test[test != 1]
  train <- setdiff(2:100, test)
  for (i in 1:n) {
    ch <- choice_h[i, ]
    stay <- c(NA, ch[-1] == ch[-100])
    gain <- c(NA, net_h[i, -100] >= 0)
    p_sw <- if (any(gain[train])) mean(stay[train][gain[train]]) else 0.5
    p_sl <- if (any(!gain[train])) mean(!stay[train][!gain[train]]) else 0.5
    p_sw <- min(max(p_sw, .01), .99)
    p_sl <- min(max(p_sl, .01), .99)
    p <- ifelse(gain, ifelse(stay, p_sw, (1 - p_sw) / 3), ifelse(stay, 1 - p_sl, p_sl / 3))
    ll[i, "wsls"] <- ll[i, "wsls"] + sum(log(p[test]))
  }
}

# Table: total held-out log-likelihood, difference from the best model with its SE,
# and how many participants each model predicts best
tot <- colSums(ll)
best <- names(which.max(tot))
tab <- data.frame(heldout_ll = round(tot, 1),
                  diff_vs_best = round(tot - tot[best], 1),
                  se = round(apply(ll, 2, function(x) sd(x - ll[, best]) * sqrt(n)), 1),
                  best_for = tabulate(apply(ll, 1, which.max), ncol(ll)))
tab[order(-tab$heldout_ll), ]

# Key pairwise comparisons (positive = first model predicts better)
pair <- function(a, b) {
  d <- ll[, a] - ll[, b]
  c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), first_better_for = sum(d > 0))
}
pairs <- rbind(`vse_split - vse` = pair("vse_split", "vse"),
               `variant - orl` = pair("variant", "orl"),
               `variant_alpha - orl_alpha` = pair("variant_alpha", "orl_alpha"),
               `vse_split - variant` = pair("vse_split", "variant"),
               `vse_split - variant_alpha` = pair("vse_split", "variant_alpha"),
               `variant - vse` = pair("variant", "vse"),
               `eef - vse` = pair("eef", "vse"))
round(pairs, 1)

# Dated predictions (see README)
# 1. VSE-split has the highest total held-out log-likelihood of the eight models
# 2. VSE-split beats VSE by more than 2 SE
pred1 <- best == "vse_split"
pred2 <- pairs["vse_split - vse", "z"] > 2
data.frame(prediction = 1:2, met = c(pred1, pred2), row.names = NULL)

saveRDS(list(ll = ll, table = tab, pairs = pairs, met = c(pred1, pred2)),
        sprintf("output/cv_mle_summary%s.rds", tag))
