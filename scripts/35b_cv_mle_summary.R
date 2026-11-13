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

# Sensitivity (added after scoring, not part of the predictions): a few participants were predicted
# catastrophically by the VSE models in single folds. Chance for a fold's ~10 held-out choices is about -14
# (10 * log(1/4)); ORL-family folds never fell below -66. VSE, VSE-split and EEF multiply their weights by a
# consistency term (up to 3^5 - 1 = 242), so a person whose training choices are almost perfectly predictable
# gets a near-certain model, and one unexpected held-out choice can score in the thousands.
# Same comparison without the participants who had any fold below -100 in any model (cutoff chosen after seeing
# the results), applied to every model alike.
fold_ll <- function(m) sapply(1:10, function(f)
  readRDS(sprintf("output/cv_mle/%s_fold%02d%s.rds", m, f, tag))$heldout_ll)
worst <- sapply(models, function(m) apply(fold_ll(m), 1, min))   # each participant's worst fold, per model
round(apply(worst, 2, min), 1)                                   # worst single fold anywhere, per model
bad <- which(apply(worst, 1, min) < -100)
bad                                                              # participants left out
ll_s <- ll[-bad, ]
n_s <- nrow(ll_s)
tot_s <- colSums(ll_s)
best_s <- names(which.max(tot_s))
tab_s <- data.frame(heldout_ll = round(tot_s, 1),
                    diff_vs_best = round(tot_s - tot_s[best_s], 1),
                    se = round(apply(ll_s, 2, function(x) sd(x - ll_s[, best_s]) * sqrt(n_s)), 1),
                    median_person = round(apply(ll_s, 2, median), 1))
tab_s[order(-tab_s$heldout_ll), ]
pair_s <- function(a, b) {
  d <- ll_s[, a] - ll_s[, b]
  c(diff = sum(d), se = sd(d) * sqrt(n_s), z = sum(d) / (sd(d) * sqrt(n_s)), first_better_for = sum(d > 0))
}
pairs_s <- rbind(`vse_split - vse` = pair_s("vse_split", "vse"),
                 `variant - orl` = pair_s("variant", "orl"),
                 `orl_alpha - orl` = pair_s("orl_alpha", "orl"),
                 `variant_alpha - orl_alpha` = pair_s("variant_alpha", "orl_alpha"),
                 `orl_alpha - vse` = pair_s("orl_alpha", "vse"),
                 `vse - orl` = pair_s("vse", "orl"),
                 `eef - vse` = pair_s("eef", "vse"))
round(pairs_s, 1)

saveRDS(list(left_out = bad, worst_fold = apply(worst, 2, min), table = tab_s, pairs = pairs_s),
        sprintf("output/cv_mle_sensitivity%s.rds", tag))
