# Exploratory: scores the pooled fits of script 54 against the unpooled fits, on the same held-out choices.
# Differences are summed over participants; SE = sd of the per-participant differences x sqrt(n) (as scripts 48, 50c).
#   1. Each model pooled against the same model unpooled.
#   2. NOVA against ORL + alpha, both pooled.
#   3. Shared mixes (as script 46): the eight with pooled members (WSLS unchanged) against the eight unpooled; five with
#      NOVA (NOVA, ORL + split + alpha, VSE-split, EEF, WSLS), pooled, against the eight unpooled.
#   4. The group means and spreads of each model's parameters (fold 1).
# Needs script 46's stacking_<data>.rds and NOVA's unpooled predictions (script 50c for Horstmann, 50d for replication).
# Usage (from the project folder): Rscript scripts/54b_pooled_score.R <data>   (horstmann, replication or replication438)
library(loo)
data_name <- commandArgs(trailingOnly = TRUE)[1]
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
d <- load_kept(data_name)
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
rows <- s$rows
n <- max(rows$subj)
stopifnot(n == nrow(d$choice), nrow(rows) == n * 99)
row_of <- function(i, trials) (i - 1) * 99 + trials - 1
nova_unpooled <- if (data_name == "horstmann") readRDS("output/hybrid/score_horstmann.rds")$lpd_h[, "hybrid"] else
  readRDS(sprintf("output/hybrid/nova_%s.rds", data_name))$nova

# Pooled held-out predictions, in the rows of script 46
fitted <- c(names(models), "nova")
pooled <- matrix(NA, nrow(rows), length(fitted), dimnames = list(NULL, fitted))
check <- 0
groups <- list()
for (mn in fitted) {
  m <- if (mn == "nova") hybrid_models$hybrid else models[[mn]]
  for (f in 1:10) {
    r <- readRDS(sprintf("output/pooled/cv/%s/%s_fold%02d.rds", data_name, mn, f))
    fit <- r$fit[order(r$fit$subj), ]
    stopifnot(identical(as.integer(fit$subj), d$kept))
    if (f == 1) groups[[mn]] <- rbind(group_mean = m$lower + (m$upper - m$lower) * plogis(r$mu), sd_logit = r$sigma)
    for (i in 1:n) {
      lp <- m$lp(unlist(fit[i, m$names]), d$choice[i, ], d$win[i, ], d$loss[i, ])
      k <- row_of(i, test_trials(f))
      pooled[k, mn] <- lp[rows$trial[k]]
      check <- max(check, abs(sum(lp[rows$trial[k]]) - fit$heldout_ll[i]))
    }
  }
}
check                                                   # refit predictions match the saved totals (should be about 0)
stopifnot(!anyNA(pooled))

unpooled <- cbind(s$lpd[, names(models)], nova = nova_unpooled)
both <- cbind(pooled, unpooled)
colnames(both) <- c(paste0("pooled.", colnames(pooled)), paste0("unpooled.", colnames(unpooled)))
per <- rowsum(both, rows$subj)
compare <- function(a, b) { x <- per[, a] - per[, b]; c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n)), better_for = sum(x > 0)) }

# 1. Pooled against unpooled, model by model
pooling_gain <- round(t(sapply(fitted, function(mn) compare(paste0("pooled.", mn), paste0("unpooled.", mn)))), 1)
pooling_gain
# 2. NOVA against ORL + alpha, both pooled
round(compare("pooled.nova", "pooled.orl_alpha"), 1)

# 3. Mixes
mix_score <- function(lp_rows, w) {
  ok <- w > 0
  a <- sweep(lp_rows[, ok, drop = FALSE], 2, log(w[ok]), "+")
  mx <- apply(a, 1, max)
  mx + log(rowSums(exp(a - mx)))
}
shared_cv <- function(L) {
  x <- rep(NA, nrow(L))
  for (f in 1:10) x[rows$fold == f] <- mix_score(L[rows$fold == f, , drop = FALSE], as.numeric(stacking_weights(L[rows$fold != f, , drop = FALSE])))
  x
}
L <- cbind(pooled, wsls = s$lpd[, "wsls"])
mixes <- cbind(eight_unpooled = s$cv_group, eight_pooled = shared_cv(L[, model_names]),
               five_nova_pooled = shared_cv(L[, c("nova", "variant_alpha", "vse_split", "eef", "wsls")]))
per <- cbind(per, rowsum(mixes, rows$subj))
mix_results <- round(rbind(eight_pooled_vs_unpooled = compare("eight_pooled", "eight_unpooled"),
                           five_nova_pooled_vs_eight_unpooled = compare("five_nova_pooled", "eight_unpooled")), 1)
mix_results

# 4. Group parameters (fold 1): group mean on the parameter's own scale, spread on the logit scale
lapply(groups, round, 3)

dir.create("output/pooled", showWarnings = FALSE)
saveRDS(list(data = data_name, pooled = pooled, per = per, pooling_gain = pooling_gain, mix_results = mix_results,
             groups = groups), sprintf("output/pooled/score_%s.rds", data_name))
