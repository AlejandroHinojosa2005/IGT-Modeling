# Exploratory: which models can the shared mix do without? For later work (EEG) that needs fewer models.
# Uses script 46's held-out predictions (no refitting). The shared mix is refit with some models left out, with the
# same procedure as script 46 (weights fit on the other nine folds of all participants, scored on the held-out fold),
# and compared with the mix of all eight: summed over participants, SE = sd of the per-participant differences
# x sqrt(n), as in script 48.
#   drop_<model>: all eight but one.
#   five_nested_out: without ORL, ORL + alpha and VSE. Each is a special case of a model that stays (ORL + split,
#         ORL + split + alpha and VSE-split: equal perseverance after wins and losses, or psi = 0), so these five keep
#         every mechanism of the eight.
#   five_cheapest: without ORL, ORL + split and VSE, the three whose removal one at a time cost least in both samples
#         (chosen after a first run of the drop_<model> comparisons).
# Each five is also compared with the best single model, ORL + alpha.
# Usage (from the project folder): Rscript scripts/52_model_subsets.R <data>   (horstmann or replication)
library(loo)
data_name <- commandArgs(trailingOnly = TRUE)[1]
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
lpd <- s$lpd
rows <- s$rows
n <- max(rows$subj)

mix_score <- function(lp_rows, w) {                     # as script 46
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
all8 <- shared_cv(lpd)                                  # the shared mix of all eight, refit here like the others
max(abs(all8 - s$cv_group))                             # against script 46's (about 0 with the same loo version)
sets <- list()
for (m in colnames(lpd)) sets[[paste0("drop_", m)]] <- setdiff(colnames(lpd), m)
sets$five_nested_out <- setdiff(colnames(lpd), c("orl", "orl_alpha", "vse"))
sets$five_cheapest <- setdiff(colnames(lpd), c("orl", "variant", "vse"))
scores <- sapply(sets, function(k) shared_cv(lpd[, k]))
per <- rowsum(scores - all8, rows$subj)
subsets <- round(t(apply(per, 2, function(x) c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n))))), 1)
subsets
fives <- c("five_nested_out", "five_cheapest")
w_five <- lapply(sets[fives], function(k) round(setNames(as.numeric(stacking_weights(lpd[, k])), k), 3))
w_five                                                  # weights of each five, from all held-out choices
five_vs_orl_alpha <- round(t(sapply(fives, function(v) {  # against the best single model (as H1 of the registration)
  x <- rowsum(scores[, v] - lpd[, "orl_alpha"], rows$subj)
  c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n)), better_for = sum(x > 0))
})), 1)
five_vs_orl_alpha

saveRDS(list(data = data_name, sets = sets, subsets = subsets, w_five = w_five, five_vs_orl_alpha = five_vs_orl_alpha,
             per = per), sprintf("output/prereg/subsets_%s.rds", data_name))
