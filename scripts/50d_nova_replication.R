# Exploratory, with predictions committed before running (predictions/nova_replication.md): NOVA on the replication
# sample (N = 440). NOVA is the hybrid model of script 50 (version "hybrid": ORL + alpha with a starting value Q0 and
# decay D of unchosen decks, net bending). Its cross-validation fits come from script 50b, run on "replication" with
# the same folds and held-out choices as the preregistered analyses; the eight models' held-out predictions and their
# shared mix come from script 46 (output/prereg/stacking_replication.rds).
# The four tests (differences summed over participants; SE = sd of the per-participant differences x sqrt(n)):
#   P1  NOVA alone does not beat ORL + alpha                      met if z < 1.645
#   P2  adding NOVA to the eight improves the shared mix           met if z > 1.645
#   P3  five models with NOVA beat the shared mix of all eight     met if z > 1.645
#   P4  five models with NOVA beat the same five with ORL + alpha  met if z > 1.645
#   (one-sided tests at alpha = .05, the rule of the preregistered tests in script 48; no correction for multiple tests)
#   (the five: NOVA or ORL + alpha, ORL + split + alpha, VSE-split, EEF, WSLS; mixes fit as script 46's shared mix,
#    weights fit on the other nine folds)
# Also described, not tested: NOVA's weight in the mix of nine, the median per-person difference from ORL + alpha, the
# five largest losses, and the fitted Q0 and D.
# Usage (from the project folder): Rscript scripts/50d_nova_replication.R [data]
#   data: replication (the default, as committed) or replication438 (without the two repeated records; the same fits,
#         an exploratory check). Saves output/hybrid/nova_<data>.rds.
library(loo)
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
data_name <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(data_name)) data_name <- "replication"
d <- load_kept(data_name)
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
rows <- s$rows
n <- max(rows$subj)
stopifnot(n == length(d$kept), nrow(d$choice) == n, nrow(rows) == n * 99)
row_of <- function(i, trials) (i - 1) * 99 + trials - 1

# NOVA's held-out predictions, in the rows of script 46
m <- hybrid_models$hybrid
fits <- lapply(1:10, function(f) {
  files <- Sys.glob(sprintf("output/hybrid/cv/%s/hybrid_fold%02d_part*.rds", cv_source(data_name), f))
  files <- files[!grepl("_quick", files)]
  if (!length(files)) stop("no NOVA fits for fold ", f)
  fit <- do.call(rbind, lapply(files, function(x) readRDS(x)$fit))
  fit <- fit[order(fit$subj), ]
  fit <- fit[fit$subj %in% d$kept, ]                    # the participants this data set uses, in order
  stopifnot(identical(as.integer(fit$subj), d$kept))
  fit
})
nova <- rep(NA, nrow(rows))
check <- 0
for (f in 1:10) for (i in 1:n) {
  lp <- m$lp(unlist(fits[[f]][i, m$names]), d$choice[i, ], d$win[i, ], d$loss[i, ])
  k <- row_of(i, test_trials(f))
  nova[k] <- lp[rows$trial[k]]
  check <- max(check, abs(sum(lp[rows$trial[k]]) - fits[[f]]$heldout_ll[i]))
}
check                                                   # refit predictions match the saved totals (should be about 0)
stopifnot(!anyNA(nova))
max(sapply(fits, function(x) max(x$train_nll - x$oa_train_nll)))   # never worse than ORL + alpha in training (<= about 0)
mean(unlist(lapply(fits, function(x) x$best_start == 21)))         # share of fits where the ORL + alpha start won

# Mixes, as script 46's shared mix
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
L <- cbind(s$lpd, nova = nova)
five_nova <- c("nova", "variant_alpha", "vse_split", "eef", "wsls")
five_oa <- c("orl_alpha", "variant_alpha", "vse_split", "eef", "wsls")
mixes <- cbind(eight = s$cv_group, nine = shared_cv(L), five_nova = shared_cv(L[, five_nova]), five_oa = shared_cv(L[, five_oa]))
per <- rowsum(cbind(L, mixes), rows$subj)
compare <- function(a, b) { x <- per[, a] - per[, b]; c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n)), better_for = sum(x > 0)) }

tests <- rbind(P1_nova_vs_orl_alpha = compare("nova", "orl_alpha"),
               P2_nine_vs_eight = compare("nine", "eight"),
               P3_five_nova_vs_eight = compare("five_nova", "eight"),
               P4_five_nova_vs_five_orl_alpha = compare("five_nova", "five_oa"))
met <- c(tests[1, "z"] < 1.645, tests[2:4, "z"] > 1.645)
predictions <- data.frame(round(tests, 1), met = met)
predictions

# Described, not tested
w_nine <- round(setNames(as.numeric(stacking_weights(L)), colnames(L)), 3)
w_nine
x <- per[, "nova"] - per[, "orl_alpha"]
list(median_per_person = round(median(x), 2), five_largest_losses = unname(round(head(sort(x), 5), 1)))
par_avg <- function(p) rowMeans(sapply(fits, function(z) z[[p]]))
round(rbind(Q0 = summary(par_avg("Q0")), D = summary(par_avg("D"))), 3)

saveRDS(list(nova = nova, per = per, predictions = predictions, w_nine = w_nine), sprintf("output/hybrid/nova_%s.rds", data_name))
