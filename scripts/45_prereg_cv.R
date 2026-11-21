# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), step 2: held-out predictions.
# Ten-fold cross-validation within each person, exactly as script 35: fold f holds out trials f, f + 10, ..., f + 90
# (trial 1 is never held out); each model is fit by maximum likelihood to the other trials (L-BFGS-B, 20 random starts,
# seed 2026 + f, the limits of the main fits); held-out outcomes still update the model, but only the held-out choices
# are scored. WSLS needs no fitting; script 46 scores it.
# The participants can be split into parts, to run as separate cluster jobs. All random starts are drawn first, in the
# same order as script 35 (participant by participant, 20 starts each), so the split does not change any start, and
# running this on "horstmann" reproduces script 35's fits exactly.
# Nothing here is scored or summarized: the held-out results are read only by the preregistered scripts 46-48.
# Usage (from the project folder):
#   Rscript scripts/45_prereg_cv.R <data> <model> <fold> <part> <parts> [quick]
#     data: replication, simulated (or horstmann, to check against script 35)
#     model: orl, variant, orl_alpha, variant_alpha, vse, vse_split, eef; fold: 1-10
#     part, parts: fit part <part> of <parts> equal blocks of participants
#     quick: the first 2 participants of the part only, for testing
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
model <- args[2]
fold <- as.integer(args[3])
part <- as.integer(args[4])
parts <- as.integer(args[5])
quick <- "quick" %in% args

source("scripts/43_prereg_data.R")
source("scripts/43b_prereg_models.R")
d <- load_set(data_name)
m <- models[[model]]
n <- nrow(d$choice)

# All random starts for this fold, in script 35's order
set.seed(2026 + fold)
starts <- lapply(1:n, function(i) lapply(1:20, function(s) runif(length(m$lower), m$lower, m$upper)))

test <- test_trials(fold)
train <- !(1:100 %in% test)
who <- split(1:n, cut(1:n, parts, labels = FALSE))[[part]]
if (quick) who <- head(who, 2)
res <- data.frame(subj = who)
res[c(m$names, "train_nll", "heldout_ll", "conv")] <- NA
heldout_lp <- matrix(NA, length(who), length(test), dimnames = list(NULL, test))
for (r in seq_along(who)) {
  i <- who[r]
  nll_train <- function(p) -sum(m$lp(p, d$choice[i, ], d$win[i, ], d$loss[i, ])[train])
  best <- list(value = Inf)
  for (st in starts[[i]]) {
    f <- optim(st, nll_train, method = "L-BFGS-B", lower = m$lower, upper = m$upper, control = list(maxit = 1000))
    if (f$value < best$value) best <- f
  }
  lp <- m$lp(best$par, d$choice[i, ], d$win[i, ], d$loss[i, ])
  res[r, m$names] <- best$par
  res$train_nll[r] <- best$value
  res$heldout_ll[r] <- sum(lp[!train])
  res$conv[r] <- best$convergence
  heldout_lp[r, ] <- lp[test]
}

out_dir <- file.path("output/prereg/cv", data_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(fit = res, heldout_lp = heldout_lp),
        sprintf("%s/%s_fold%02d_part%d%s.rds", out_dir, model, fold, part, if (quick) "_quick" else ""))
table(res$conv)                                         # optimizer codes (0 = converged)
