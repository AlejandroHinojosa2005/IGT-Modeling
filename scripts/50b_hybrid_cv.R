# Exploratory: held-out predictions of the hybrid model (script 50), with the same cross-validation as every other
# model (scripts 35 and 45). Ten folds within each person: fold f holds out trials f, f + 10, ..., f + 90 (trial 1 is
# never held out); each version is fit by maximum likelihood to the other trials (L-BFGS-B, the limits of script 50);
# held-out outcomes still update the model, but only the held-out choices are scored.
# Starts: 20 random starts (seed 2050 + f, drawn participant by participant before any fitting, so splitting the
# participants into parts does not change any start), plus one more: the person's ORL + alpha fit for the same fold
# (script 35 for Horstmann, script 45 otherwise), with Q0 = 0 and D = 1. At that point the net-bending versions are
# exactly ORL + alpha, so their training fit can never be worse than ORL + alpha's.
# Usage (from the project folder):
#   Rscript scripts/50b_hybrid_cv.R <data> <version> <fold> <part> <parts> [quick]
#     data: horstmann or replication
#     version: orl_alpha_start, orl_alpha_decay, hybrid, hybrid_sep; fold: 1-10
#     part, parts: fit part <part> of <parts> equal blocks of participants
#     quick: the first 2 participants of the part only, for testing
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
version <- args[2]
fold <- as.integer(args[3])
part <- as.integer(args[4])
parts <- as.integer(args[5])
quick <- "quick" %in% args

source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
d <- load_set(data_name)
m <- hybrid_models[[version]]
stopifnot(!is.null(m))
n <- nrow(d$choice)

# All random starts for this fold, participant by participant
set.seed(2050 + fold)
starts <- lapply(1:n, function(i) lapply(1:20, function(s) runif(length(m$lower), m$lower, m$upper)))

# The ORL + alpha fits for this fold, one row per participant in order
oa_files <- if (data_name == "horstmann") sprintf("output/cv_mle/orl_alpha_fold%02d.rds", fold) else
  Sys.glob(sprintf("output/prereg/cv/%s/orl_alpha_fold%02d_part*.rds", data_name, fold))
oa_files <- oa_files[!grepl("_quick", oa_files)]
oa <- do.call(rbind, lapply(oa_files, function(x) { r <- readRDS(x); if (is.data.frame(r)) r else r$fit }))
oa <- oa[order(oa$subj), ]
stopifnot(identical(as.integer(oa$subj), 1:n))
nested <- function(i) {                                  # ORL + alpha's parameters, then Q0 = 0 and/or D = 1
  p <- unlist(oa[i, models$orl_alpha$names])
  extra <- c(Q0 = 0, D = 1)[intersect(c("Q0", "D"), m$names)]
  pmin(pmax(c(p, extra), m$lower), m$upper)
}

test <- test_trials(fold)
train <- !(1:100 %in% test)
who <- if (parts == 1) 1:n else split(1:n, cut(1:n, parts, labels = FALSE))[[part]]
if (quick) who <- head(who, 2)
res <- data.frame(subj = who)
res[c(m$names, "train_nll", "heldout_ll", "conv", "oa_train_nll", "best_start")] <- NA
heldout_lp <- matrix(NA, length(who), length(test), dimnames = list(NULL, test))
for (r in seq_along(who)) {
  i <- who[r]
  nll_train <- function(p) -sum(m$lp(p, d$choice[i, ], d$win[i, ], d$loss[i, ])[train])
  best <- list(value = Inf)
  all_starts <- c(starts[[i]], list(nested(i)))
  for (s in seq_along(all_starts)) {
    f <- optim(all_starts[[s]], nll_train, method = "L-BFGS-B", lower = m$lower, upper = m$upper, control = list(maxit = 1000))
    if (f$value < best$value) { best <- f; best$start <- s }
  }
  lp <- m$lp(best$par, d$choice[i, ], d$win[i, ], d$loss[i, ])
  res[r, m$names] <- best$par
  res$train_nll[r] <- best$value
  res$heldout_ll[r] <- sum(lp[!train])
  res$conv[r] <- best$convergence
  res$oa_train_nll[r] <- oa$train_nll[i]
  res$best_start[r] <- best$start                       # 21 = the ORL + alpha start
  heldout_lp[r, ] <- lp[test]
}

out_dir <- file.path("output/hybrid/cv", data_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(fit = res, heldout_lp = heldout_lp),
        sprintf("%s/%s_fold%02d_part%d%s.rds", out_dir, version, fold, part, if (quick) "_quick" else ""))
table(res$conv)                                         # optimizer codes (0 = converged)
if (version != "hybrid_sep") max(res$train_nll - res$oa_train_nll)   # should be <= 0 (never worse than ORL + alpha)
