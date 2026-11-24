# Exploratory: the models fit with a group prior ("pooled"), in the same ten-fold cross-validation as every other
# fit (scripts 35, 45 and 50b: fold f holds out trials f, f + 10, ..., f + 90; held-out outcomes still update the
# model, but only held-out choices are scored).
# Method: empirical Bayes by expectation-maximization with a Laplace approximation, the standard shortcut for
# hierarchical fits in computational psychiatry (Huys et al.; the emfit toolbox). Each parameter is mapped to an
# unbounded scale, x = logit of its position within its limits, and given a normal group prior N(mu, sigma^2):
#   E-step: each person's fit maximizes the log-likelihood of their training trials + log N(x; mu, sigma^2), and the
#           Hessian at that point gives the uncertainty v of each parameter (Laplace approximation);
#   M-step: mu = mean(x), sigma^2 = mean(x^2 + v) - mu^2, over participants;
# repeated until mu and sigma change by less than 0.01 (at most 30 rounds). People whose choices pin a parameter down
# keep values near their own fit; people whose choices barely constrain it are pulled toward the group.
# The prior is estimated from training trials only (every fit in it is a fit to training trials), so held-out choices
# never inform it. Starting points: the unpooled fit for the same fold (scripts 35, 45, 50b), moved just inside the
# limits. WSLS needs no fitting and is not refit.
# Usage (from the project folder): Rscript scripts/54_pooled_cv.R <data> <model> <fold> [quick]
#   data: horstmann, replication or replication438 (without the two repeated records); model: orl, variant, orl_alpha, variant_alpha, vse, vse_split, eef, nova
#   quick: the first 12 participants and 3 rounds, for testing
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
model <- args[2]
fold <- as.integer(args[3])
quick <- "quick" %in% args

source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")                      # also loads script 43b (the eight models)
src <- cv_source(data_name)                              # replication438 starts from the replication fits
d <- load_set(src)
m <- if (model == "nova") hybrid_models$hybrid else models[[model]]
stopifnot(!is.null(m))
n <- nrow(d$choice)
who <- which(keep_of(data_name, d$study))               # rows of the source data set (all, except replication438)
if (quick) who <- head(who, 12)

# The unpooled fits for this fold
unpooled_files <- if (model == "nova") {
  Sys.glob(sprintf("output/hybrid/cv/%s/hybrid_fold%02d_part*.rds", src, fold))
} else if (src == "horstmann") {
  sprintf("output/cv_mle/%s_fold%02d.rds", model, fold)
} else Sys.glob(sprintf("output/prereg/cv/%s/%s_fold%02d_part*.rds", src, model, fold))
unpooled_files <- unpooled_files[!grepl("_quick", unpooled_files)]
start <- do.call(rbind, lapply(unpooled_files, function(x) { r <- readRDS(x); if (is.data.frame(r)) r else r$fit }))
start <- start[order(start$subj), ]
stopifnot(identical(as.integer(start$subj), 1:n))

# Unbounded scale: x = logit((p - lower) / (upper - lower))
lo <- m$lower; hi <- m$upper; K <- length(lo)
to_x <- function(p) qlogis(pmin(pmax((p - lo) / (hi - lo), 1e-4), 1 - 1e-4))
to_p <- function(x) lo + (hi - lo) * plogis(x)

test <- test_trials(fold)
train <- !(1:100 %in% test)
nll <- function(x, i) -sum(m$lp(to_p(x), d$choice[i, ], d$win[i, ], d$loss[i, ])[train])

X <- t(sapply(who, function(i) to_x(unlist(start[i, m$names]))))
mu <- colMeans(X)
sig2 <- pmax(apply(X, 2, var), 1e-4)
V <- matrix(NA, length(who), K)
conv <- rep(NA, length(who))
rounds <- if (quick) 3 else 30
for (r in 1:rounds) {
  for (k in seq_along(who)) {                           # E-step: each person's best fit under the group prior
    i <- who[k]
    obj <- function(x) nll(x, i) + sum((x - mu)^2 / (2 * sig2))
    f <- optim(X[k, ], obj, method = "BFGS", control = list(maxit = 500))
    X[k, ] <- f$par
    conv[k] <- f$convergence
    H <- tryCatch(optimHess(f$par, obj), error = function(e) diag(1 / sig2))
    v <- tryCatch(diag(solve(H)), error = function(e) 1 / pmax(diag(H), 1e-6))
    V[k, ] <- pmin(pmax(v, 0), 100)
  }
  mu_new <- colMeans(X)                                 # M-step: the group prior
  sig2_new <- pmax(colMeans(X^2 + V) - mu_new^2, 1e-4)
  change <- max(abs(mu_new - mu), abs(sqrt(sig2_new) - sqrt(sig2)))
  mu <- mu_new; sig2 <- sig2_new
  cat(sprintf("round %d: largest change %.4f\n", r, change))
  if (change < 0.01) break
}

# Held-out predictions with the pooled fits
res <- data.frame(subj = who)
res[c(m$names, "train_nll", "heldout_ll", "conv")] <- NA
heldout_lp <- matrix(NA, length(who), length(test), dimnames = list(NULL, test))
for (k in seq_along(who)) {
  i <- who[k]
  p <- to_p(X[k, ])
  lp <- m$lp(p, d$choice[i, ], d$win[i, ], d$loss[i, ])
  res[k, m$names] <- p
  res$train_nll[k] <- -sum(lp[train])
  res$heldout_ll[k] <- sum(lp[!train])
  res$conv[k] <- conv[k]
  heldout_lp[k, ] <- lp[test]
}
out_dir <- file.path("output/pooled/cv", data_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(fit = res, heldout_lp = heldout_lp, mu = setNames(mu, m$names), sigma = setNames(sqrt(sig2), m$names),
             rounds = r, last_change = change),
        sprintf("%s/%s_fold%02d%s.rds", out_dir, model, fold, if (quick) "_quick" else ""))
table(res$conv)                                         # optimizer codes (0 = converged)
round(rbind(group_mean = to_p(mu), sd_on_logit_scale = sqrt(sig2)), 3)
