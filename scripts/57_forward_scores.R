# Forward predictions, for the comparison with Centaur (centaur/score_centaur.py). Centaur predicts each choice only
# from the trials before it, and was trained on other people. Our other scores (scripts 45, 50b, 54) fit each person to
# 90% of their own trials, later ones included, and score the rest, which uses information Centaur never gets. Here
# each choice t (t = 2, ..., 100; trial 1 is never scored) is predicted from a fit to that person's trials 1 to t - 1
# only, under a group prior estimated from the OTHER sample (script 54's pooled fits, averaged over its ten folds:
# the Horstmann prior for replication participants, the replication prior for Horstmann participants). Nothing about
# the person being scored, other than their earlier trials, informs the prediction.
# The fit maximizes log-likelihood(trials 1 to t - 1) + log prior, as in script 54's E-step, starting from the fit for
# the previous trial (the prior mean at t = 2).
# WSLS: P(stay | gain) and P(shift | loss) from the earlier trials, with one added stay and one added shift (so 1/2
# before any data); a shift is spread evenly over the other three decks, as in script 43b.
# Usage (from the project folder):
#   Rscript scripts/57_forward_scores.R <data> <model> [part=P parts=N] [every=k] [who=file.csv] [prior=NAME] [quick]
#     data: horstmann or replication
#     model: orl, variant, orl_alpha, variant_alpha, vse, vse_split, eef, nova, wsls
#     part=P parts=N: run the P-th of N groups of participants (to run several at once)
#     every=k: refit every k trials instead of every trial (faster; predictions still use only earlier trials)
#     who=file.csv: only the participants in the file (columns sample and subj; rows with split = train are skipped),
#                   e.g. centaur/psych101_igt_match.csv, to score only the participants Centaur did not train on
#     prior=NAME: take the group prior from script 54's fits on NAME instead of the other sample (e.g.
#                 prior=replication438 for Horstmann); the output then goes to output/forward/<data>_prior-NAME/
#     quick: the first 3 participants, for testing (and the prior from whatever folds of the other sample exist)
# Output: output/forward/<data>/<model>_partP.rds: lp, the log-probability of each scored choice (participants x
# trials 2-100, rows named by participant), and each person's parameters after their last fit.
# Cost: 99 fits per person and model, about 5 seconds per person and model on one core (NOVA about 10).
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
model <- args[2]
opt <- function(key, default) {
  a <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(a)) sub("^[^=]*=", "", a[1]) else default
}
part <- as.integer(opt("part", 1))
parts <- as.integer(opt("parts", 1))
every <- as.integer(opt("every", 1))
who_file <- opt("who", "")
quick <- "quick" %in% args

source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")                       # also loads script 43b (the eight models)
d <- load_set(data_name)
n <- nrow(d$choice)
other <- opt("prior", if (data_name == "horstmann") "replication" else "horstmann")

who <- 1:n
if (who_file != "") {
  w <- read.csv(who_file)
  if ("split" %in% names(w)) w <- w[w$split != "train", ]
  who <- sort(unique(w$subj[w$sample == data_name]))
}
if (quick) who <- head(who, 3)
if (parts > 1) who <- split(who, cut(seq_along(who), parts, labels = FALSE))[[part]]

lp <- matrix(NA, length(who), 99, dimnames = list(who, 2:100))
par_last <- NULL

if (model == "wsls") {
  for (k in seq_along(who)) {
    i <- who[k]
    ch <- d$choice[i, ]
    stay <- c(NA, ch[-1] == ch[-100])
    gain <- c(NA, (d$win[i, -100] + d$loss[i, -100]) >= 0)   # gain[t]: was trial t - 1's net outcome >= 0?
    for (t in 2:100) {
      prev <- seq_len(t - 1)[-1]                          # the earlier transitions (trials 2 to t - 1)
      g <- prev[gain[prev]]
      s <- prev[!gain[prev]]
      p_sw <- (sum(stay[g]) + 1) / (length(g) + 2)          # P(stay | gain)
      p_sl <- (sum(!stay[s]) + 1) / (length(s) + 2)         # P(shift | loss)
      p <- if (gain[t]) { if (stay[t]) p_sw else (1 - p_sw) / 3 } else { if (stay[t]) 1 - p_sl else p_sl / 3 }
      lp[k, t - 1] <- log(p)
    }
  }
} else {
  m <- if (model == "nova") hybrid_models$hybrid else models[[model]]
  stopifnot(!is.null(m))
  lo <- m$lower; hi <- m$upper
  to_x <- function(p) qlogis(pmin(pmax((p - lo) / (hi - lo), 1e-4), 1 - 1e-4))   # as script 54
  to_p <- function(x) lo + (hi - lo) * plogis(x)

  # The group prior, from the other sample's pooled fits (script 54)
  files <- sprintf("output/pooled/cv/%s/%s_fold%02d.rds", other, model, 1:10)
  files <- files[file.exists(files)]
  if (!length(files) || (length(files) < 10 && !quick)) stop("run script 54 on ", other, " first")
  pri <- lapply(files, readRDS)
  mu <- colMeans(do.call(rbind, lapply(pri, function(r) r$mu)))
  sig2 <- colMeans(do.call(rbind, lapply(pri, function(r) r$sigma^2)))

  par_last <- matrix(NA, length(who), length(lo), dimnames = list(who, m$names))
  for (k in seq_along(who)) {
    i <- who[k]
    ch <- d$choice[i, ]; w <- d$win[i, ]; l <- d$loss[i, ]
    x <- mu
    for (t in 2:100) {
      if ((t - 2) %% every == 0) {                        # refit to trials 1 to t - 1
        obj <- function(z) -sum(m$lp(to_p(z), ch[1:(t - 1)], w[1:(t - 1)], l[1:(t - 1)])) + sum((z - mu)^2 / (2 * sig2))
        f <- tryCatch(optim(x, obj, method = "BFGS", control = list(maxit = 500)), error = function(e) NULL)
        if (!is.null(f) && is.finite(f$value)) x <- f$par
      }
      lp[k, t - 1] <- m$lp(to_p(x), ch[1:t], w[1:t], l[1:t])[t]   # predict choice t
    }
    par_last[k, ] <- to_p(x)
  }
}

out_dir <- file.path("output/forward", if (opt("prior", "") == "") data_name else paste0(data_name, "_prior-", other))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(data = data_name, model = model, who = who, lp = lp, par_last = par_last, prior_from = other,
             every = every),
        sprintf("%s/%s_part%d%s.rds", out_dir, model, part, if (quick) "_quick" else ""))
stopifnot(!anyNA(lp))
round(c(people = length(who), mean_lp_per_choice = mean(lp), chance = log(1 / 4)), 3)
