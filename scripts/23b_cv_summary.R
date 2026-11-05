# Combine the cross-validation folds and compare ORL vs. the variant.
# Positive differences favor the variant (wl).
quick <- FALSE   # set to FALSE once the full Kodiak runs are in
tag <- if (quick) "_quick" else ""
folds <- if (quick) 1 else 1:10

read_fold <- function(model, fold) {
  readRDS(sprintf("output/cv/cv_%s_fold%02d%s.rds", model, fold, tag))
}

# 1. Fold 0: per-trial PSIS-LOO
loo_orl <- read_fold("orl", 0)$loo
loo_wl <- read_fold("wl", 0)$loo
loo_orl
loo_wl
loo_cmp <- loo::loo_compare(list(orl = loo_orl, wl = loo_wl))
loo_cmp

# 2. Folds 1-10: held-out trials
test_orl <- do.call(rbind, lapply(folds, function(f) read_fold("orl", f)$test))
test_wl <- do.call(rbind, lapply(folds, function(f) read_fold("wl", f)$test))
both <- merge(test_orl, test_wl, by = c("subj", "trial"), suffixes = c(".orl", ".wl"))
nrow(both)   # quick: 100; full: 16200 (every trial held out once)
both$diff <- both$elpd.wl - both$elpd.orl

# Sum within each person, SE across people (a person's trials aren't independent)
per_subj <- tapply(both$diff, both$subj, sum)
cv_diff <- c(elpd_diff = sum(per_subj), se = sqrt(length(per_subj)) * sd(per_subj))
cv_diff
sum(per_subj > 0)   # people better predicted by the variant
binom.test(sum(per_subj > 0), length(per_subj), alternative = "greater")

# 3. Convergence across every run
runs <- expand.grid(model = c("orl", "wl"), fold = c(0, folds), stringsAsFactors = FALSE)
runs$max_rhat <- NA
runs$divergent <- NA
for (r in 1:nrow(runs)) {
  x <- read_fold(runs$model[r], runs$fold[r])
  runs$max_rhat[r] <- x$max_rhat
  runs$divergent[r] <- sum(x$diagnostics$num_divergent)
}
runs

saveRDS(list(loo_cmp = loo_cmp, cv_diff = cv_diff, per_subj = per_subj, runs = runs),
        sprintf("output/cv_summary%s.rds", tag))

# Variant advantage per fold (compare with ORL's max_rhat in runs)
both$fold <- (both$trial - 1) %% 10 + 1
round(tapply(both$diff, both$fold, sum), 1)




