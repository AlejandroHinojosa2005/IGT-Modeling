# Scores the hierarchical cross-validation (scripts/55_stan_fit.R, folds 1-10) on the same held-out choices as every
# other analysis. Each held-out choice is scored by its posterior predictive probability (the mean over posterior
# draws of p(choice), computed by script 55_stan_fit.R), so the scores are comparable with the maximum-likelihood ones.
#   1. Each model fit hierarchically against the same model unpooled (scripts 35/45, 50b) and pooled (script 54).
#   2. NOVA against ORL + alpha, both hierarchical.
#   3. Shared mixes (script 46's method): the eight with hierarchical members (WSLS unchanged) against the eight unpooled
#      and the eight pooled.
#   4. Each fit's diagnostics: divergent transitions, the largest R-hat and the smallest bulk effective sample size.
# Needs scripts 46 and 54b to have been run for the same data set (the unpooled and pooled predictions).
# Usage (from the project folder): Rscript scripts/55c_stan_cv_score.R <data>   (horstmann or replication)
library(loo)
data_name <- commandArgs(trailingOnly = TRUE)[1]
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
pooled_score <- readRDS(sprintf("output/pooled/score_%s.rds", data_name))
rows <- s$rows
n <- max(rows$subj)
stopifnot(nrow(rows) == n * 99)
row_of <- function(i, trials) (i - 1) * 99 + trials - 1

# The Stan models' names for the R models
stan_of <- c(orl = "orl", variant = "orl_split", orl_alpha = "orl_alpha", variant_alpha = "orl_split_alpha",
             vse = "vse", vse_split = "vse_split", eef = "eef", nova = "nova")
hier <- matrix(NA, nrow(rows), length(stan_of), dimnames = list(NULL, names(stan_of)))
diagnostics <- NULL
for (mn in names(stan_of)) for (f in 1:10) {
  file <- sprintf("output/stan/cv/%s/%s_fold%02d.rds", data_name, stan_of[[mn]], f)
  if (!file.exists(file)) stop("missing ", file, " (scripts/55_stan_fit.R ", data_name, " ", stan_of[[mn]], " ", f, ")")
  r <- readRDS(file)
  test <- test_trials(f)
  stopifnot(nrow(r$heldout_lp) == n, identical(as.integer(colnames(r$heldout_lp)), as.integer(test)))
  for (i in 1:n) {
    k <- row_of(i, test)
    stopifnot(all(rows$trial[k] == test), all(rows$subj[k] == i))
    hier[k, mn] <- r$heldout_lp[i, ]
  }
  diagnostics <- rbind(diagnostics, data.frame(model = mn, fold = f, t(r$diagnostics)))
}
stopifnot(!anyNA(hier))

unpooled <- cbind(s$lpd[, names(models)], nova = if (data_name == "horstmann")
  readRDS("output/hybrid/score_horstmann.rds")$lpd_h[, "hybrid"] else readRDS("output/hybrid/nova_replication.rds")$nova)
pooled <- pooled_score$pooled[, names(stan_of)]
all_lp <- cbind(hier, unpooled[, names(stan_of)], pooled)
colnames(all_lp) <- c(paste0("hier.", names(stan_of)), paste0("unpooled.", names(stan_of)), paste0("pooled.", names(stan_of)))

# Shared mixes, as script 46
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
wsls <- s$lpd[, "wsls"]
mixes <- cbind(eight_hier = shared_cv(cbind(hier[, names(models)], wsls = wsls)),
               eight_pooled = shared_cv(cbind(pooled[, names(models)], wsls = wsls)),
               eight_unpooled = s$cv_group)
per <- rowsum(cbind(all_lp, mixes), rows$subj)
compare <- function(a, b) {
  x <- per[, a] - per[, b]
  c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n)), better_for = sum(x > 0))
}

# 1. Hierarchical against unpooled and pooled, model by model
cat("Hierarchical against unpooled\n")
round(t(sapply(names(stan_of), function(mn) compare(paste0("hier.", mn), paste0("unpooled.", mn)))), 1)
cat("Hierarchical against pooled\n")
round(t(sapply(names(stan_of), function(mn) compare(paste0("hier.", mn), paste0("pooled.", mn)))), 1)
# 2. NOVA against ORL + alpha, both hierarchical
round(compare("hier.nova", "hier.orl_alpha"), 1)
# 3. Mixes
round(rbind(eight_hier_vs_unpooled = compare("eight_hier", "eight_unpooled"),
            eight_hier_vs_pooled = compare("eight_hier", "eight_pooled"),
            hier_nova_vs_eight_unpooled = compare("hier.nova", "eight_unpooled")), 1)
# 4. Diagnostics: fits with any divergent transition or an R-hat above 1.05
diagnostics$flag <- diagnostics$divergent > 0 | diagnostics$max_rhat > 1.05
aggregate(cbind(flagged = flag, divergent = divergent, max_rhat = max_rhat) ~ model, diagnostics,
          function(x) round(max(x), 3))

saveRDS(list(data = data_name, hier = hier, per = per, diagnostics = diagnostics),
        sprintf("output/stan/cv_score_%s.rds", data_name))
