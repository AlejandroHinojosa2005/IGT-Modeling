# Scores the hierarchical cross-validation (scripts/55_stan_fit.R, folds 1-10) on the same held-out choices as every
# other analysis. Each held-out choice is scored by its posterior predictive probability (the mean over posterior
# draws of p(choice), computed by script 55_stan_fit.R), so the scores are comparable with the maximum-likelihood ones.
# For each model and fold it uses the rerun of the convergence rule when there is one (scripts/55_stan_common.R).
#   1. Each model fit hierarchically against the same model unpooled (scripts 35/45, 50b) and pooled (script 54).
#   2. NOVA against ORL + alpha, both hierarchical.
#   3. Shared mixes (script 46's method): the eight with hierarchical members (WSLS unchanged) against the eight unpooled
#      and the eight pooled; hierarchical NOVA alone against the eight unpooled.
#   4. Each fit's diagnostics; which models count as not tested; the people with an R-hat above 1.05 in any fit of
#      each model.
# Every test is reported twice: with everyone (this version decides), and without the people who have an R-hat above
# 1.05 in any fit of the models in that test. The registered predictions are scored by script 55d.
# Needs scripts 46 and 54b to have been run for the same data set (the unpooled and pooled predictions).
# Usage (from the project folder): Rscript scripts/55c_stan_cv_score.R <data> [preliminary]
#   data: horstmann, replication, or replication438 (the replication fits without the second record of each repeated
#         pair; nothing is refit, those two people's scores are dropped)
#   preliminary: score even if some fits that did not pass have not been rerun yet (for a first look only)
library(loo)
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
preliminary <- "preliminary" %in% args
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
source("scripts/55_stan_common.R")
d <- load_kept(data_name)                                # $kept: each participant's row in the data set that was fit
src <- cv_source(data_name)
s <- readRDS(sprintf("output/prereg/stacking_%s.rds", data_name))
pooled_score <- readRDS(sprintf("output/pooled/score_%s.rds", data_name))
rows <- s$rows
n <- max(rows$subj)
stopifnot(n == length(d$kept), nrow(rows) == n * 99, identical(lapply(1:10, stan_test_trials), lapply(1:10, test_trials)))
row_of <- function(i, trials) (i - 1) * 99 + trials - 1

# The Stan models' names for the R models
stan_of <- c(orl = "orl", variant = "orl_split", orl_alpha = "orl_alpha", variant_alpha = "orl_split_alpha",
             vse = "vse", vse_split = "vse_split", eef = "eef", nova = "nova")
hier <- matrix(NA, nrow(rows), length(stan_of), dimnames = list(NULL, names(stan_of)))
diagnostics <- NULL
bad_people <- list()
for (mn in names(stan_of)) {
  bad <- rep(FALSE, n)
  for (f in 1:10) {
    file <- stan_final_file(src, stan_of[[mn]], f)
    if (is.na(file)) stop("missing ", stan_file(src, stan_of[[mn]], f), " (scripts/55_stan_fit.R ", src, " ", stan_of[[mn]], " ", f, ")")
    r <- readRDS(file)
    rerun <- grepl("_long", file)
    if (!rerun && !stan_pass(r$diagnostics) && !preliminary)
      stop(file, " did not pass and has not been rerun yet (scripts/55_stan_fit.R ", src, " ", stan_of[[mn]], " ", f, " long)")
    test <- test_trials(f)
    n_fit <- nrow(r$heldout_lp)
    stopifnot(identical(as.integer(colnames(r$heldout_lp)), as.integer(test)), identical(as.integer(r$fit$subj), seq_len(n_fit)),
              all(d$kept <= n_fit), data_name != src || n_fit == n)
    for (i in 1:n) {
      k <- row_of(i, test)
      stopifnot(all(rows$trial[k] == test), all(rows$subj[k] == i))
      hier[k, mn] <- r$heldout_lp[d$kept[i], ]
    }
    bad <- bad | r$fit$max_rhat[d$kept] > 1.05
    diagnostics <- rbind(diagnostics, data.frame(model = mn, fold = f, rerun = rerun, pass = stan_pass(r$diagnostics),
                                                 unusable = stan_unusable(r$diagnostics), t(r$diagnostics)))
  }
  bad_people[[mn]] <- which(bad)
}
stopifnot(!anyNA(hier))
not_tested <- tapply(diagnostics$unusable, diagnostics$model, any)[names(stan_of)]

nova_unpooled <- if (data_name == "horstmann") readRDS("output/hybrid/score_horstmann.rds")$lpd_h[, "hybrid"] else
  readRDS(sprintf("output/hybrid/nova_%s.rds", data_name))$nova
unpooled <- cbind(s$lpd[, names(models)], nova = nova_unpooled)
pooled <- pooled_score$pooled[, names(stan_of)]
stopifnot(nrow(unpooled) == nrow(rows), nrow(pooled) == nrow(rows))
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
stopifnot(nrow(per) == n, identical(names(models), stan_learning_models))
sc <- list(per = per, bad_people = bad_people)
compare <- function(a, b, converged_only = FALSE) stan_compare(sc, a, b, converged_only)
both <- function(a, b) rbind(everyone = compare(a, b), converged = compare(a, b, TRUE))

# 1. Hierarchical against unpooled and pooled, model by model
cat("Hierarchical against unpooled (everyone)\n")
print(round(t(sapply(names(stan_of), function(mn) compare(paste0("hier.", mn), paste0("unpooled.", mn)))), 1))
cat("Hierarchical against pooled (everyone)\n")
print(round(t(sapply(names(stan_of), function(mn) compare(paste0("hier.", mn), paste0("pooled.", mn)))), 1))
# 2. NOVA against ORL + alpha, both hierarchical
print(round(both("hier.nova", "hier.orl_alpha"), 1))
# 3. Mixes
print(round(rbind(eight_hier_vs_unpooled = compare("eight_hier", "eight_unpooled"),
                  eight_hier_vs_pooled = compare("eight_hier", "eight_pooled"),
                  hier_nova_vs_eight_unpooled = compare("hier.nova", "eight_unpooled")), 1))
# 4. Diagnostics: the fit used for each model and fold (the rerun where there is one)
print(aggregate(cbind(reruns = rerun, not_passed = !pass, unusable = unusable, divergent = divergent, max_rhat = max_rhat,
                      max_rhat_group = max_rhat_group, minutes = minutes) ~ model, diagnostics,
                function(x) round(max(x, na.rm = TRUE), 3), na.action = na.pass))
cat("Models that count as not tested:", if (any(not_tested)) names(stan_of)[not_tested] else "none", "\n")
cat("People with an R-hat above 1.05 in any fit, by model:\n")
print(sapply(bad_people, length))
if (preliminary) cat("PRELIMINARY: some fits that did not pass have not been rerun yet\n")

saveRDS(list(data = data_name, preliminary = preliminary, hier = hier, per = per, diagnostics = diagnostics,
             not_tested = not_tested, bad_people = bad_people, kept = d$kept),
        sprintf("output/stan/cv_score_%s.rds", data_name))
