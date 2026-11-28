# The registered predictions of the addendum to U86SY (predictions/hierarchical.md), scored from script 55c's output
# (output/stan/cv_score_<data>.rds). "Unpooled" = each person fit alone (U86SY); "pooled" = the group-prior fits of
# script 54; "hierarchical" = scripts 55. Differences are summed over participants; SE = sd of the per-participant
# differences x sqrt(n); one-sided tests at alpha = .05, with no correction for multiple tests (the rule of U86SY).
#   S1a  hierarchical NOVA beats unpooled NOVA                                                    met if z > 1.645
#   S1b  hierarchical ORL + alpha beats unpooled ORL + alpha                                      met if z > 1.645
#   S2   hierarchical NOVA beats hierarchical ORL + alpha                                         met if z > 1.645
#   S3a  hierarchical NOVA beats pooled NOVA                                                      met if z > 1.645
#   S3b  hierarchical ORL + alpha beats pooled ORL + alpha                                        met if z > 1.645
#   S4   the shared mix of the eight with hierarchical members is worse than that of the eight unpooled
#                                                                                                 met if z < -1.645
#   S5   hierarchical NOVA alone does worse than the shared mix of the eight unpooled             met if z < -1.645
# The replication sample (N = 440) decides; horstmann and replication438 give the same table, reported beside it.
# Convergence rule (scripts/55_stan_common.R): a prediction that uses a model that counts as not tested is "not
# tested", neither met nor failed. Each prediction is also scored without the people with an R-hat above 1.05 in any
# fit of its models; that version is reported and does not decide.
# Described, not tested: the same comparisons for the other six models, and the group parameters (fits to all trials).
# Usage (from the project folder): Rscript scripts/55d_stan_predictions.R <data>
source("scripts/55_stan_common.R")
data_name <- commandArgs(trailingOnly = TRUE)[1]
sc <- readRDS(sprintf("output/stan/cv_score_%s.rds", data_name))
if (sc$preliminary) cat("PRELIMINARY: scored before every rerun was done\n")
cat("N =", nrow(sc$per), "\n")

predictions <- list(S1a = c("hier.nova", "unpooled.nova", ">"),
                    S1b = c("hier.orl_alpha", "unpooled.orl_alpha", ">"),
                    S2 = c("hier.nova", "hier.orl_alpha", ">"),
                    S3a = c("hier.nova", "pooled.nova", ">"),
                    S3b = c("hier.orl_alpha", "pooled.orl_alpha", ">"),
                    S4 = c("eight_hier", "eight_unpooled", "<"),
                    S5 = c("hier.nova", "eight_unpooled", "<"))
result <- function(z, direction) isTRUE(if (direction == ">") z > 1.645 else z < -1.645)
table <- do.call(rbind, lapply(names(predictions), function(id) {
  p <- predictions[[id]]
  everyone <- stan_compare(sc, p[1], p[2])
  converged <- stan_compare(sc, p[1], p[2], converged_only = TRUE)
  untested <- any(sc$not_tested[c(stan_uses(p[1]), stan_uses(p[2]))])
  data.frame(prediction = id, compared = paste(p[1], "vs", p[2]), met_if = paste("z", p[3], if (p[3] == ">") "1.645" else "-1.645"),
             t(round(everyone, 1)),
             result = if (untested) "not tested" else if (result(everyone[["z"]], p[3])) "met" else "not met",
             n_converged = converged[["n"]], z_converged = round(converged[["z"]], 1),
             result_converged = if (untested) "not tested" else if (result(converged[["z"]], p[3])) "met" else "not met")
}))
print(table, row.names = FALSE)

# Described, not tested: each of the other six models, hierarchical against unpooled and against pooled
others <- setdiff(names(sc$not_tested), c("nova", "orl_alpha"))
described <- rbind(
  t(sapply(others, function(mn) stan_compare(sc, paste0("hier.", mn), paste0("unpooled.", mn)))),
  t(sapply(others, function(mn) stan_compare(sc, paste0("hier.", mn), paste0("pooled.", mn)))))
rownames(described) <- c(paste0(others, "_hier_vs_unpooled"), paste0(others, "_hier_vs_pooled"))
cat("Described, not tested\n")
print(round(described, 1))

# Group parameters from the fits to all trials (the rerun where there is one), where they are done
fit_source <- if (data_name == "replication438") "replication" else data_name
stan_of <- c(orl = "orl", variant = "orl_split", orl_alpha = "orl_alpha", variant_alpha = "orl_split_alpha",
             vse = "vse", vse_split = "vse_split", eef = "eef", nova = "nova")
groups <- list()
for (mn in names(stan_of)) {
  f <- stan_final_file(fit_source, stan_of[[mn]], NA)
  if (!is.na(f)) groups[[mn]] <- readRDS(f)$group
}
if (length(groups)) { cat("Group means (fits to all trials)\n"); print(groups) }

saveRDS(list(data = data_name, preliminary = sc$preliminary, predictions = table, described = described, groups = groups),
        sprintf("output/stan/predictions_%s.rds", data_name))
