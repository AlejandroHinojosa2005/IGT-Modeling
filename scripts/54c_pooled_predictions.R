# Exploratory, with predictions committed before running (predictions/pooled_replication.md): the pooled fits of
# script 54 on the replication sample (N = 440), as scored by script 54b (output/pooled/score_<data>.rds).
# The four tests (differences summed over participants; SE = sd of the per-participant differences x sqrt(n)):
#   Q1  pooling improves NOVA (pooled against unpooled NOVA)                                      met if z > 1.645
#   Q2  pooled NOVA beats pooled ORL + alpha                                                      met if z > 1.645
#   Q3  the shared mix of the eight pooled models is worse than that of the eight unpooled        met if z < -1.645
#   Q4  five pooled models with NOVA do not beat the shared mix of the eight unpooled             met if z < 1.645
#   Q5  pooled NOVA alone does worse than the shared mix of the eight unpooled                    met if z < -1.645
#   (one-sided tests at alpha = .05, the rule of the preregistered tests in script 48; no correction for multiple tests)
# Also described, not tested: pooled NOVA against the pooled five, each model's gain from pooling, and the group
# parameters.
# Usage (from the project folder): Rscript scripts/54c_pooled_predictions.R <data>
#   replication for the tests; horstmann gives the same numbers for the discovery sample
data_name <- commandArgs(trailingOnly = TRUE)[1]
r <- readRDS(sprintf("output/pooled/score_%s.rds", data_name))
per <- r$per
n <- nrow(per)
compare <- function(a, b) {
  x <- per[, a] - per[, b]
  c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n)), better_for = sum(x > 0))
}

tests <- rbind(Q1_pooling_improves_nova = compare("pooled.nova", "unpooled.nova"),
               Q2_pooled_nova_vs_pooled_orl_alpha = compare("pooled.nova", "pooled.orl_alpha"),
               Q3_eight_pooled_vs_eight_unpooled = compare("eight_pooled", "eight_unpooled"),
               Q4_five_nova_pooled_vs_eight_unpooled = compare("five_nova_pooled", "eight_unpooled"),
               Q5_pooled_nova_vs_eight_unpooled = compare("pooled.nova", "eight_unpooled"))
met <- c(tests[1:2, "z"] > 1.645, tests[3, "z"] < -1.645, tests[4, "z"] < 1.645, tests[5, "z"] < -1.645)
predictions <- data.frame(round(tests, 1), met = met)
cat("N =", n, "\n")
predictions

# Described, not tested
described <- round(rbind(pooled_nova_vs_five_nova_pooled = compare("pooled.nova", "five_nova_pooled")), 1)
described
r$pooling_gain
lapply(r$groups, round, 3)

saveRDS(list(data = data_name, predictions = predictions, described = described),
        sprintf("output/pooled/predictions_%s.rds", data_name))
