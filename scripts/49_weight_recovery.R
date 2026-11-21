# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), exploratory: are personal mixing weights recovered?
# The no-switching players of script 44 were made from known weights (w_mix: 0.7 x each Horstmann participant's own
# stacking weights + 0.3 x the shared weights), two players per participant. Script 46 then estimated every simulated
# player's own weights from their held-out choices, without knowing w_mix. This compares the two:
#   1. For each model, the correlation across players between true and estimated weight.
#   2. For each player, the distance between true and estimated weights (half the summed absolute difference: 0 when
#      they match, 1 when they share nothing), against two baselines: everyone given the shared weights, and the
#      shrinkage estimate (lambda x own + (1 - lambda) x shared, lambda the average chosen by script 46).
#   3. Whether the model with the largest estimated weight is the one with the largest true weight.
#   4. The same comparisons for model families (ORL: ORL, variant, ORL + alpha, variant + alpha; VSE: VSE, VSE-split;
#      EEF; WSLS), because models within a family make similar predictions and can trade weight.
#   5. Reliability: the two players made from the same participant have the same true weights, so the correlation
#      between their estimated weights shows how much of the estimate is the person and how much is noise.
#   6. Prediction: the held-out choices scored with the true weights (the best any estimate could do) against the
#      shared, shrinkage and personal mixes of script 46, so the loss from estimating the weights is in the same
#      units as H2. The same differences in the replication sample are printed beside them.
# No refitting and no random numbers. Usage (from the project folder): Rscript scripts/49_weight_recovery.R
st <- readRDS("output/prereg/stacking_simulated.rds")
sim <- readRDS("output/prereg/null_sim_data.rds")
model_names <- colnames(st$w_person)
true <- sim$w_mix[sim$source, ]                              # 324 players x 8 models
stopifnot(identical(colnames(true), model_names), nrow(true) == nrow(st$w_person))
est <- st$w_person
shared <- matrix(st$w_group, nrow(est), length(model_names), byrow = TRUE, dimnames = list(NULL, model_names))
lam <- mean(st$lambda_chosen)
shrunk <- lam * est + (1 - lam) * shared

dist <- function(a, b) rowSums(abs(a - b)) / 2
families <- list(ORL = c("orl", "variant", "orl_alpha", "variant_alpha"), VSE = c("vse", "vse_split"),
                 EEF = "eef", WSLS = "wsls")
fam <- function(w) sapply(families, function(f) rowSums(w[, f, drop = FALSE]))
twin_r <- function(w) sapply(seq_len(ncol(w)), function(k) cor(w[seq(1, nrow(w), 2), k], w[seq(2, nrow(w), 2), k]))
summarize <- function(true, est, shrunk, shared) {
  list(per_model = round(rbind(r_true_est = diag(cor(true, est)), r_true_shrunk = diag(cor(true, shrunk)),
                               twin_r = setNames(twin_r(est), colnames(est)),
                               true_sd = apply(true, 2, sd)), 2),
       distance = round(c(personal = mean(dist(true, est)), shrinkage = mean(dist(true, shrunk)),
                          shared = mean(dist(true, shared))), 3),
       personal_closer_than_shared = mean(dist(true, est) < dist(true, shared)),
       top_match = round(c(personal = mean(max.col(true, "first") == max.col(est, "first")),
                           shrinkage = mean(max.col(true, "first") == max.col(shrunk, "first")),
                           shared = mean(max.col(true, "first") == max.col(shared, "first"))), 2))
}
lam
by_model <- summarize(true, est, shrunk, shared)
by_family <- summarize(fam(true), fam(est), fam(shrunk), fam(shared))
by_model
by_family

# 6. Prediction with the true weights vs. the estimates (summed held-out log-likelihood; higher is better)
true_score <- log(rowSums(exp(st$lpd) * true[st$rows$subj, ]))
r <- readRDS("output/prereg/stacking_replication.rds")
prediction <- round(rbind(
  simulated = c(true_weights = sum(true_score - st$cv_group), shrinkage = sum(st$cv_shrink - st$cv_group),
                personal = sum(st$cv_person - st$cv_group), lambda = mean(st$lambda_chosen), N = nrow(true)),
  replication = c(true_weights = NA, shrinkage = sum(r$cv_shrink - r$cv_group),
                  personal = sum(r$cv_person - r$cv_group), lambda = mean(r$lambda_chosen), N = max(r$rows$subj))), 2)
prediction                                               # each column minus the shared mix
saveRDS(list(lambda = lam, by_model = by_model, by_family = by_family, prediction = prediction),
        "output/prereg/weight_recovery.rds")
