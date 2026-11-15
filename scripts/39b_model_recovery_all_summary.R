# Summarizes script 39: when each model generated the data, which model wins on AIC and on BIC?
# Rows = the model that generated the data, columns = the model that won (share of the 162 synthetic participants).
models <- c("orl", "variant", "orl_alpha", "variant_alpha", "vse", "vse_split", "eef", "wsls")
k <- c(orl = 5, variant = 6, orl_alpha = 6, variant_alpha = 7, vse = 5, vse_split = 6, eef = 4, wsls = 2)

runs <- lapply(models, function(g) readRDS(sprintf("output/model_recovery_all/%s.rds", g))$res)
names(runs) <- models
nll <- lapply(runs, function(r) as.matrix(r[, paste0("nll.", models)]))
aic <- lapply(nll, function(x) sweep(2 * x, 2, 2 * k, "+"))
bic <- lapply(nll, function(x) sweep(2 * x, 2, k * log(100), "+"))

winners <- function(crit) {
  t(sapply(crit, function(x) table(factor(models[apply(x, 1, which.min)], levels = models)) / nrow(x)))
}
round(winners(aic), 2)   # AIC: share of synthetic participants each model won
round(winners(bic), 2)   # BIC

# Summed over the 162 synthetic participants: which model has the lowest total?
summed <- function(crit) sapply(crit, function(x) models[which.min(colSums(x))])
rbind(AIC = summed(aic), BIC = summed(bic))

# Convergence codes other than 0, per generating model (52 = line search stopped early)
sapply(runs, function(r) sum(as.matrix(r[, paste0("conv.", models)]) != 0, na.rm = TRUE))

saveRDS(list(aic = winners(aic), bic = winners(bic), summed = rbind(AIC = summed(aic), BIC = summed(bic))),
        "output/model_recovery_all_summary.rds")
