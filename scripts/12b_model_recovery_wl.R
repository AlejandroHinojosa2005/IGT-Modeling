# 12b_model_recovery_wl.R (09b copied and edited as needed)
# Summarizes model recovery: which model wins when each model generated the data?

results_mr <- readRDS("output/model_recovery_wl.rds")
gens <- c("orl", "wl")
results_mr$gen <- factor(results_mr$gen, levels = gens)

aic_cols <- c("aic.orl", "aic.wl")
bic_cols <- c("bic.orl", "bic.wl")

aic_win <- factor(gens[apply(results_mr[, aic_cols], 1, which.min)], levels = gens)
bic_win <- factor(gens[apply(results_mr[, bic_cols], 1, which.min)], levels = gens)

# Per-participant winners: rows = model that generated the data, columns = model that won
table(generated = results_mr$gen, won_AIC = aic_win)
table(generated = results_mr$gen, won_BIC = bic_win)

# Same tables as proportions of each row
round(prop.table(table(results_mr$gen, aic_win), 1), 2)
round(prop.table(table(results_mr$gen, bic_win), 1), 2)

# Summed scores for each generating model
aggregate(results_mr[, c(aic_cols, bic_cols)], by = list(generated = results_mr$gen), FUN = sum)