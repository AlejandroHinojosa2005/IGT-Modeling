# Summarizes the repeated recovery (script 26) across replicates.
# Single-replicate results for comparison: false positives 20% (AIC) and 3% (BIC), scripts 12/12b;
# variant parameter recovery r = .90 (betaPw) and .83 (betaPl), script 13.
quick <- FALSE
rep_ids <- if (quick) 1:2 else 1:10
reps <- lapply(rep_ids, function(r)
  readRDS(sprintf("output/recovery_reps/rep_%02d%s.rds", r, if (quick) "_quick" else "")))

spread <- function(m) round(rbind(mean = colMeans(m), min = apply(m, 2, min), max = apply(m, 2, max)), 3)

# 1. Model recovery, per replicate
mr <- t(sapply(reps, function(r) {
  s <- r$scores
  wl_aic <- s$aic.wl < s$aic.orl    # variant wins this simulated person (AIC)
  wl_bic <- s$bic.wl < s$bic.orl
  c(false_pos_aic = mean(wl_aic[s$gen == "orl"]),   # ORL generated the data, variant won
    false_pos_bic = mean(wl_bic[s$gen == "orl"]),
    hit_aic = mean(wl_aic[s$gen == "wl"]),          # variant generated the data, variant won
    hit_bic = mean(wl_bic[s$gen == "wl"]),
    # summed scores: positive = the summed criterion picked the generating model
    orl_gen_aic_margin = sum(s$aic.wl[s$gen == "orl"]) - sum(s$aic.orl[s$gen == "orl"]),
    wl_gen_aic_margin = sum(s$aic.orl[s$gen == "wl"]) - sum(s$aic.wl[s$gen == "wl"]),
    orl_gen_bic_margin = sum(s$bic.wl[s$gen == "orl"]) - sum(s$bic.orl[s$gen == "orl"]),
    wl_gen_bic_margin = sum(s$bic.orl[s$gen == "wl"]) - sum(s$bic.wl[s$gen == "wl"]))
}))
round(mr, 3)
spread(mr)

# 2. Variant parameter recovery, per replicate (Pearson r, then Spearman)
par_wl <- t(sapply(reps, function(r) diag(cor(r$true$wl, r$recovered$wl))))
par_wl_s <- t(sapply(reps, function(r) diag(cor(r$true$wl, r$recovered$wl, method = "spearman"))))
spread(par_wl)
spread(par_wl_s)

# 3. The win/loss split (betaPw - betaPl), per replicate
split <- t(sapply(reps, function(r) {
  td <- r$true$wl[, "betaPw"] - r$true$wl[, "betaPl"]
  rd <- r$recovered$wl[, "betaPw"] - r$recovered$wl[, "betaPl"]
  c(r = cor(td, rd), rho = cor(td, rd, method = "spearman"),
    sign_agree = mean(sign(td) == sign(rd)), tilt = mean(rd - td))
}))
spread(split)

# 4. ORL parameter recovery, per replicate (Pearson r)
par_orl <- t(sapply(reps, function(r) diag(cor(r$true$orl, r$recovered$orl))))
spread(par_orl)

saveRDS(list(model_recovery = mr, par_wl = par_wl, par_wl_spearman = par_wl_s, split = split, par_orl = par_orl),
        sprintf("output/recovery_reps_summary%s.rds", if (quick) "_quick" else ""))
