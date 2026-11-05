# Compares the hierarchical (hBayesDM, script 16) and maximum-likelihood estimates;
# convergence and LOOIC reliability.
hb_pvl <- readRDS("output/hbayesdm_pvl_horstmann.rds")
hb_vpp <- readRDS("output/hbayesdm_vpp_horstmann.rds")
hb_orl <- readRDS("output/hbayesdm_orl_horstmann.rds")
mle_pvl <- readRDS("output/pvl_delta_fit_horstmann.rds")
mle_vpp <- readRDS("output/vpp_fit_horstmann.rds")
mle_orl <- readRDS("output/orl_fit_horstmann.rds")

# Rank agreement (Spearman) between the two methods, per parameter
# (hBayesDM's ORL "K" is on the same 0-5 scale as my Kprime)
round(diag(cor(mle_pvl[, c("A", "alpha", "cons", "lambda")], hb_pvl$ind[, c("A", "alpha", "cons", "lambda")], method = "spearman")), 2)
round(diag(cor(mle_vpp[, c("A", "alpha", "cons", "lambda", "epP", "epN", "K", "w")],
               hb_vpp$ind[, c("A", "alpha", "cons", "lambda", "epP", "epN", "K", "w")], method = "spearman")), 2)
round(diag(cor(mle_orl[, c("Arew", "Apun", "Kprime", "betaF", "betaP")],
               hb_orl$ind[, c("Arew", "Apun", "K", "betaF", "betaP")], method = "spearman")), 2)   # A_rew, A_pun, K', betaF, betaP

# Convergence: parameters with Rhat above 1.1
rownames(hb_pvl$rhat)[hb_pvl$rhat$Rhat > 1.1]
rownames(hb_vpp$rhat)[hb_vpp$rhat$Rhat > 1.1]
rownames(hb_orl$rhat)[hb_orl$rhat$Rhat > 1.1]

# LOOIC reliability: participants with Pareto k above 0.7 (of 162)
sum(hb_pvl$looic$LOOIC$diagnostics$pareto_k > 0.7)
sum(hb_vpp$looic$LOOIC$diagnostics$pareto_k > 0.7)
sum(hb_orl$looic$LOOIC$diagnostics$pareto_k > 0.7)
readRDS("output/hbayesdm_looic_horstmann.rds")







