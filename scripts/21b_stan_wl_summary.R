# Summary of the hierarchical variant fit (script 21)
s <- readRDS("output/stan_wl_horstmann.rds")
s$diagnostics                       # want 0 divergences
max(s$summary$rhat)                 # want <= 1.1
min(s$summary$ess_bulk)             # want above ~400

# Group-level asymmetry: is perseveration after gains stronger than after losses, across the population?
d <- s$group_draws$mu_betaPw - s$group_draws$mu_betaPl
mean(d > 0)                         # posterior probability that the difference is positive
quantile(d, c(.025, .5, .975))      # median & 95% credible interval

# Agreement with the maximum-likelihood fits (Spearman), and the share with betaP_win > betaP_loss
mle <- readRDS("output/orl_wl_fit_horstmann.rds")
round(diag(cor(mle[, colnames(s$ind)], s$ind, method = "spearman")), 2)
mean(s$ind[, "betaPw"] > s$ind[, "betaPl"])




