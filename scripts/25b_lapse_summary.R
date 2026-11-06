# Summarizes the lapse-term check (script 25)
quick <- FALSE
tag <- if (quick) "_quick" else ""
orl <- readRDS("output/orl_fit_horstmann.rds")      # without a lapse term
wl <- readRDS("output/orl_wl_fit_horstmann.rds")
lo <- readRDS(sprintf("output/lapse/lapse_orl%s.rds", tag))   # with a lapse term
lw <- readRDS(sprintf("output/lapse/lapse_wl%s.rds", tag))
n <- nrow(lw)
orl <- orl[1:n, ]
wl <- wl[1:n, ]

# 1. How large is the lapse rate?
summary(lo$eps)
summary(lw$eps)
sum(lw$eps > 0.05)                 # people with more than 5% random choices (variant)

# 2. Does the win/loss asymmetry survive? (without a lapse term: 121 of 162)
sum(lw$betaPw > lw$betaPl)
binom.test(sum(lw$betaPw > lw$betaPl), n, alternative = "greater")
c(without = median(wl$betaPw - wl$betaPl), with = median(lw$betaPw - lw$betaPl))
cor(wl$betaPw - wl$betaPl, lw$betaPw - lw$betaPl, method = "spearman")   # do people keep their split?

# 3. Does the variant still beat ORL once both have a lapse term?
c(orl_eps_aic = sum(2 * lo$nll + 2 * 6), wl_eps_aic = sum(2 * lw$nll + 2 * 7))
c(orl_eps_bic = sum(2 * lo$nll + 6 * log(100)), wl_eps_bic = sum(2 * lw$nll + 7 * log(100)))
sum(lo$nll - lw$nll > 1)               # variant wins per person, AIC
sum(lo$nll - lw$nll > log(100) / 2)    # variant wins per person, BIC

# 4. Is the lapse term worth its extra parameter? (summed AIC without minus with; positive = lapse helps)
c(orl = sum(2 * orl$nll + 10) - sum(2 * lo$nll + 12),
  wl = sum(2 * wl$nll + 12) - sum(2 * lw$nll + 14))
c(orl = sum(2 * orl$nll + 5 * log(100)) - sum(2 * lo$nll + 6 * log(100)),
  wl = sum(2 * wl$nll + 6 * log(100)) - sum(2 * lw$nll + 7 * log(100)))

saveRDS(list(eps_orl = summary(lo$eps), eps_wl = summary(lw$eps),
             split_count = sum(lw$betaPw > lw$betaPl),
             aic = c(orl_eps = sum(2 * lo$nll + 12), wl_eps = sum(2 * lw$nll + 14)),
             bic = c(orl_eps = sum(2 * lo$nll + 6 * log(100)), wl_eps = sum(2 * lw$nll + 7 * log(100)))),
        sprintf("output/lapse_summary%s.rds", tag))
