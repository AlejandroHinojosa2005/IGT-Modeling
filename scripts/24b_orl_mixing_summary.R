# Which hierarchical ORL parameters mix poorly (script 24), and do the chains disagree
# about K and betaP together (a trade-off between how fast the trace fades and how strong it is)?
quick <- FALSE
tag <- if (quick) "_quick" else ""
orl <- readRDS("output/orl_fit_horstmann.rds")    # single-person ORL fits, for comparison
wl <- readRDS("output/orl_wl_fit_horstmann.rds")  # single-person variant fits

for (fold in c(0, 8)) {
  m <- readRDS(sprintf("output/orl_mixing_fold%02d%s.rds", fold, tag))
  s <- m$summary
  cat("\n========== Fold", fold, "==========\n")
  cat("Group level (mu_pr / sigma: 1 = Arew, 2 = Apun, 3 = K', 4 = betaF, 5 = betaP):\n")
  print(s[grepl("mu_pr|sigma", s$variable), c("variable", "mean", "rhat", "ess_bulk")], digits = 3, row.names = FALSE)

  cat("\nIndividual parameters with R-hat above 1.05:\n")
  bad <- s[!grepl("mu_pr|sigma", s$variable) & s$rhat > 1.05, ]
  print(bad[order(-bad$rhat), c("variable", "mean", "rhat", "ess_bulk")], digits = 3, row.names = FALSE)

  cat("\nChain means of the K and betaP group parameters (do the chains disagree?):\n")
  g <- m$group[, c("mu_pr[3]", "sigma[3]", "mu_pr[5]", "sigma[5]")]
  print(round(apply(g, 2, function(x) tapply(x, m$chain, mean)), 2))

  cat("\nChain means of K' and betaP for the 5 worst-mixing participants:\n")
  for (i in m$worst_subj) {
    v <- c(paste0("K[", i, "]"), paste0("betaP[", i, "]"))
    cm <- round(apply(m$worst[, v], 2, function(x) tapply(x, m$chain, mean)), 2)
    cat("\nParticipant", i, "| single-person ORL: K' =", round(orl$Kprime[i], 2), " betaP =", round(orl$betaP[i], 2),
        "| variant: betaPw =", round(wl$betaPw[i], 2), " betaPl =", round(wl$betaPl[i], 2), "\n")
    print(cm)
  }
}
