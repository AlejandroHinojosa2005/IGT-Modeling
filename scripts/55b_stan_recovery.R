# Parameter recovery for the hierarchical Stan models: simulate players with known parameters, fit them with
# scripts/55_stan_fit.R's model, and compare. True parameters: 60 Horstmann participants' maximum-likelihood fits
# (the full-data fits used by script 44; for NOVA, its fold-1 fits from script 50b), seed 2055. Each simulated player
# makes 100 choices; outcomes are cards drawn from the chosen deck's cards in the Horstmann data (deck_draws, as in
# scripts 31, 39 and 44). Reported per parameter: Spearman correlation of true and recovered (posterior mean) values,
# and the share of true values inside the 95% posterior interval.
# Usage (from the project folder): Rscript scripts/55b_stan_recovery.R <model>   (orl, orl_split, ..., nova, eef)
library(cmdstanr)
library(posterior)
model <- commandArgs(trailingOnly = TRUE)[1]
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")                      # also loads 43b: step functions, deck_draws

r_name <- c(orl = "orl", orl_split = "variant", orl_alpha = "orl_alpha", orl_split_alpha = "variant_alpha",
            nova = "nova", vse = "vse", vse_split = "vse_split", eef = "eef")[[model]]
fit_files <- list(orl = "output/orl_fit_horstmann.rds", variant = "output/orl_wl_fit_horstmann.rds",
                  orl_alpha = "output/alpha/alpha_orl.rds", variant_alpha = "output/alpha/alpha_wl.rds",
                  vse = "output/vse/fit_vse.rds", vse_split = "output/vse/fit_vse_split.rds", eef = "output/vse/fit_eef.rds")

# NOVA's step-by-step version (the same equations as hybrid_lp, net bending), checked below
nova_steps <- function(p) list(                          # p: Arew, Apun, Kprime, betaF, betaP, alpha, Q0, D
  start = function() list(Ev = rep(p[7], 4), EF = rep(0, 4), PS = rep(0, 4)),
  probs = function(s) softmax4(s$Ev + s$EF * p[4] + s$PS * p[5]),
  update = function(s, j, win, loss) {
    K <- 3^p[3] - 1; x <- (win + loss) / 100; sg <- sign(x); u <- sg * abs(x)^p[6]
    ac <- if (x >= 0) p[1] else p[2]; au <- if (x >= 0) p[2] else p[1]
    s$Ev[-j] <- p[8] * s$Ev[-j]
    s$Ev[j] <- s$Ev[j] + ac * (u - s$Ev[j])
    s$EF[j] <- s$EF[j] + ac * (sg - s$EF[j])
    s$EF[-j] <- s$EF[-j] + au * (-sg / 3 - s$EF[-j])
    s$PS <- s$PS / (1 + K); s$PS[j] <- 1 / (1 + K)
    s
  })

if (model == "nova") {
  f1 <- do.call(rbind, lapply(Sys.glob("output/hybrid/cv/horstmann/hybrid_fold01_part*.rds"), function(x) readRDS(x)$fit))
  truth_all <- as.matrix(f1[order(f1$subj), hybrid_models$hybrid$names])
  steps <- nova_steps
} else {
  truth_all <- as.matrix(readRDS(fit_files[[r_name]])[, models[[r_name]]$names])
  steps <- steps_of[[r_name]]
}
set.seed(2055)
pick <- sample(nrow(truth_all), 60)
truth <- truth_all[pick, ]
N <- nrow(truth)

simulate <- function(p) {
  st <- steps(p); s <- st$start()
  choice <- win <- loss <- integer(100)
  for (t in 1:100) {
    j <- sample(4, 1, prob = st$probs(s))
    k <- sample(nrow(deck_draws[[j]]), 1)
    choice[t] <- j; win[t] <- deck_draws[[j]][k, 1]; loss[t] <- deck_draws[[j]][k, 2]
    s <- st$update(s, j, win[t], loss[t])
  }
  list(choice = choice, win = win, loss = loss)
}
sims <- lapply(1:N, function(i) simulate(truth[i, ]))
if (model == "nova") {                                   # the step version matches hybrid_lp
  st <- nova_steps(truth[1, ]); s <- st$start(); lp_step <- numeric(100)
  for (t in 1:100) { pr <- st$probs(s); lp_step[t] <- log(pr[sims[[1]]$choice[t]]); s <- st$update(s, sims[[1]]$choice[t], sims[[1]]$win[t], sims[[1]]$loss[t]) }
  stopifnot(max(abs(lp_step - hybrid_lp(truth[1, ], sims[[1]]$choice, sims[[1]]$win, sims[[1]]$loss))) < 1e-10)
}

stan_names <- sub("^Kprime$", "Kp", colnames(truth))
stan_data <- list(N = N, T = 100, Tsubj = rep(100L, N),
                  choice = t(sapply(sims, `[[`, "choice")), win = t(sapply(sims, `[[`, "win")), loss = t(sapply(sims, `[[`, "loss")),
                  train = matrix(1L, N, 100), n_test = 1L, test_t = matrix(1L, N, 1))
mod <- cmdstan_model(sprintf("stan/hier/igt_%s.stan", model))
fit <- mod$sample(data = stan_data, chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                  adapt_delta = 0.95, seed = 2055, refresh = 0, show_messages = FALSE)
summ <- fit$summary(variables = stan_names, "mean", ~quantile(.x, c(.025, .975)), "rhat")
recovery <- t(sapply(seq_along(stan_names), function(k) {
  rows <- match(sprintf("%s[%d]", stan_names[k], 1:N), summ$variable)
  c(spearman = cor(truth[, k], summ$mean[rows], method = "spearman"),
    covered = mean(truth[, k] >= summ$`2.5%`[rows] & truth[, k] <= summ$`97.5%`[rows]))
}))
rownames(recovery) <- stan_names
round(recovery, 2)
c(divergent = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent), max_rhat = round(max(summ$rhat), 3))
dir.create("output/stan/recovery", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(truth = truth, recovered = summ, recovery = recovery), sprintf("output/stan/recovery/%s.rds", model))
