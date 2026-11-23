# Checks the hierarchical Stan models before any fit:
#   1. every model file compiles;
#   2. the Stan likelihood functions give the same log-probability for every trial as the R functions used for all
#      maximum-likelihood fits (scripts 43b and 50), for random parameters within the limits and every Horstmann
#      participant (stan/hier/check_lp.stan, run with the fixed_param sampler). Differences should be about 1e-12.
# Usage (from the project folder): Rscript scripts/55a_stan_check.R
library(cmdstanr)
library(posterior)
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
d <- load_set("horstmann")
n <- nrow(d$choice)

# 1. Compile every model
stan_models <- c("orl", "orl_split", "orl_alpha", "orl_split_alpha", "nova", "vse", "vse_split", "eef")
for (m in stan_models) cmdstan_model(sprintf("stan/hier/igt_%s.stan", m))
cat("all", length(stan_models), "models compile\n")

# 2. Likelihoods: Stan against R
check <- cmdstan_model("stan/hier/check_lp.stan")
set.seed(55)
u <- function(lo, hi) runif(1, lo, hi)
worst <- c()
for (m in stan_models) {
  par9 <- matrix(0, n, 9); lp_r <- matrix(NA, n, 100)
  for (i in 1:n) {
    A <- c(u(.01, 1), u(.01, 1), u(0, 5), u(-10, 10), u(-10, 10), u(-10, 10), u(.01, 2), u(-5, 5), u(0, 1))
    V <- c(u(0, 1), u(0, 1), u(0, 1), u(-10, 10), u(0, 5), u(-10, 10), 0, 0, 0)
    E <- c(u(0, 1), u(0, 1), u(-5, 5), u(0, 5), 0, 0, 0, 0, 0)
    ch <- d$choice[i, ]; w <- d$win[i, ]; l <- d$loss[i, ]
    if (m == "orl")             { p <- c(A[1:5], A[5], 1, 0, 1);    r <- models$orl$lp(A[1:5], ch, w, l) }
    if (m == "orl_split")       { p <- c(A[1:6], 1, 0, 1);          r <- models$variant$lp(A[1:6], ch, w, l) }
    if (m == "orl_alpha")       { p <- c(A[1:5], A[5], A[7], 0, 1); r <- models$orl_alpha$lp(c(A[1:5], A[7]), ch, w, l) }
    if (m == "orl_split_alpha") { p <- c(A[1:7], 0, 1);             r <- models$variant_alpha$lp(A[1:7], ch, w, l) }
    if (m == "nova")            { p <- c(A[1:5], A[5], A[7:9]);     r <- hybrid_models$hybrid$lp(c(A[1:5], A[7:9]), ch, w, l) }
    if (m == "vse")             { p <- c(V[1:5], 0, 0, 0, 0);       r <- models$vse$lp(V[1:5], ch, w, l) }
    if (m == "vse_split")       { p <- V;                           r <- models$vse_split$lp(V[1:6], ch, w, l) }
    if (m == "eef")             { p <- E;                           r <- models$eef$lp(E[1:4], ch, w, l) }
    par9[i, ] <- p; lp_r[i, ] <- r
  }
  family <- if (m %in% c("vse", "vse_split")) 2L else if (m == "eef") 3L else 1L
  fit <- check$sample(data = list(N = n, T = 100, choice = matrix(as.integer(d$choice), n, 100), win = d$win,
                                  loss = d$loss, family = family, par = par9),
                      fixed_param = TRUE, chains = 1, iter_sampling = 1, refresh = 0, show_messages = FALSE)
  dr <- fit$draws("lp", format = "draws_matrix")
  lp_s <- matrix(NA, n, 100)
  for (i in 1:n) lp_s[i, ] <- as.numeric(dr[1, sprintf("lp[%d,%d]", i, 1:100)])
  worst[m] <- max(abs(lp_s - lp_r))
}
signif(worst, 3)                                        # largest difference per model (should be about 1e-12)
stopifnot(all(worst < 1e-8))
cat("Stan and R likelihoods agree for every model\n")
