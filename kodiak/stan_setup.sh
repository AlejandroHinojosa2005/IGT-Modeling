#!/bin/bash
# One-time setup for the hierarchical Stan fits on Kodiak (kodiak/stan_cv.sh). Run it on the login node, from the
# project folder, after git pull:   bash kodiak/stan_setup.sh
# 1. installs the R packages cmdstanr and posterior, and the CmdStan compiler (into your home folder; about 10 minutes);
# 2. compiles the eight models once, so the array jobs find them ready instead of all compiling at the same time;
# 3. checks that each model's likelihood equals the R functions' (scripts/55a_stan_check.R; minutes).
# If CmdStan fails to build, a newer compiler is usually the fix: module avail gcc, then module load gcc/<version>.
set -e
module load R/4.5.2
Rscript -e 'install.packages(c("cmdstanr", "posterior"), repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))'
Rscript -e 'cmdstanr::install_cmdstan(cores = 4)'
Rscript -e 'library(cmdstanr); for (f in Sys.glob("stan/hier/*.stan")) invisible(cmdstan_model(f)); cat("all models compiled\n")'
Rscript scripts/55a_stan_check.R
