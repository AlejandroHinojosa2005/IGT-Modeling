# Exports the trial-level IGT data of the discovery (Horstmann, N = 162) and replication (N = 440) samples to
# data-raw/igt_trials.csv, for the Centaur benchmark (centaur/). Same participants, order and trials as script 43
# (the two 150-trial studies cut to their first 100 trials). The file stays in data-raw/, which is never pushed.
# Usage (from the project folder): Rscript scripts/53_export_trials.R
source("scripts/43_prereg_data.R")
out <- do.call(rbind, lapply(c("horstmann", "replication"), function(s) {
  d <- load_set(s)
  n <- nrow(d$choice)
  data.frame(sample = s, subj = rep(1:n, each = 100), study = rep(d$study, each = 100), trial = rep(1:100, n),
             choice = as.vector(t(d$choice)), win = as.vector(t(d$win)), loss = as.vector(t(d$loss)))
}))
stopifnot(nrow(out) == (162 + 440) * 100, !anyNA(out))
write.csv(out, "data-raw/igt_trials.csv", row.names = FALSE)
table(out$sample[out$trial == 1])
