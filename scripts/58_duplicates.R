# Exploratory check, not part of the registration: the preregistered tests without the two repeated records.
# Two records in the replication sample repeat another participant trial for trial (choices, wins and losses on all
# 100 trials): Kjome's row 17 repeats row 2, and Wood's row 299 repeats row 298 (scripts/43_prereg_data.R,
# duplicate_rows). Found November 24, 2026, when matching our participants to Psych-101. The registered results
# (script 48, N = 440) stand as registered; this reports the same tests on the 438 left when the second record of each
# pair is dropped. Every model is fit to each person separately, so nothing is refit: scripts 46 and 47 are rerun on
# "replication438" (the shared mix and the switching models are estimated across participants).
# Also checks that the two pairs are the only records in either sample that repeat another one.
# Needs: scripts 46 and 47 run on replication and replication438.
# Usage (from the project folder): Rscript scripts/58_duplicates.R
source("scripts/43_prereg_data.R")

# 1. The repeated records: every record compared with every other, in both samples
d <- load_set("replication")
h <- load_set("horstmann")
key <- function(x) apply(cbind(x$choice, x$win, x$loss), 1, paste, collapse = ",")
k_all <- c(paste0("replication:", seq_len(nrow(d$choice))), paste0("horstmann:", seq_len(nrow(h$choice))))
groups <- split(k_all, c(key(d), key(h)))
repeats <- unname(Filter(function(g) length(g) > 1, groups))
repeats
stopifnot(length(repeats) == 2,
          setequal(unlist(repeats), paste0("replication:", c(duplicate_of, duplicate_rows))))
d$study[c(duplicate_of, duplicate_rows)]

# 2. The registered tests, as script 48 computes them, on 440 (as registered) and on 438
per_of <- function(name) {
  a <- readRDS(sprintf("output/prereg/stacking_%s.rds", name))
  b <- readRDS(sprintf("output/prereg/switching_%s.rds", name))
  list(per = cbind(a$per, b$cv), st = a, sw = b)
}
compare <- function(per, x, y) {
  d <- per[, x] - per[, y]
  n <- length(d)
  c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better_for = sum(d > 0), N = n)
}
hyps <- list(H1 = c("shared_mix", "orl_alpha"), H1b = c("shared_mix", "orl_alpha_hedged"), H2 = c("shrink_mix", "shared_mix"),
             H3 = c("switching", "blend"), H4 = c("switching", "shuffled"), H5 = c("personal_switching", "shrink_mix"),
             H7 = c("switching", "drift_curve"))
test_all <- function(per) {
  tab <- as.data.frame(t(sapply(hyps, function(h) compare(per, h[1], h[2]))))
  tab$met <- tab$z > 1.645
  tab
}
show <- function(tab) { tab[, 1:3] <- round(tab[, 1:3], 1); tab }
registered <- test_all(per_of("replication")$per)
deduplicated <- test_all(per_of("replication438")$per)
cat("As registered (N = 440):\n"); show(registered)
cat("Without the two repeated records (N = 438):\n"); show(deduplicated)
cat("Change in z (438 minus 440):\n"); round(deduplicated$z - registered$z, 2)
cat("Same tests met:", identical(registered$met, deduplicated$met), "\n")

saveRDS(list(repeats = repeats, registered = registered, deduplicated = deduplicated),
        "output/prereg/duplicates.rds")
