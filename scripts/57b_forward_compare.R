# Compares forward predictions (each choice predicted only from the trials before it): our models (script 57), the
# shared mixes of them, and Centaur (centaur/score_centaur.py), on the same participants and choices.
# Mixes: the eight models (with WSLS), and the five with NOVA (NOVA, ORL + split + alpha, VSE-split, EEF, WSLS), with
# the shared-mix weights fit in the OTHER sample (script 46's stacking of its held-out predictions; NOVA's from scripts
# 50c and 50d), so no weight is fit to the people being scored. Each choice's mixed probability is the weighted mean
# of the members' probabilities.
# Each comparison is a difference in log-probability summed over the scored choices of each person, then summed over
# people; SE = sd of the per-person differences x sqrt(n), as in scripts 48 and 50d.
# Usage (from the project folder):
#   Rscript scripts/57b_forward_compare.R <data> [centaur_scores.csv] [prior=NAME] [quick]
#     data: horstmann, replication, or replication438 (the replication sample's forward predictions without the two
#           repeated records; they use the Horstmann prior, so nothing is refit)
#     prior=NAME: the forward predictions made with script 57's prior=NAME, with the mix weights fit on NAME
#                 (e.g. horstmann prior=replication438)
#     centaur_scores.csv: the output of score_centaur.py (columns sample, subj, split, trial, deck, lp_token, lp_4);
#                         lp_4, Centaur's probability renormalized over the four decks, is compared with our models
#     quick: use script 57's quick output, for testing
library(loo)
args <- commandArgs(trailingOnly = TRUE)
data_name <- args[1]
quick <- "quick" %in% args
prior <- sub("^prior=", "", grep("^prior=", args, value = TRUE)[1])
centaur_file <- setdiff(args[-1], c("quick", grep("^prior=", args, value = TRUE)))
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
other <- if (!is.na(prior)) prior else if (data_name == "horstmann") "replication" else "horstmann"
fw_dir <- if (!is.na(prior)) paste0(data_name, "_prior-", prior) else cv_source(data_name)
tag <- if (!is.na(prior)) paste0(data_name, "_prior-", prior) else data_name
if (length(centaur_file)) tag <- paste0(tag, "_centaur")       # keep the comparison without Centaur (all participants) as it is

# Our models' forward predictions (all parts), on the participants every model has been scored on
fw_names <- c(names(models), "nova", "wsls")
load_fw <- function(mn) {
  files <- Sys.glob(sprintf("output/forward/%s/%s_part*.rds", fw_dir, mn))
  files <- files[grepl("_quick", files) == quick]
  if (!length(files)) return(NULL)
  do.call(rbind, lapply(files, function(f) readRDS(f)$lp))
}
FW <- Filter(Negate(is.null), setNames(lapply(fw_names, load_fw), fw_names))
stopifnot(length(FW) > 0)
subj <- Reduce(intersect, lapply(FW, rownames))
subj <- subj[order(as.integer(subj))]
if (data_name == "replication438") subj <- setdiff(subj, as.character(duplicate_rows))
FW <- lapply(FW, function(x) x[subj, , drop = FALSE])
cat("models:", names(FW), "\npeople:", length(subj), "\n")

# Shared mixes, with the other sample's weights
s_other <- readRDS(sprintf("output/prereg/stacking_%s.rds", other))
nova_other <- if (other == "horstmann") readRDS("output/hybrid/score_horstmann.rds")$lpd_h[, "hybrid"] else
  readRDS(sprintf("output/hybrid/nova_%s.rds", other))$nova
L_other <- cbind(s_other$lpd, nova = nova_other)
mix_lp <- function(members) {
  w <- as.numeric(stacking_weights(L_other[, members]))
  a <- simplify2array(lapply(seq_along(members), function(k) FW[[members[k]]] + log(w[k])))   # people x trials x models
  mx <- apply(a, 1:2, max)
  mx + log(apply(exp(a - array(mx, dim(a))), 1:2, sum))
}
mix_sets <- list(mix_eight = model_names, mix_five_nova = c("nova", "variant_alpha", "vse_split", "eef", "wsls"))
for (mn in names(mix_sets)) if (all(mix_sets[[mn]] %in% names(FW))) FW[[mn]] <- mix_lp(mix_sets[[mn]])

# Centaur, on the same people and choices
scored <- matrix(TRUE, length(subj), 99)
if (length(centaur_file)) {
  cz <- read.csv(centaur_file[1])
  cz <- cz[cz$sample == cv_source(data_name) & cz$trial >= 2 & cz$trial <= 100 & as.character(cz$subj) %in% subj, ]
  C <- matrix(NA, length(subj), 99, dimnames = list(subj, 2:100))
  C[cbind(match(as.character(cz$subj), subj), cz$trial - 1)] <- cz$lp_4
  scored <- !is.na(C)
  keep <- rowSums(scored) > 0
  FW <- lapply(FW, function(x) x[keep, , drop = FALSE])
  scored <- scored[keep, , drop = FALSE]
  FW$centaur <- C[keep, , drop = FALSE]
  cat("people Centaur scored:", sum(keep), "\n")
}

# Per-person totals over the scored choices, and the comparisons
per <- sapply(FW, function(x) rowSums(ifelse(scored, x, 0)))
per <- matrix(per, ncol = length(FW), dimnames = list(rownames(FW[[1]]), names(FW)))
n <- nrow(per)
choices <- sum(scored)
summary_tab <- data.frame(total = colSums(per), per_choice = colSums(per) / choices,
                          probability = exp(colSums(per) / choices))
summary_tab <- rbind(summary_tab, chance = c(log(1 / 4) * choices, log(1 / 4), 1 / 4))
round(summary_tab[order(-summary_tab$total), ], 3)
compare <- function(a, b) {
  x <- per[, a] - per[, b]
  c(diff = sum(x), se = sd(x) * sqrt(n), z = sum(x) / (sd(x) * sqrt(n)), better_for = sum(x > 0))
}
ref <- if ("centaur" %in% colnames(per)) "centaur" else if ("mix_eight" %in% colnames(per)) "mix_eight" else colnames(per)[1]
against_ref <- round(t(sapply(setdiff(colnames(per), ref), function(a) compare(a, ref))), 1)
cat("each against", ref, "\n")
against_ref

saveRDS(list(data = data_name, per = per, summary = summary_tab, against = against_ref, reference = ref),
        sprintf("output/forward/compare_%s%s.rds", tag, if (quick) "_quick" else ""))
