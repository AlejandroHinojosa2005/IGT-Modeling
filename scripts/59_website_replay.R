# For the website's "Watch a real player" (alejandrohinojosa.com/research.html): the 40 participants of Maia's study,
# which is not in Psych-101 at all, so Centaur never saw any of it. For each person: their 100 choices and outcomes,
# and for every choice from trial 2 on, the probability each model gave the deck they chose, as scored in the paper:
# the eight models' forward predictions (script 57: refit after every choice, Horstmann prior), the shared mix of the
# eight (script 57b's weights, fit on Horstmann), and Centaur's (centaur/score_centaur.py). Probabilities are stored in
# thousandths. The choices are public data (Maia & McClelland, 2004, in the many-labs pool of Steingroever et al., 2015,
# CC BY-SA 4.0); nothing here comes from Psych-101's gated test set.
# Needs: script 57 on replication, scripts 46 (horstmann) and 50c, and centaur/centaur70b_scores.csv.
# Usage (from the project folder): Rscript scripts/59_website_replay.R   (writes output/website/realplayers.json)
library(loo)
source("scripts/43_prereg_data.R")
source("scripts/50_hybrid_model.R")
d <- load_set("replication")
who <- which(d$study == "Maia")

fw <- function(m) {
  x <- do.call(rbind, lapply(Sys.glob(sprintf("output/forward/replication/%s_part*.rds", m)), function(f) readRDS(f)$lp))
  x[as.character(who), , drop = FALSE]
}
FW <- setNames(lapply(model_names, fw), model_names)

# The shared mix of the eight, with the weights fit on Horstmann (as script 57b)
s_h <- readRDS("output/prereg/stacking_horstmann.rds")
w <- as.numeric(stacking_weights(s_h$lpd[, model_names]))
a <- simplify2array(lapply(seq_along(model_names), function(k) FW[[model_names[k]]] + log(w[k])))
mx <- apply(a, 1:2, max)
FW$mix <- mx + log(apply(exp(a - array(mx, dim(a))), 1:2, sum))

cz <- read.csv("centaur/centaur70b_scores.csv")
cz <- cz[cz$sample == "replication" & cz$subj %in% who & cz$trial >= 2, ]
C <- matrix(NA, length(who), 99, dimnames = list(as.character(who), 2:100))
C[cbind(match(cz$subj, who), cz$trial - 1)] <- cz$lp_4
stopifnot(!anyNA(C))
FW$centaur <- C

# Check against script 57b's numbers for these 40 (probability per choice, geometric mean)
round(sapply(FW, function(x) exp(mean(x))), 3)

players <- lapply(seq_along(who), function(k) {
  i <- who[k]
  list(c = paste(d$choice[i, ], collapse = ""), w = d$win[i, ], l = d$loss[i, ],
       p = lapply(FW, function(x) as.integer(pmax(1, round(1000 * exp(x[k, ]))))))
})
# JSON written by hand (no extra package)
arr <- function(x) paste0("[", paste(x, collapse = ","), "]")
one <- function(pl) sprintf('{"c":"%s","w":%s,"l":%s,"p":{%s}}', pl$c, arr(pl$w), arr(pl$l),
                            paste(sprintf('"%s":%s', names(pl$p), sapply(pl$p, arr)), collapse = ","))
json <- sprintf('{"source":"Maia & McClelland (2004), many-labs IGT pool (Steingroever et al., 2015), CC BY-SA 4.0","note":"p: probability (in thousandths) each model gave the deck chosen, trials 2-100","players":[%s]}',
                paste(sapply(players, one), collapse = ",\n"))
dir.create("output/website", recursive = TRUE, showWarnings = FALSE)
writeLines(json, "output/website/realplayers.json")
