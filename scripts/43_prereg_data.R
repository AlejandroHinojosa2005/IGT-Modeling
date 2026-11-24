# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY, registered November 20, 2026): the data sets.
# load_set(name) returns one data set as matrices (participants x 100 trials) of choices (1-4), wins and losses
# (losses stored as negative numbers), with each participant's study label and payoff scheme.
#   "horstmann":     the discovery sample (N = 162). Used only to check that the new scripts reproduce scripts 35-38.
#   "replication":   every other study in the pool except Fridberg (95-trial sessions). For the two 150-trial studies,
#                    only the first 100 trials (N = 440).
#   "simulated":     the no-switching data made by script 44 (N = 324).
#   "replication342": the replication sample without the two 150-trial studies (robustness check 5).
#   "replication438": the replication sample without the two records that repeat another participant (see
#                     duplicate_rows below; an exploratory check, not part of the registration).
#   "sim_switch", "sim_types": exploratory data WITH switching, made by script 44b (N = 162 each).
# Payoff schemes as described in the data paper (Steingroever et al., 2015): 1 = the original IGT, 2 = the variant with
# a constant loss on deck C and shuffled outcomes, 3 = the Bechara and Damasio scheme in which the decks change over blocks.
scheme_of <- c(Horstmann = 2, Kjome = 3, Maia = 1, Premkumar = 3, SteingroverInPrep = 2,
               Steingroever2011 = 2, Wetzels = 2, Wood = 3, Worthy = 1)

load_set <- function(name) {
  if (name == "simulated") return(readRDS("output/prereg/null_sim_data.rds"))
  if (name %in% c("sim_switch", "sim_types")) return(readRDS(sprintf("output/prereg/%s_data.rds", name)))
  e <- new.env()
  load("data-raw/IGTdata.rdata", envir = e)
  s100 <- e$index_100[, "Study"]
  s150 <- e$index_150[, "Study"]
  if (name == "horstmann") {
    k <- s100 == "Horstmann"
    d <- list(choice = e$choice_100[k, ], win = e$wi_100[k, ], loss = e$lo_100[k, ], study = s100[k])
  } else if (name %in% c("replication", "replication342")) {
    k <- s100 != "Horstmann"
    d <- list(choice = rbind(e$choice_100[k, ], e$choice_150[, 1:100]),
              win = rbind(e$wi_100[k, ], e$wi_150[, 1:100]),
              loss = rbind(e$lo_100[k, ], e$lo_150[, 1:100]),
              study = c(s100[k], s150))
  } else stop("unknown data set: ", name)
  d[c("choice", "win", "loss")] <- lapply(d[c("choice", "win", "loss")], function(x) unname(as.matrix(x)))
  d$study <- unname(d$study)
  d$scheme <- unname(scheme_of[d$study])
  d
}

# Two records in the replication sample repeat another participant trial for trial (choices, wins and losses on all
# 100 trials): Kjome's participant in row 17 repeats row 2, and Wood's in row 299 repeats row 298 (rows of
# load_set("replication")). Found November 24, 2026, when matching our participants to Psych-101, which keeps one
# record of each pair. "replication438" leaves out the second record of each pair.
duplicate_rows <- c(Kjome = 17L, Wood = 299L)
duplicate_of <- c(Kjome = 2L, Wood = 298L)

# Which participants of the replication sample a data set's analysis uses (all of them, the 342 in 100-trial studies,
# or the 438 without the repeated records)
keep_of <- function(name, study) {
  k <- rep(TRUE, length(study))
  if (name == "replication342") k <- !study %in% c("Steingroever2011", "Wetzels")
  if (name == "replication438") {
    stopifnot(length(study) == 440, identical(unname(study[duplicate_rows]), names(duplicate_rows)))
    k[duplicate_rows] <- FALSE
  }
  k
}
# Where a data set's fits are (the subsets reuse the replication fits wherever each person is fit on their own)
cv_source <- function(name) if (name %in% c("replication342", "replication438")) "replication" else name

# A data set's participants as its analyses use them: the source data set without the participants keep_of leaves out
# (for replication438, after checking that the two records left out do repeat their pairs). $kept: their rows in the
# source data set.
load_kept <- function(name) {
  d <- load_set(cv_source(name))
  k <- keep_of(name, d$study)
  if (name == "replication438") for (j in seq_along(duplicate_rows)) for (v in c("choice", "win", "loss"))
    stopifnot(identical(d[[v]][duplicate_rows[j], ], d[[v]][duplicate_of[j], ]))
  d[c("choice", "win", "loss")] <- lapply(d[c("choice", "win", "loss")], function(x) x[k, , drop = FALSE])
  d$study <- d$study[k]
  d$scheme <- d$scheme[k]
  d$kept <- which(k)
  d
}
