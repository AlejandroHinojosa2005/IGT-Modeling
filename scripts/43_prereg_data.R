# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY, registered November 20, 2026): the data sets.
# load_set(name) returns one data set as matrices (participants x 100 trials) of choices (1-4), wins and losses
# (losses stored as negative numbers), with each participant's study label and payoff scheme.
#   "horstmann":     the discovery sample (N = 162). Used only to check that the new scripts reproduce scripts 35-38.
#   "replication":   every other study in the pool except Fridberg (95-trial sessions). For the two 150-trial studies,
#                    only the first 100 trials (N = 440).
#   "simulated":     the no-switching data made by script 44 (N = 324).
#   "replication342": the replication sample without the two 150-trial studies (robustness check 5).
# Payoff schemes as described in the data paper (Steingroever et al., 2015): 1 = the original IGT, 2 = the variant with
# a constant loss on deck C and shuffled outcomes, 3 = the Bechara and Damasio scheme in which the decks change over blocks.
scheme_of <- c(Horstmann = 2, Kjome = 3, Maia = 1, Premkumar = 3, SteingroverInPrep = 2,
               Steingroever2011 = 2, Wetzels = 2, Wood = 3, Worthy = 1)

load_set <- function(name) {
  if (name == "simulated") return(readRDS("output/prereg/null_sim_data.rds"))
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

# Which participants of the replication sample a data set's analysis uses (all of them, or the 342 in 100-trial studies)
keep_of <- function(name, study) if (name == "replication342") !study %in% c("Steingroever2011", "Wetzels") else rep(TRUE, length(study))
# Where a data set's cross-validation fits are (replication342 reuses the replication fits: they are fit per person)
cv_source <- function(name) if (name == "replication342") "replication" else name
