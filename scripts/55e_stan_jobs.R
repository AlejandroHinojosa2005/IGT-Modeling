# Lists the hierarchical fits (scripts/55_stan_fit.R) still to run for one data set, one per line, in the order
# scripts/55_queue.sh runs them: "<model> <fold>" for a first fit, "<model> <fold> long" for the rerun of a fit that
# did not pass the convergence rule (scripts/55_stan_common.R). Prints nothing when every fit is done.
# Order: the cross-validation folds of NOVA and ORL + alpha (predictions S1-S3 and S5), then the other models' folds
# (S4), then the fits to all trials (described only). A rerun takes its first fit's place in the order.
# Usage (from the project folder):
#   Rscript scripts/55e_stan_jobs.R <data>                the fits left
#   Rscript scripts/55e_stan_jobs.R summary <data> ...    what is done (scripts/55_queue.sh status)
source("scripts/55_stan_common.R")
args <- commandArgs(trailingOnly = TRUE)
run_order <- c("nova", "orl_alpha", "orl", "orl_split", "orl_split_alpha", "vse", "vse_split", "eef")
stopifnot(setequal(run_order, stan_models))
all_fits <- rbind(data.frame(model = rep(run_order, each = 10), fold = rep(as.character(1:10), length(run_order))),
                  data.frame(model = run_order, fold = "full"))

# The state of every fit of one data set
fit_state <- function(data) {
  st <- all_fits
  st$first_done <- st$passed <- st$rerun_done <- st$unusable <- FALSE
  st$minutes <- NA
  for (j in seq_len(nrow(st))) {
    fold <- if (st$fold[j] == "full") NA else as.integer(st$fold[j])
    first <- stan_file(data, st$model[j], fold)
    if (!file.exists(first)) next
    dg <- readRDS(first)$diagnostics
    st$first_done[j] <- TRUE
    st$passed[j] <- stan_pass(dg)
    st$minutes[j] <- dg[["minutes"]]
    rerun <- stan_file(data, st$model[j], fold, long = TRUE)
    if (file.exists(rerun)) {
      st$rerun_done[j] <- TRUE
      dg <- readRDS(rerun)$diagnostics
    }
    st$unusable[j] <- !st$passed[j] && stan_unusable(dg)
  }
  st
}

if (length(args) >= 1 && args[1] == "summary") {
  for (data in args[-1]) {
    st <- fit_state(data)
    cv <- st$fold != "full"
    need <- st$first_done & !st$passed
    cat(sprintf("%s: %d of 80 cross-validation fits and %d of 8 fits to all trials done; %d did not pass, %d of them rerun%s\n",
                data, sum(st$first_done & cv), sum(st$first_done & !cv), sum(need), sum(need & st$rerun_done),
                if (any(st$unusable & st$rerun_done)) sprintf(" (%d still failing after the rerun)", sum(st$unusable & st$rerun_done)) else ""))
    done <- st[st$first_done & cv, ]
    if (nrow(done)) {
      mins <- tapply(done$minutes, done$model, function(x) round(median(x, na.rm = TRUE)))
      left <- tapply(!st$first_done[cv], st$model[cv], sum)
      cat("  folds done and median minutes per fit, by model:",
          paste(sprintf("%s %d/10 (%s min)", names(left), 10 - left, ifelse(is.na(mins[names(left)]), "-", mins[names(left)])),
                collapse = ", "), "\n")
    }
  }
} else {
  data <- args[1]
  st <- fit_state(data)
  todo <- ifelse(!st$first_done, paste(st$model, st$fold),
                 ifelse(!st$passed & !st$rerun_done, paste(st$model, st$fold, "long"), NA))
  todo <- todo[!is.na(todo)]
  if (length(todo)) writeLines(todo)
}
