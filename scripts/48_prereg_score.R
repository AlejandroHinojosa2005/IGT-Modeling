# Preregistered replication (https://doi.org/10.17605/OSF.IO/U86SY), step 6: scores every prediction as registered.
# For each comparison, d_i is participant i's summed held-out log-likelihood under the first model minus the second.
# Reported: sum of d, SE = sd(d) x sqrt(N), z = sum(d) / SE, and how many participants have d_i > 0.
# Each hypothesis is met if its one-sided z > 1.645 (alpha = .05); no correction for multiple tests.
# H6 compares the mean gain per participant in the real Horstmann data (scripts 37-38; SE = sd / sqrt(162)) with the
# mean in the simulated data (SE = sd / sqrt(324)): z = (m_real - m_sim) / sqrt(SE_real^2 + SE_sim^2).
# Then the robustness checks (registration, Q24) and the exploratory descriptions (Q25), all labeled.
# Needs scripts 46 and 47 run on replication, simulated and replication342.
# Usage (from the project folder): Rscript scripts/48_prereg_score.R
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

repl <- per_of("replication")
sim <- per_of("simulated")
model_names <- names(repl$st$w_group)

# 1. The predictions, in the replication sample (N = 440)
main <- test_all(repl$per)
show(main)

# 2. H6: switching gains per participant, real Horstmann data vs. simulated data with no switching
h_sw <- readRDS("output/switching.rds")$cv
h_ps <- readRDS("output/switching_personal.rds")$cv
gains <- function(sw, ps) cbind(H6a = sw[, "switching"] - sw[, "blend"], H6b = sw[, "switching"] - sw[, "shuffled"],
                                H6c = ps[, "personal_switching"] - ps[, "shrink_mix"])
real <- gains(h_sw, h_ps)
simg <- gains(sim$per, sim$per)
round(colSums(real), 1)                                  # the registered Horstmann gains: +99.6, +83.5, +32.3
h6 <- as.data.frame(t(sapply(colnames(real), function(k) {
  mr <- mean(real[, k]); sr <- sd(real[, k]) / sqrt(nrow(real))
  ms <- mean(simg[, k]); ss <- sd(simg[, k]) / sqrt(nrow(simg))
  c(m_real = mr, se_real = sr, m_sim = ms, se_sim = ss, z = (mr - ms) / sqrt(sr^2 + ss^2))
})))
h6$met <- h6$z > 1.645
round(h6[, 1:5], 3); h6$met

# 3. How the results are read (the rule registered under Hypotheses)
met <- c(setNames(main$met, rownames(main)), setNames(h6$met, rownames(h6)))
met
rule <- c(
  mixing = if (all(met[c("H1", "H1b", "H2")])) "Mixing (H1, H1b, H2) is supported." else "Mixing is not fully supported: not all of H1, H1b and H2 were met.",
  switching = if (all(met[c("H3", "H4", "H7", "H6a", "H6b")])) "The results are consistent with switching between strategies (H3, H4 and H7 met in the new data; H6a and H6b met)." else "The results are not called consistent with switching: not all of H3, H4, H7, H6a and H6b were met.",
  drift = if (!met["H7"]) "Drifting weights predicted as well as switching (H7 not met): the balance between strategies slides over the session rather than switching." else "",
  simulation = if (!all(met[c("H6a", "H6b")])) "The simulated data showed switching gains as large as the real data (H6a or H6b not met): the switching model captures how the models' relative accuracy changes over a session, not people switching strategies." else "",
  personal = if (all(met[c("H5", "H6c")])) "Personal switching: H5 and H6c met." else "Personal switching: not both of H5 and H6c met.")
writeLines(rule[rule != ""])

# 4. Robustness checks (registration, Q24)
# 4.1 Without participants who have any fold below -100 in any model (rule of Hinojosa, 2026)
lpd <- repl$st$lpd
fold_sums <- rowsum(lpd, interaction(repl$st$rows$subj, repl$st$rows$fold, lex.order = TRUE))
worst <- apply(array(fold_sums, c(10, nrow(repl$per), ncol(lpd))), 2, min)
bad <- which(worst < -100)
length(bad)                                              # participants left out
r1 <- test_all(repl$per[!seq_len(nrow(repl$per)) %in% bad, , drop = FALSE])
show(r1)
# 4.2 The mixes vs. the best single model in the replication sample, if that is not ORL + alpha
best_single <- names(which.max(colSums(repl$per[, model_names])))
best_single
r2 <- if (best_single != "orl_alpha") round(rbind(shared = compare(repl$per, "shared_mix", best_single),
                                                  personal = compare(repl$per, "personal_mix", best_single),
                                                  shrinkage = compare(repl$per, "shrink_mix", best_single)), 1)
r2
# 4.3 The shared mix vs. the best hedged single model, if that is not ORL + alpha
best_hedged <- names(which.max(colSums(repl$per[, paste0(model_names, "_hedged")])))
best_hedged
r3 <- if (best_hedged != "orl_alpha_hedged") round(compare(repl$per, "shared_mix", best_hedged), 1)
r3
# 4.4 Switching vs. the straight-line drift (k = 0), as in script 42
r4 <- round(compare(repl$per, "switching", "drift_line"), 1)
r4
# 4.5 Without the two 150-trial studies (N = 342): scripts 46-47 rerun on that subset
r5 <- test_all(per_of("replication342")$per)
show(r5)
# 4.6 By study and by payoff scheme, described without tests: sum of d and participants with d > 0
by_group <- function(g) lapply(hyps, function(h) {
  d <- repl$per[, h[1]] - repl$per[, h[2]]
  cbind(N = tapply(d, g, length), diff = round(tapply(d, g, sum), 1), better_for = tapply(d > 0, g, sum))
})
r6_study <- by_group(repl$st$study)
r6_scheme <- by_group(repl$st$scheme)
r6_study
r6_scheme

# 5. Exploratory (registration, Q25), described only
round(repl$st$w_group, 3)                                 # shared weights
round(colMeans(repl$st$w_person), 3)                      # personal weight profiles: average weight per model
table(factor(model_names[apply(repl$st$w_person, 1, which.max)], levels = model_names))   # each person's largest weight
repl$st$lambda_chosen                                     # shrinkage lambda chosen in each fold
round(repl$sw$delta_fit[, "switching"], 3)                # switching: delta chosen in each group
round(repl$sw$personal_fits, 3)                           # personal switching: delta and lambda in each group
round(rbind(early = repl$sw$drift_all$we, late = repl$sw$drift_all$wl), 3)   # drift fit on everyone
round(repl$sw$drift_all$k, 2)                             # its curve
round(repl$sw$k_fit, 2)                                   # k chosen in each group
round(rbind(drift_switching_vs_switching = compare(repl$per, "drift_switching", "switching"),
            drift_switching_vs_drift = compare(repl$per, "drift_switching", "drift_curve"),
            types_vs_blend = compare(repl$per, "types", "blend")), 1)

saveRDS(list(main = main, h6 = h6, real_gains = real, sim_gains = simg, met = met, rule = rule,
             robustness = list(left_out = bad, without_left_out = r1, best_single = best_single, mixes_vs_best = r2,
                               best_hedged = best_hedged, shared_vs_best_hedged = r3, switching_vs_line = r4,
                               without_150 = r5, by_study = r6_study, by_scheme = r6_scheme)),
        "output/prereg/results.rds")
