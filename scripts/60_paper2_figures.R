# The figures for paper 2 (Mixing and Switching Between Models of the Iowa Gambling Task), drawn from saved results:
# nothing is refit and no random numbers are drawn, so every run gives the same figures.
#   Figure 1  Seldon: held-out predictions, then three ways to combine them (a diagram, filled in with one participant)
#   Figure 2  The registered predictions: Horstmann (exploratory) next to the replication (registered)
#   Figure 3  The shared mix against ORL + alpha, person by person
#   Figure 4  Switching without switching: simulated players with and without switching
#   Figure 5  VSE's exploration bonus is ORL's perseverance with a negative weight
#   Figure 6  NOVA: more useful in the mix than alone
#   Figure 7  How predictable are choices? Our models against Centaur
# Each figure is 6.5 inches wide (full page width) at 300 dpi, with text at about 9 points, as script 40 does for
# paper 1. Greek letters and minus signs are drawn with R's plotmath, so they come out the same on any computer.
# Written to output/figures/paper2/ and copied to paper/figures/paper2/ for paper/paper2.qmd, which gives each figure its
# number, title and note (APA style puts these in the text, not in the image). The numbers each figure shows are printed
# as it is drawn, to check against the logs.
#
# Needs (all saved by earlier scripts; see REPLICATION.md):
#   output/prereg/stacking_*.rds and switching_*.rds for replication, simulated and sim_switch, stacking_horstmann.rds
#     and switching_sim_types.rds (scripts 46 and 47), sim_switch_data.rds (script 44b)
#   output/switching.rds, switching_personal.rds, drift_curved.rds (scripts 37, 38, 42b: the Horstmann results the
#     registration reports)
#   output/hybrid/score_horstmann.rds, nova_replication.rds (scripts 50c, 50d)
#   output/pooled/score_horstmann.rds, score_replication.rds (script 54b)
#   For Figure 7: output/forward/ (scripts 57 and 57b, the latter run with the Centaur scores), the Centaur scores
#     centaur/centaur70b_scores.csv (centaur/run_centaur.sh), and data-raw/IGTdata.rdata. The Centaur scores are not in
#     the repository (Psych-101's test set may not be redistributed); without them Figure 7 is skipped.
#
# Usage, from the project folder:
#   in a terminal:   Rscript scripts/60_paper2_figures.R
#   in Positron:     open the IGT-Modeling folder (File > Open Folder), then in the R console:
#                    source("scripts/60_paper2_figures.R")
#                    (the console starts in the project folder, so the paths below work as they are)
source("scripts/43_prereg_data.R")
out_dir <- "output/figures/paper2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Look and shared helpers --------------------------------------------------------------------------------------
ink <- "#33312E"; muted <- "#8C877F"; rule <- "#DAD6CE"; band <- "#F1EEE8"
brown <- "#6B4423"; brown_light <- "#C4A484"; brown_dark <- "#2E1C0E"; brick <- "#A9533F"
model_names <- c("orl", "variant", "orl_alpha", "variant_alpha", "vse", "vse_split", "eef", "wsls")
model_label <- list(orl = "ORL", variant = "Variant", orl_alpha = quote("ORL +" ~ alpha),
                    variant_alpha = quote("Variant +" ~ alpha), vse = "VSE", vse_split = "VSE-split", eef = "EEF",
                    wsls = "WSLS", nova = "NOVA")
labs <- function(x) as.expression(unname(x))                   # a list of labels (strings or plotmath) for mtext/axis
# Model families share a hue: ORL and its extensions in browns, VSE in olives, EEF in slate, WSLS in grey
model_col <- c(orl = "#E3CCA6", variant = "#BF9A68", orl_alpha = "#8A5A30", variant_alpha = "#4A2E18",
               vse = "#B7C189", vse_split = "#66763A", eef = "#5E7C8E", wsls = "#B3AEA6", nova = brick)

figure <- function(file, height, draw) {
  png(file.path(out_dir, file), width = 6.5, height = height, units = "in", res = 300)
  on.exit(dev.off())
  draw()
}
look <- function() par(cex = 0.75, las = 1, tcl = -0.3, mgp = c(2, 0.45, 0), col.axis = ink, col.lab = ink, fg = ink)
start <- function(mat = matrix(1), widths = rep(1, ncol(mat)), heights = rep(1, nrow(mat))) {
  layout(mat, widths = widths, heights = heights)             # layout() resets the text size, so set it after
  look()
}
panel_label <- function(letter, title = "", line = 0.6) {
  mtext(bquote(bold(.(letter)) ~ ~ .(title)), side = 3, line = line, adj = 0, cex = 0.8, col = ink)
}
# A signed number with a true minus sign, and an optional note after it, as plotmath
num <- function(x, digits = 2, after = NULL) {
  s <- formatC(abs(x), format = "f", digits = digits)
  e <- if (x < 0) bquote(symbol("-") * .(s)) else bquote("+" * .(s))
  if (is.null(after)) e else bquote(.(e) ~ ~ ~ .(after))
}
ticks <- function(at, digits = 0, lead0 = TRUE) {                # tick labels with true minus signs
  labs(lapply(at, function(v) {
    s <- formatC(abs(v), format = "f", digits = digits)
    if (!lead0) s <- sub("^0\\.", ".", s)
    if (v < 0) bquote(symbol("-") * .(s)) else s
  }))
}
# A comparison as registered (script 48): d_i = person i's summed held-out log-likelihood under a minus under b;
# SE = sd(d) x sqrt(N), z = sum(d) / SE
compare <- function(a, b = 0) {
  d <- a - b; n <- length(d)
  c(diff = sum(d), se = sd(d) * sqrt(n), z = sum(d) / (sd(d) * sqrt(n)), better = sum(d > 0), n = n, mean = mean(d))
}
per_of <- function(name) {
  a <- readRDS(sprintf("output/prereg/stacking_%s.rds", name))
  b <- readRDS(sprintf("output/prereg/switching_%s.rds", name))
  list(per = cbind(a$per, b$cv), st = a, sw = b)
}
# A forest of z values: Horstmann (open, above) and the replication (filled, below), each with z +/- 1.96 (the 95%
# interval of the difference, in SE units); labels on the left, the replication's numbers on the right
forest <- function(zh, zr, left, tests, xlim, threshold = 1.645, both_sides = FALSE, codes = NULL) {
  k <- length(zr); y <- k:1
  plot(NA, xlim = xlim, ylim = c(0.5, k + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
  abline(h = y, col = band, lwd = 9)
  abline(v = 0, col = rule)
  abline(v = if (both_sides) c(-threshold, threshold) else threshold, lty = 2, col = muted)
  off <- 0.17
  segments(zh - 1.96, y + off, zh + 1.96, y + off, col = brown_light, lwd = 1.4)
  points(zh, y + off, pch = 21, bg = "white", col = brown_light, cex = 1.1, lwd = 1.4)
  segments(zr - 1.96, y - off, zr + 1.96, y - off, col = brown, lwd = 1.4)
  points(zr, y - off, pch = 19, col = brown, cex = 1.1)
  at <- seq(ceiling(xlim[1] / 2) * 2, xlim[2], 2)
  axis(1, at = at, labels = ticks(at), cex.axis = 0.9)
  mtext(quote(italic(z) ~ "(with its 95% interval)"), side = 1, line = 1.7, cex = 0.75)
  mtext(labs(left), side = 2, at = y, line = 0.5, adj = 1, cex = 0.75)
  if (!is.null(codes)) mtext(codes, side = 2, at = y, line = par("mar")[2] - 0.3, adj = 0, cex = 0.75, font = 2)
  # the replication's mean gain per person and how many people it favoured, in two columns with their headings
  mtext(labs(lapply(tests[, "mean"], num)), side = 4, at = y, line = 0.5, adj = 0, cex = 0.75)
  mtext(sprintf("%d of %d", tests[, "better"], tests[, "n"]), side = 4, at = y, line = 4.1, adj = 0, cex = 0.75)
  mtext(c("Gain", "Better for"), side = 4, at = k + 0.62, line = c(0.5, 4.1), adj = 0, cex = 0.7, col = muted)
}
forest_key <- function(y) {
  legend(par("usr")[1], y, c("Horstmann (N = 162, exploratory)", "Replication (N = 440)"), pch = c(21, 19),
         col = c(brown_light, brown), pt.bg = "white", pt.lwd = 1.4, horiz = TRUE, bty = "n", cex = 0.85, xpd = NA,
         yjust = 0, x.intersp = 0.6, text.width = NA)
}

# ---- The data -----------------------------------------------------------------------------------------------------
repl <- per_of("replication")
# Horstmann, as reported in the registration: the mixes from script 46 run on horstmann (it reproduces script 36 and
# adds the hedged models of script 42), switching from script 37, personal switching from 38, drifting weights from 42b
hors <- readRDS("output/prereg/stacking_horstmann.rds")$per
h_sw <- readRDS("output/switching.rds")
h_ps <- readRDS("output/switching_personal.rds")$cv
h_dr <- readRDS("output/drift_curved.rds")$cv
stopifnot(max(abs(h_sw$cv[, "orl_alpha"] - hors[, "orl_alpha"])) < 1e-8,
          max(abs(h_ps[, "shrink_mix"] - hors[, "shrink_mix"])) < 1e-8)
hors <- cbind(hors, h_sw$cv[, c("blend", "switching", "shuffled")], personal_switching = h_ps[, "personal_switching"],
              drift_curve = h_dr[, "drift_curve"])

# ---- Figure 1: Seldon, the procedure ------------------------------------------------------------------------------
# A diagram of the procedure, filled in with one replication participant so that the bars and colours are real: the
# one with the median expected number of switches. The switching model fit to everyone (script 47's fit_all) is run
# forward and backward over each person's held-out predictions, giving the probability that each model is in charge
# of each choice given all of that person's choices, and the expected number of switches.
n_r <- nrow(repl$per)
L <- aperm(array(repl$st$lpd, c(99, n_r, 8)), c(2, 1, 3))     # people x trials (2-100) x models, as script 47
w_fit <- repl$sw$fit_all$w; delta <- repl$sw$fit_all$delta
smooth_hmm <- function(Li, w, delta) {                          # one person: trials x models
  Tn <- nrow(Li); E <- exp(Li - apply(Li, 1, max))
  a <- matrix(NA, Tn, 8); b <- matrix(1, Tn, 8); pred <- w
  for (t in 1:Tn) { x <- pred * E[t, ]; a[t, ] <- x / sum(x); pred <- (1 - delta) * a[t, ] + delta * w }
  for (t in (Tn - 1):1) { v <- E[t + 1, ] * b[t + 1, ]; x <- (1 - delta) * v + delta * sum(w * v); b[t, ] <- x / sum(x) }
  g <- a * b
  stay <- sapply(1:(Tn - 1), function(t) {                    # P(the same model is in charge of choices t and t + 1)
    v <- E[t + 1, ] * b[t + 1, ]
    sum(a[t, ] * ((1 - delta) + delta * w) * v) / sum(a[t, ] * ((1 - delta) * v + delta * sum(w * v)))
  })
  list(in_charge = g / rowSums(g), switches = sum(1 - stay))
}
hmm_all <- lapply(1:n_r, function(i) smooth_hmm(L[i, , ], w_fit, delta))
switches <- sapply(hmm_all, `[[`, "switches")
ex <- order(switches)[ceiling(n_r / 2)]
lam <- unique(repl$st$lambda_chosen)
stopifnot(length(lam) == 1)                                     # the same lambda in every fold
w_personal <- lam * repl$st$w_person[ex, ] + (1 - lam) * repl$st$w_group
fold_shown <- 3                                                 # the fold drawn as held out: trials 3, 13, ..., 93
stopifnot(all(repl$st$rows$fold[repl$st$rows$subj == ex & repl$st$rows$trial %in% seq(fold_shown, 100, 10)] == fold_shown))
P_shown <- exp(L[ex, seq(fold_shown, 100, 10) - 1, ])           # 10 choices x 8 models
in_charge_top <- max.col(hmm_all[[ex]]$in_charge, "first")      # the most probable model in charge of each choice
cat("\nFigure 1: participant", ex, "(", repl$st$study[ex], "), expected switches", round(switches[ex], 2),
    "(median", round(median(switches), 2), "); delta", round(delta, 3), "; lambda", lam, "\n")
print(round(rbind(shared = repl$st$w_group, personal = w_personal), 3))

figure("fig1_seldon.png", 3.9, function() {
  look(); par(mar = c(0, 0, 0, 0))
  plot.new(); plot.window(xlim = c(0, 130), ylim = c(1, 79), xaxs = "i", yaxs = "i")
  panel <- function(x0, y0, x1, y1) rect(x0, y0, x1, y1, col = "#F8F6F1", border = rule)
  head <- function(x, y, num, txt) text(x, y, bquote(bold(.(num)) ~ ~ .(txt)), adj = c(0, 1), cex = 0.98)
  say <- function(x, y, lines, cex = 0.8, grey_from = Inf) {   # lines of text, each a string or plotmath; the
    for (k in seq_along(lines))                                 # lines naming the tests in grey
      text(x, y - (k - 1) * 2.7, lines[[k]], adj = c(0, 1), cex = cex, col = if (k >= grey_from) "#6F6A63" else ink)
  }
  stack <- function(x0, x1, y0, y1, w) {                        # a bar split by the models' weights
    e <- x0 + (x1 - x0) * cumsum(c(0, w))
    rect(e[-9], y0, e[-1], y1, col = model_col[model_names], border = "white", lwd = 0.6)
  }
  # 1: hold out each fold in turn
  panel(1, 2, 35, 78)
  head(2.5, 76.5, "1", "Hold out choices")
  pitch <- 2.95; sz <- 2.5; gx <- 3.6; gy <- 68
  text(gx + (1:10 - 1) * pitch + sz / 2, gy + 1.6, 1:10, cex = 0.62, col = muted)
  for (r in 0:9) for (cc in 1:10) {
    t <- 10 * r + cc
    fill <- if (t == 1) "white" else if (cc == fold_shown) brown else "#E4D8C4"
    rect(gx + (cc - 1) * pitch, gy - r * pitch - sz, gx + (cc - 1) * pitch + sz, gy - r * pitch, col = fill, border = rule,
         lwd = 0.5)
  }
  say(2.5, 36, list("100 choices, 10 to a row; each", "column is a fold. Fit each model to", "nine folds, predict the tenth (dark),",
                    "and repeat for all ten."))
  # 2: every model predicts the held-out choices
  panel(39, 2, 76, 78)
  head(40.5, 76.5, "2", "Predict them")
  yr <- 68 - (0:7) * 5.3
  for (m in 1:8) {
    rect(40.8, yr[m] - 1.2, 42.3, yr[m] + 1.2, col = model_col[model_names[m]], border = NA)
    text(43.2, yr[m], model_label[[model_names[m]]], adj = c(0, 0.5), cex = 0.78)
    segments(55, yr[m] - 2, 75, yr[m] - 2, col = rule, lwd = 0.6)
    xb <- 55.3 + (0:9) * 2
    rect(xb, yr[m] - 2, xb + 1.4, yr[m] - 2 + 4 * P_shown[, m], col = brown, border = NA)
  }
  say(40.5, 26.5, list("Each bar: the probability a model", "gave the deck chosen, on the ten", "choices of the dark fold."))
  # 3: three ways to combine the predictions
  dstr <- sub("^0", "", sprintf("%.3f", delta))
  panel(80, 57, 129, 78)
  head(81.5, 76.5, "3a", "Shared mix")
  stack(82, 127, 68.5, 72, repl$st$w_group)
  say(81.5, 66.5, list("One set of weights for everyone.", quote("Tested against ORL +" ~ alpha ~ "(H1, H1b)")), grey_from = 2)
  panel(80, 33, 129, 54)
  head(81.5, 52.5, "3b", "Personal weights")
  stack(82, 127, 44.5, 48, w_personal)
  say(81.5, 42.5, list(bquote(.(sub("^0", "", lam)) %*% "own weights" + .(sub("^0", "", 1 - lam)) %*% "shared weights"),
                       "Tested against the shared mix (H2)"), grey_from = 2)
  panel(80, 2, 129, 30)
  head(81.5, 28.5, "3c", "Switching")
  xe <- seq(82, 127, length.out = 100)
  rect(xe[-100], 20.5, xe[-1], 24, col = model_col[model_names[in_charge_top]], border = NA)
  say(81.5, 18.5, list("One model in charge at a time (bar),",
                       bquote("redrawn after each choice with probability" ~ delta * "."),
                       "Tested against a fixed blend (H3), shuffled", "trials (H4) and drifting weights (H7);",
                       "personal switching (H5)"), grey_from = 3)
  # arrows
  arrows(35.4, 50, 38.6, 50, length = 0.05, col = ink, lwd = 1)
  for (y in c(67.5, 43.5, 16)) arrows(76.4, 50, 79.6, y, length = 0.05, col = ink, lwd = 1)
})

# ---- Figure 2: the registered predictions -------------------------------------------------------------------------
hyps <- list(
  H1  = list("shared_mix", "orl_alpha", quote("Shared mix vs. ORL +" ~ alpha)),
  H1b = list("shared_mix", "orl_alpha_hedged", quote("Shared mix vs. ORL +" ~ alpha * ", hedged")),
  H2  = list("shrink_mix", "shared_mix", "Personal weights vs. shared mix"),
  H3  = list("switching", "blend", "Switching vs. a fixed blend"),
  H4  = list("switching", "shuffled", "Switching vs. shuffled trials"),
  H5  = list("personal_switching", "shrink_mix", "Personal switching vs. weights"),
  H7  = list("switching", "drift_curve", "Switching vs. drifting weights"))
tests_h <- t(sapply(hyps, function(h) compare(hors[, h[[1]]], hors[, h[[2]]])))
tests_r <- t(sapply(hyps, function(h) compare(repl$per[, h[[1]]], repl$per[, h[[2]]])))
cat("\nFigure 2 (Horstmann, then replication):\n")
print(round(cbind(tests_h[, c("diff", "z", "better")], tests_r[, c("diff", "z", "better", "mean")]), 2))

figure("fig2_predictions.png", 3.5, function() {
  start()
  par(mar = c(3, 15, 3, 8.6))
  forest(tests_h[, "z"], tests_r[, "z"], lapply(hyps, `[[`, 3), tests_r, c(0, 15), codes = names(hyps))
  text(1.645, 0.45, quote(italic(z) == 1.645), pos = 4, cex = 0.8, col = muted, offset = 0.2)
  forest_key(7.75)
})

# ---- Figure 3: the shared mix against ORL + alpha, person by person ------------------------------------------------
# x: the probability ORL + alpha gave each choice (geometric mean over the person's 99 held-out choices);
# y: how many times as likely the shared mix made the choices (the same ratio for the mix over ORL + alpha), log scale
fig3 <- lapply(list(Horstmann = hors, Replication = repl$per), function(p)
  cbind(x = exp(p[, "orl_alpha"] / 99), ratio = exp((p[, "shared_mix"] - p[, "orl_alpha"]) / 99)))
cat("\nFigure 3: median ratio and range\n")
print(sapply(fig3, function(f) round(c(median = median(f[, "ratio"]), min = min(f[, "ratio"]), max = max(f[, "ratio"]),
                                       better = sum(f[, "ratio"] > 1)), 3)))

figure("fig3_by_person.png", 3.2, function() {
  start(matrix(1:2, 1))
  par(mar = c(3.1, 3.6, 2.6, 0.8))
  for (j in seq_along(fig3)) {
    f <- fig3[[j]]; better <- f[, "ratio"] > 1
    plot(NA, xlim = c(0.2, 0.9), ylim = c(0.92, 2), log = "y", axes = FALSE, xlab = "", ylab = "")
    abline(h = 1, col = muted, lty = 2)
    points(f[, "x"], f[, "ratio"], pch = 16, cex = 0.65,
           col = ifelse(better, adjustcolor(brown, 0.55), adjustcolor(muted, 0.7)))
    axis(1, at = seq(0.2, 0.9, 0.1), labels = sub("^0", "", sprintf("%.1f", seq(0.2, 0.9, 0.1))), cex.axis = 0.9)
    axis(2, at = c(0.95, 1, 1.1, 1.25, 1.5, 2), labels = c(".95", "1", "1.1", "1.25", "1.5", "2"), cex.axis = 0.9)
    mtext(quote("Probability ORL +" ~ alpha ~ "gave the choice made"), side = 1, line = 1.7, cex = 0.72)
    mtext("Times as likely under the mix", side = 2, line = 2.4, cex = 0.72, las = 0)
    text(0.9, 2, sprintf("Mix better for %d of %d", sum(better), length(better)), adj = c(1, 1), cex = 0.85)
    panel_label(LETTERS[j], sprintf("%s (N = %d)", names(fig3)[j], nrow(f)))
  }
})

# ---- Figure 4: switching without switching ------------------------------------------------------------------------
# A: switching gains per person (the comparisons of H6) in the real data and in simulated players who mix on every
#    choice and never switch (script 44) or who switch at the Horstmann rate (script 44b). B: the switching rate chosen
#    in each of the ten cross-validation groups (script 47), against the simulated players' true rate.
sim <- per_of("simulated"); sws <- per_of("sim_switch"); swt <- readRDS("output/prereg/switching_sim_types.rds")
gain_sets <- list(`Horstmann (162)` = hors, `Replication (440)` = repl$per, `Simulated, no switching (324)` = sim$per,
                  `Simulated, switching (162)` = sws$per)
set_col <- c(brown_light, brown, "#9AA7B0", "#4F6575")
gain_defs <- list(H6a = list("switching", "blend", "vs. a fixed blend"),
                  H6b = list("switching", "shuffled", "vs. shuffled trials"),
                  H6c = list("personal_switching", "shrink_mix", "personal switching"))
gains <- lapply(gain_defs, function(g) t(sapply(gain_sets, function(p) {
  d <- p[, g[[1]]] - p[, g[[2]]]; se <- sd(d) / sqrt(length(d))
  c(mean = mean(d), lo = mean(d) - 1.96 * se, hi = mean(d) + 1.96 * se)
})))
rates <- list(`Simulated, no switching` = list(sim$sw$delta_fit[, "switching"], 1),
              `Simulated, one model only` = list(swt$delta_fit[, "switching"], 0),
              `Simulated, switching` = list(sws$sw$delta_fit[, "switching"],
                                            readRDS("output/prereg/sim_switch_data.rds")$delta),
              Horstmann = list(h_sw$delta_fit[, "switching"], NA), Replication = list(repl$sw$delta_fit[, "switching"], NA))
rate_col <- c("#9AA7B0", "#7C8C98", "#4F6575", brown_light, brown)
cat("\nFigure 4: gains per person (mean, 95% interval)\n"); print(lapply(gains, round, 3))
print(sapply(rates, function(r) round(c(mean_found = mean(r[[1]]), true = r[[2]]), 3)))

figure("fig4_simulations.png", 4.4, function() {
  start(matrix(c(1, 2, 3, 4, 4, 4), 2, byrow = TRUE), widths = c(0.526, 0.237, 0.237), heights = c(1.2, 1))
  xr <- c(-0.15, 0.9)
  for (b in seq_along(gains)) {
    g <- gains[[b]]; k <- nrow(g); y <- k:1
    par(mar = c(3, if (b == 1) 13.4 else 0.8, 3.4, 0.8))
    plot(NA, xlim = xr, ylim = c(0.5, k + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    abline(h = y, col = band, lwd = 9); abline(v = 0, col = muted, lty = 2)
    segments(g[, "lo"], y, g[, "hi"], y, col = set_col, lwd = 1.6)
    points(g[, "mean"], y, pch = 19, col = set_col, cex = 1.05)
    axis(1, at = c(0, 0.4, 0.8), labels = c("0", ".4", ".8"), cex.axis = 0.9)
    if (b == 1) mtext(rownames(g), side = 2, at = y, line = 0.5, adj = 1, cex = 0.75)
    if (b == 2) mtext("Switching gain per person (log-likelihood)", side = 1, line = 1.7, cex = 0.75)
    if (b == 1) panel_label("A", "Switching gains: real data and simulated players", line = 2.1)
    mtext(bquote(bold(.(names(gains)[b])) ~ ~ .(gain_defs[[b]][[3]])), side = 3, line = 0.4, adj = 0, cex = 0.72)
  }
  par(mar = c(3.2, 13.4, 2.4, 0.8))
  k <- length(rates); y <- k:1
  plot(NA, xlim = log10(c(0.001, 1.6)), ylim = c(0.5, k + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
  abline(h = y, col = band, lwd = 9)
  for (j in seq_len(k)) {
    r <- rates[[j]]
    points(log10(pmax(r[[1]], 0.001)), y[j] + seq(-0.2, 0.2, length.out = length(r[[1]])), pch = 16, cex = 0.6,
           col = rate_col[j])
    if (!is.na(r[[2]])) {
      xt <- log10(max(r[[2]], 0.001))
      segments(xt, y[j] - 0.36, xt, y[j] + 0.36, col = ink, lwd = 1.6)
      if (r[[2]] == 0) text(xt, y[j], "true rate 0", pos = 4, cex = 0.75, offset = 0.3)
    }
  }
  mtext(names(rates), side = 2, at = y, line = 0.5, adj = 1, cex = 0.75)
  axis(1, at = log10(c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)), labels = c(".001", ".003", ".01", ".03", ".1", ".3", "1"),
       cex.axis = 0.9)
  mtext(quote("Switching rate" ~ delta ~ "chosen in each of the ten groups (log scale); black line: true rate"),
        side = 1, line = 1.7, cex = 0.75)
  panel_label("B", "The rate of switching found")
})

# ---- Figure 5: VSE's exploration bonus is ORL's perseverance with a negative weight --------------------------------
# A deck chosen l trials ago, against a deck not chosen for a long time:
#   VSE:  (3^c - 1)(Explore - phi) = -(3^c - 1) phi (1 - eta)^(l - 1)
#   ORL:  beta_P PS = beta_P (1 + K)^(-l), the same curve when 1 / (1 + K) = 1 - eta and beta_P = -(3^c - 1) phi (1 + K)
# (script 56 checks this on every choice of both samples). Shown in units of (3^c - 1) phi.
etas <- c(0.1, 0.3, 0.6)
figure("fig5_identity.png", 2.8, function() {
  start()
  par(mar = c(3.1, 3.8, 0.8, 13.2))
  lags <- 1:15; cols <- c(brown_light, brown, brown_dark)
  plot(NA, xlim = c(1, 15), ylim = c(-1.03, 0.03), xlab = "", ylab = "", axes = FALSE)
  abline(h = 0, col = rule)
  for (j in seq_along(etas)) {
    eta <- etas[j]; K <- eta / (1 - eta)
    lines(lags, -(1 - eta)^(lags - 1), col = adjustcolor(cols[j], 0.35), lwd = 7)        # VSE
    points(lags, -(1 + K) * (1 + K)^(-lags), pch = 19, cex = 0.6, col = cols[j])          # ORL
  }
  axis(1, at = c(1, 5, 10, 15), cex.axis = 0.9)
  axis(2, at = c(-1, -0.5, 0), labels = expression(-1, -.5, 0), cex.axis = 0.9)
  mtext("Trials since the deck was chosen", side = 1, line = 1.7, cex = 0.75)
  mtext("Pull toward the deck", side = 2, line = 2.4, cex = 0.75, las = 0)
  x0 <- grconvertX(grconvertX(par("usr")[2], "user", "inches") + 0.15, "inches", "user")
  legend(x0, 0.05, c("VSE: exploration bonus", "ORL: perseverance,", "    negative weight"), lwd = c(7, NA, NA),
         pch = c(NA, 19, NA), pt.cex = 0.7, col = c(adjustcolor(brown, 0.35), brown, NA), bty = "n", cex = 0.8, xpd = NA,
         seg.len = 1.6)
  legend(x0, -0.5, as.expression(lapply(etas, function(e) bquote(eta == .(sub("^0", "", e)) * "," ~ ~ K == .(sprintf("%.2f", e / (1 - e)))))),
         lwd = 2.5, col = cols, bty = "n", cex = 0.8, xpd = NA, seg.len = 1.6)
})

# ---- Figure 6: NOVA -----------------------------------------------------------------------------------------------
nh <- readRDS("output/hybrid/score_horstmann.rds"); nr <- readRDS("output/hybrid/nova_replication.rds")
ph <- readRDS("output/pooled/score_horstmann.rds")$per; pr <- readRDS("output/pooled/score_replication.rds")$per
# The codes are those of the predictions written before the replication was scored: P1-P4 in
# predictions/nova_replication.md, Q1-Q5 in predictions/pooled_replication.md
nova_tests <- list(
  list(quote("NOVA vs. ORL +" ~ alpha * ", each fit alone"), nh$per[, "hybrid"] - nh$per[, "orl_alpha"],
       nr$per[, "nova"] - nr$per[, "orl_alpha"]),
  list(quote("NOVA vs. ORL +" ~ alpha * ", group prior"), ph[, "pooled.nova"] - ph[, "pooled.orl_alpha"],
       pr[, "pooled.nova"] - pr[, "pooled.orl_alpha"]),
  list("Mix with NOVA vs. mix of eight", nh$per[, "eight_plus_hybrid"] - nh$per[, "shared_mix"],
       nr$per[, "nine"] - nr$per[, "eight"]),
  list("Five with NOVA vs. mix of eight", nh$per[, "five_hybrid"] - nh$per[, "shared_mix"],
       nr$per[, "five_nova"] - nr$per[, "eight"]),
  list("NOVA (group prior) vs. mix of eight", ph[, "pooled.nova"] - ph[, "eight_unpooled"],
       pr[, "pooled.nova"] - pr[, "eight_unpooled"]))
nt_h <- t(sapply(nova_tests, function(x) compare(x[[2]])))
nt_r <- t(sapply(nova_tests, function(x) compare(x[[3]])))
w9 <- rbind(Horstmann = nh$w_nine, Replication = nr$w_nine)
colnames(w9)[colnames(w9) == "hybrid"] <- "nova"
cat("\nFigure 6 (Horstmann, then replication):\n")
print(round(cbind(nt_h[, c("diff", "z")], nt_r[, c("diff", "z", "better", "mean")]), 2))
print(round(w9, 3))

figure("fig6_nova.png", 5.2, function() {
  start(matrix(1:2, 2), heights = c(1.15, 1))
  # A: weights in the mix of nine, one row per model, aligned with B
  par(mar = c(2.6, 15, 2.6, 8.6))
  ord <- c("nova", model_names)
  k <- length(ord); y <- k:1
  plot(NA, xlim = c(0, 0.4), ylim = c(0.45, k + 0.55), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  rect(0, y + 0.02, w9["Horstmann", ord], y + 0.36, col = brown_light, border = NA)
  rect(0, y - 0.36, w9["Replication", ord], y - 0.02, col = brown, border = NA)
  mtext(labs(model_label[ord]), side = 2, at = y, line = 0.5, adj = 1, cex = 0.75,
        col = ifelse(ord == "nova", brick, ink), font = ifelse(ord == "nova", 2, 1))
  mtext(sprintf("%s   %s", sub("^0", "", sprintf("%.2f", w9["Horstmann", ord])), sub("^0", "", sprintf("%.2f", w9["Replication", ord]))),
        side = 4, at = y, line = 0.5, adj = 0, cex = 0.75)
  mtext("Horstmann, replication", side = 4, at = k + 0.62, line = 0.5, adj = 0, cex = 0.7, col = muted)
  axis(1, at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c("0", ".1", ".2", ".3", ".4"), cex.axis = 0.9)
  mtext("Weight in the shared mix of nine", side = 1, line = 1.6, cex = 0.75)
  legend("bottomright", c("Horstmann", "Replication"), fill = c(brown_light, brown), border = NA, bty = "n", cex = 0.8)
  panel_label("A", "NOVA gets the largest weight")
  # B: the comparisons
  par(mar = c(3, 15, 3.4, 8.6))
  forest(nt_h[, "z"], nt_r[, "z"], lapply(nova_tests, `[[`, 1), nt_r, c(-11, 9), both_sides = TRUE,
         codes = c("P1", "Q2", "P2", "P3", "Q5"))
  panel_label("B", "NOVA alone and in the mix", line = 2)
  forest_key(5.62)
})

# ---- Figure 7: our models against Centaur --------------------------------------------------------------------------
# Forward predictions (each choice from the trials before it only; script 57) for the 104 participants Centaur never
# trained on: Psych-101's held-out participants from both samples, and Maia and McClelland's study, which is not in
# Psych-101 at all. A: the probability each predictor gave the choice made (geometric mean over choices, as script 57b).
# B: Centaur's advantage over NOVA by what happened on the trial before, averaged over people.
cz_file <- "centaur/centaur70b_scores.csv"
if (!file.exists(cz_file)) {
  message("Figure 7 skipped: ", cz_file, " not found (run centaur/run_centaur.sh; see centaur/RENTING.md)")
} else {
  cz <- read.csv(cz_file)
  show <- c("orl_alpha", "nova", "mix_eight", "centaur")
  groups <- list(); edge <- list()
  cells <- c("After a gain, stayed", "After a gain, switched", "After a loss, stayed", "After a loss, switched")
  for (smp in c("horstmann", "replication")) {
    per <- readRDS(sprintf("output/forward/compare_%s_centaur.rds", smp))$per
    who <- unique(cz[cz$sample == smp, c("subj", "split")])
    split_of <- setNames(who$split, who$subj)[rownames(per)]
    for (sp in c("test", "absent")) {
      k <- rownames(per)[split_of == sp]
      if (!length(k)) next
      name <- if (smp == "horstmann") "Horstmann, held out" else if (sp == "test") "Replication, held out" else
        "Maia & McClelland, not in Psych-101"
      groups[[name]] <- c(chance = 0.25, exp(colSums(per[k, show, drop = FALSE]) / (99 * length(k))), n = length(k))
    }
    # per-trial log-probabilities of NOVA (script 57) and Centaur, for the situations in B
    d <- load_set(smp)
    nova_lp <- do.call(rbind, lapply(Sys.glob(sprintf("output/forward/%s/nova_part*.rds", smp)), function(f) readRDS(f)$lp))
    for (s in rownames(per)) {
      i <- as.integer(s)
      ci <- cz[cz$sample == smp & cz$subj == i & cz$trial >= 2, ]
      ci <- ci[order(ci$trial), ]
      stopifnot(nrow(ci) == 99, all(ci$deck == d$choice[i, 2:100]))
      lost <- (d$win[i, ] + d$loss[i, ])[1:99] < 0
      stayed <- d$choice[i, 2:100] == d$choice[i, 1:99]
      cell <- factor(cells[1 + 2 * lost + !stayed], levels = cells)
      edge[[paste(smp, s)]] <- tapply(ci$lp_4 - nova_lp[s, ], cell, mean)
    }
  }
  G <- do.call(rbind, groups)
  E <- do.call(rbind, edge)
  edge_tab <- t(apply(E, 2, function(x) { x <- x[!is.na(x)]; se <- sd(x) / sqrt(length(x))
    c(mean = mean(x), lo = mean(x) - 1.96 * se, hi = mean(x) + 1.96 * se, people = length(x)) }))
  cat("\nFigure 7: probability of the choice made, and Centaur minus NOVA by situation\n")
  print(round(G, 3)); print(round(edge_tab, 3))
  group_col <- c(brown_light, brown, brown_dark)
  figure("fig7_centaur.png", 3.3, function() {
    start(matrix(1:2, 1), widths = c(0.54, 0.46))
    par(mar = c(4, 3.6, 2.6, 0.6))
    xs <- 1:5
    plot(NA, xlim = c(0.6, 5.4), ylim = c(0.2, 0.46), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
    abline(h = seq(0.25, 0.45, 0.05), col = band)
    for (j in seq_len(nrow(G))) {
      v <- G[j, c("chance", show)]
      lines(xs, v, col = group_col[j], lwd = 1.2)
      points(xs, v, pch = 19, col = group_col[j], cex = 0.9)
    }
    xl <- list("Chance", quote("ORL +" ~ alpha), "NOVA", "Mix of\neight", "Centaur")
    for (j in xs) mtext(xl[[j]], side = 1, at = j, line = 0.5, cex = 0.72, padj = 1)
    axis(2, at = seq(0.2, 0.45, 0.05), labels = sub("^0", "", sprintf("%.2f", seq(0.2, 0.45, 0.05))), cex.axis = 0.9)
    mtext("Probability given to the choice made", side = 2, line = 2.4, cex = 0.72, las = 0)
    legend("bottomright", sprintf("%s (%d)", rownames(G), G[, "n"]), col = group_col, pch = 19, lwd = 1.2, bty = "n",
           cex = 0.72, inset = c(0, 0.02))
    panel_label("A", "Each choice predicted from the ones before it")
    par(mar = c(4, 9, 2.6, 0.8))
    k <- nrow(edge_tab); y <- k:1
    plot(NA, xlim = c(-0.1, 0.9), ylim = c(0.5, k + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    abline(h = y, col = band, lwd = 9); abline(v = 0, col = muted, lty = 2)
    segments(edge_tab[, "lo"], y, edge_tab[, "hi"], y, col = brown, lwd = 1.6)
    points(edge_tab[, "mean"], y, pch = 19, col = brown)
    mtext(rownames(edge_tab), side = 2, at = y, line = 0.5, adj = 1, cex = 0.72)
    axis(1, at = c(0, 0.4, 0.8), labels = c("0", ".4", ".8"), cex.axis = 0.9)
    mtext("Centaur minus NOVA,\nlog-probability per choice", side = 1, line = 2.5, cex = 0.72)
    panel_label("B", "Where Centaur gains")
  })
}
# The paper (paper/paper2.qmd, not in the repository) includes the figures from paper/figures/paper2/, where each gets
# its number, title and note
if (dir.exists("paper")) {
  dir.create("paper/figures/paper2", showWarnings = FALSE, recursive = TRUE)
  file.copy(list.files(out_dir, full.names = TRUE), "paper/figures/paper2", overwrite = TRUE)
}
cat("\nFigures written to", out_dir, if (dir.exists("paper")) "and copied to paper/figures/paper2", "\n")
