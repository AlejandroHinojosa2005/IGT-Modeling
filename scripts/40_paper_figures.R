# Paper versions of the parameter recovery figures (scripts 06 and 13 made the originals).
# Uses the saved recovery results; nothing is refit. Sized to print at full page width
# (6.5 inches, 300 dpi), so the text comes out at about 9 points.
dir.create("output/figures/paper", showWarnings = FALSE, recursive = TRUE)
dot_col <- adjustcolor("#2F5E50", 0.5)   # pine, matching the paper's headings

# Parameter names as they appear in the paper's equations
labels <- list(
  pvl_delta = expression(A, alpha, c, lambda),
  vpp = expression(A, alpha, c, lambda, epsilon[P], epsilon[N], K, w),
  orl = expression(A[rew], A[pun], K * "'", beta[F], beta[P]),
  orl_wl = expression(A[rew], A[pun], K * "'", beta[F], beta[P * "," * win], beta[P * "," * loss]))
# Each parameter's fitting limits (scripts 02, 03, 04, 10); the axes show the full allowed range
limits <- list(
  pvl_delta = list(c(0, 1), c(0, 2), c(0, 5), c(0, 10)),
  vpp = list(c(0, 1), c(0, 2), c(0, 5), c(0, 10), c(-5, 5), c(-5, 5), c(0, 1), c(0, 1)),
  orl = list(c(0, 1), c(0, 1), c(0, 5), c(-10, 10), c(-10, 10)),
  orl_wl = list(c(0, 1), c(0, 1), c(0, 5), c(-10, 10), c(-10, 10), c(-10, 10)))
rows_cols <- list(pvl_delta = c(1, 4), vpp = c(2, 4), orl = c(2, 3), orl_wl = c(2, 3))
height <- c(pvl_delta = 2.1, vpp = 3.9, orl = 4.5, orl_wl = 4.5)   # inches

for (m in names(labels)) {
  x <- readRDS(sprintf("output/%s_recovery.rds", m))
  png(sprintf("output/figures/paper/recovery_%s.png", m), width = 6.5, height = height[[m]], units = "in", res = 300)
  if (m == "orl") {
    layout(matrix(c(1, 1, 2, 2, 3, 3,
                    0, 4, 4, 5, 5, 0), nrow = 2, byrow = TRUE))   # five panels: center the second row
  } else {
    par(mfrow = rows_cols[[m]])
  }
  par(mar = c(3.1, 3.3, 1.8, 0.6), mgp = c(2.1, 0.55, 0), cex = 0.75, las = 1, tcl = -0.3)
  for (k in seq_len(ncol(x$true))) {
    lim <- limits[[m]][[k]]                                         # same range on both axes
    plot(x$true[, k], x$recovered[, k], xlim = lim, ylim = lim, pch = 16, col = dot_col,
         xlab = "True", ylab = "Recovered")
    abline(0, 1, lty = 2, col = "grey40")                           # perfect recovery
    r <- sub("^(-?)0\\.", "\\1.", sprintf("%.2f", cor(x$true[, k], x$recovered[, k], use = "complete.obs")))
    title(main = bquote(.(labels[[m]][[k]]) * ",  " * italic(r) == .(r)), font.main = 1, line = 0.6)
  }
  dev.off()
}
