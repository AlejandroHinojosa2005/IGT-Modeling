# Checks the correspondence used in paper 2 (and stated in script 50): VSE's exploration term is ORL's perseverance
# with a negative weight.
#
# In VSE (scripts 31 and 43b), after each choice the chosen deck's exploration value is reset (to 0; in VSE-split, to
# psi after a net loss) and every other deck's moves toward phi at rate eta (the parameter VSE calls alpha):
#   Explore[k] <- Explore[k] + eta * (phi - Explore[k])
# Writing G = phi - Explore: G is set to phi (or phi - psi after a loss) for the chosen deck, shrinks by (1 - eta) for
# every other deck, and starts at 0 (Explore starts at phi). ORL's perseverance PS is set to 1 / (1 + K) for the chosen
# deck (PSw after a gain, PSl after a loss, in the split version), shrinks by 1 / (1 + K) for every other deck, and
# starts at 0. So when 1 - eta = 1 / (1 + K),
#   G = (1 + K) * (phi * PSw + (phi - psi) * PSl),   and   C * Explore = C * phi - C * G,
# where C is VSE's consistency (3^beta - 1). C * phi is the same for every deck and cancels in the softmax, so VSE's
# exploration term is ORL's perseverance with weights
#   betaPw = -C * phi * (1 + K),   betaPl = C * (psi - phi) * (1 + K)
# (VSE has psi = 0, so betaPw = betaPl = -C * phi * (1 + K): one perseverance weight, negative when phi > 0.)
#
# Part 1: for random VSE and VSE-split parameters and every participant in both samples, VSE's log-probability of
# each choice (script 43b's function, unchanged) equals that of the same model with its exploration term replaced by
# ORL's perseverance (the lines of script 43b's ORL function) at the weights above. Differences should be about 1e-12.
# Part 2: the starting value matters. If every deck's exploration value starts at a value other than phi (here 1),
# the two differ from the second choice until every deck has been chosen once, and are identical from then on.
# Usage (from the project folder): Rscript scripts/56_identity_check.R
source("scripts/43_prereg_data.R")
source("scripts/43b_prereg_models.R")

# VSE with its exploration term replaced by ORL's perseverance; everything else is VSE's (script 43b)
vse_as_orl_lp <- function(par6, choice, win, loss) {
  theta <- par6[1]; Delta <- par6[2]; eta <- par6[3]; phi <- par6[4]; C <- 3^par6[5] - 1; psi <- par6[6]
  K <- 1 / (1 - eta) - 1                                  # so that 1 / (1 + K) = 1 - eta
  betaPw <- -C * phi * (1 + K)
  betaPl <- C * (psi - phi) * (1 + K)
  Exploit <- c(0, 0, 0, 0)
  PSw <- PSl <- c(0, 0, 0, 0)
  lp <- rep(NA, length(choice))
  for (t in seq_along(choice)) {
    j <- choice[t]
    V <- Exploit * C + PSw * betaPw + PSl * betaPl
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- Exploit * Delta
    Exploit[j] <- Exploit[j] + v
    x <- (win[t] + loss[t]) / 100                         # ORL's perseverance, as in script 43b
    PSw <- PSw / (1 + K)
    PSl <- PSl / (1 + K)
    if (x >= 0) { PSw[j] <- 1 / (1 + K); PSl[j] <- 0 } else { PSl[j] <- 1 / (1 + K); PSw[j] <- 0 }
  }
  lp
}

# VSE (script 43b's function) with the exploration values starting at `start` instead of phi
vse_start_lp <- function(par6, choice, win, loss, start) {
  theta <- par6[1]; Delta <- par6[2]; alpha <- par6[3]; phi <- par6[4]; C <- 3^par6[5] - 1; psi <- par6[6]
  Exploit <- c(0, 0, 0, 0)
  Explore <- rep(start, 4)
  lp <- rep(NA, length(choice))
  for (t in 1:length(choice)) {
    j <- choice[t]
    V <- (Exploit + Explore) * C
    lp[t] <- (V - max(V) - log(sum(exp(V - max(V)))))[j]
    g <- win[t] / 100
    l <- abs(loss[t]) / 100
    v <- (if (g > 0) g^theta else 0) - (if (l > 0) l^theta else 0)
    Exploit <- Exploit * Delta
    Exploit[j] <- Exploit[j] + v
    Explore[-j] <- Explore[-j] + alpha * (phi - Explore[-j])
    Explore[j] <- if (g - l < 0) psi else 0
  }
  lp
}

# Random parameters within VSE-split's limits (script 31); eta kept below .99 so that K is finite
draw <- function() { p <- runif(6, lower_vse_split, upper_vse_split); p[3] <- runif(1, 0, .99); p }

# Part 1
set.seed(56)
worst <- c(vse = 0, vse_split = 0, start_phi = 0)
n_people <- 0
for (data_name in c("horstmann", "replication")) {
  d <- load_set(data_name)
  for (i in 1:nrow(d$choice)) {
    ch <- d$choice[i, ]; w <- d$win[i, ]; l <- d$loss[i, ]
    p <- draw()
    worst["vse_split"] <- max(worst["vse_split"], abs(models$vse_split$lp(p, ch, w, l) - vse_as_orl_lp(p, ch, w, l)))
    worst["start_phi"] <- max(worst["start_phi"], abs(vse_start_lp(p, ch, w, l, start = p[4]) - models$vse_split$lp(p, ch, w, l)))
    p[6] <- 0                                             # VSE: psi = 0
    worst["vse"] <- max(worst["vse"], abs(models$vse$lp(p[1:5], ch, w, l) - vse_as_orl_lp(p, ch, w, l)))
    n_people <- n_people + 1
  }
}
cat("participants:", n_people, "\n")
signif(worst, 3)                                          # largest difference in log-probability (about 1e-12)
stopifnot(all(worst < 1e-8))
cat("VSE's exploration term equals ORL's perseverance at betaPw = -C phi (1 + K), betaPl = C (psi - phi) (1 + K)\n")

# Part 2: exploration values starting at 1 instead of phi
d <- load_set("horstmann")
set.seed(57)
before <- after <- 0
for (i in 1:nrow(d$choice)) {
  ch <- d$choice[i, ]; w <- d$win[i, ]; l <- d$loss[i, ]
  all4 <- which(sapply(1:100, function(t) length(unique(ch[1:t])) == 4))[1]   # the choice that completes the set
  if (is.na(all4) || all4 == 100) next
  p <- draw()
  diff <- abs(vse_start_lp(p, ch, w, l, start = 1) - vse_as_orl_lp(p, ch, w, l))
  before <- max(before, diff[1:all4])
  after <- max(after, diff[(all4 + 1):100])
}
signif(c(until_every_deck_chosen = before, afterwards = after), 3)
