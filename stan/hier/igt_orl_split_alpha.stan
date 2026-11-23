/*
  ORL + split + alpha: hierarchical Bayesian version, for the IGT-Modeling project (Hinojosa, 2026).
  Written by scripts/55_make_stan.py (do not edit by hand; edit the generator so every model keeps the same likelihood).
  Structure after hBayesDM's igt_orl.stan (CCS-Lab; Ahn, Haines & Zhang, 2017): group means and spreads on a probit
  scale, non-centered person-level parameters, bounded parameters through Phi_approx, unbounded ones directly.
  The likelihood is identical to the R functions used for every maximum-likelihood fit (checked by
  scripts/55a_stan_check.R).
  train[i, t] = 1 for the trials that count toward the fit (all of them for a full fit; for cross-validation, all but
  the held-out fold). lp_test holds each held-out trial's log-probability, for held-out scoring.
*/
functions {
  // Per-trial log-probability of each choice (names cannot end in _lp in Stan), exactly as the R functions of scripts 43b (orl_family_lp, vse_family_lp,
  // eef_lp) and 50 (hybrid_lp, net bending). Outcomes are divided by 100; win >= 0, loss <= 0.

  // ORL family and NOVA. Q0 = 0 and D = 1 give ORL + split + alpha; alpha = 1 also gives ORL + split;
  // betaPw = betaPl gives ORL (+ alpha). Kp is K' (K = 3^K' - 1).
  vector orl_nova_trials(array[] int choice, array[] real win, array[] real loss, int Tn,
                     real Arew, real Apun, real Kp, real betaF, real betaPw, real betaPl,
                     real alpha, real Q0, real D) {
    vector[Tn] lp;
    vector[4] ev = rep_vector(Q0, 4);
    vector[4] ef = rep_vector(0, 4);
    vector[4] psw = rep_vector(0, 4);
    vector[4] psl = rep_vector(0, 4);
    real K = pow(3, Kp) - 1;
    for (t in 1:Tn) {
      int j = choice[t];
      lp[t] = categorical_logit_lpmf(j | ev + ef * betaF + psw * betaPw + psl * betaPl);
      real x = (win[t] + loss[t]) / 100;
      real sgn = x > 0 ? 1.0 : (x < 0 ? -1.0 : 0.0);
      real u = sgn * pow(abs(x), alpha);
      real ac = x >= 0 ? Arew : Apun;
      real au = x >= 0 ? Apun : Arew;
      for (k in 1:4) {
        if (k != j) {
          ev[k] = D * ev[k];
          ef[k] += au * (-sgn / 3 - ef[k]);
        }
      }
      ev[j] += ac * (u - ev[j]);
      ef[j] += ac * (sgn - ef[j]);
      psw /= 1 + K;
      psl /= 1 + K;
      if (x >= 0) {
        psw[j] = 1 / (1 + K);
        psl[j] = 0;
      } else {
        psl[j] = 1 / (1 + K);
        psw[j] = 0;
      }
    }
    return lp;
  }

  // VSE (psi = 0) and VSE-split. Exploration weights start at phi (the igt-toolbox convention); consistency
  // C = 3^beta - 1; wins and losses bent separately.
  vector vse_trials(array[] int choice, array[] real win, array[] real loss, int Tn,
                real theta, real Delta, real alpha, real phi, real beta, real psi) {
    vector[Tn] lp;
    vector[4] exploit = rep_vector(0, 4);
    vector[4] explore = rep_vector(phi, 4);
    real C = pow(3, beta) - 1;
    for (t in 1:Tn) {
      int j = choice[t];
      lp[t] = categorical_logit_lpmf(j | (exploit + explore) * C);
      real g = win[t] / 100;
      real l = abs(loss[t]) / 100;
      real v = (g > 0 ? pow(g, theta) : 0.0) - (l > 0 ? pow(l, theta) : 0.0);
      exploit *= Delta;
      exploit[j] += v;
      for (k in 1:4) {
        if (k != j) explore[k] += alpha * (phi - explore[k]);
      }
      explore[j] = (g - l < 0) ? psi : 0.0;
    }
    return lp;
  }

  // EEF (Yang, Xie & Wang, 2025). Exploitation weights start at w_ini, as in script 31.
  vector eef_trials(array[] int choice, array[] real win, array[] real loss, int Tn,
                real theta, real lambda, real phi, real beta) {
    vector[Tn] lp;
    vector[4] exploit = [0.0184, 0.0020, -0.0050, -0.0155]';
    vector[4] explore = rep_vector(0, 4);
    real C = pow(3, beta) - 1;
    for (t in 1:Tn) {
      int j = choice[t];
      lp[t] = categorical_logit_lpmf(j | (exploit + explore) * C);
      real g = win[t] / 100;
      real l = abs(loss[t]) / 100;
      real v = (g > 0 ? pow(g, theta) : 0.0) - (l > 0 ? pow(l, theta) : 0.0);
      exploit *= 1 - lambda;
      exploit[j] += v;
      for (k in 1:4) {
        if (k != j) explore[k] = lambda * explore[k] + (1 - lambda) * phi;
      }
      explore[j] = 0;
    }
    return lp;
  }
}
data {
  int<lower=1> N;                                  // participants
  int<lower=1> T;                                  // trials (the longest session)
  array[N] int<lower=1, upper=T> Tsubj;            // trials per participant
  array[N, T] int<lower=1, upper=4> choice;        // decks chosen (padding past Tsubj: any deck)
  array[N, T] real win;                            // each trial's gain (>= 0)
  array[N, T] real loss;                           // each trial's loss (<= 0)
  array[N, T] int<lower=0, upper=1> train;         // 1 = the trial counts toward the fit
  int<lower=1> n_test;                             // held-out trials per participant (1 for a full fit)
  array[N, n_test] int<lower=1, upper=T> test_t;   // which trials
}
parameters {
  vector[7] mu_pr;                                // group means (probit scale for bounded parameters)
  vector<lower=0>[7] sigma;                       // group spreads
  vector[N] Arew_pr;
  vector[N] Apun_pr;
  vector[N] Kp_pr;
  vector[N] betaF_pr;
  vector[N] betaPw_pr;
  vector[N] betaPl_pr;
  vector[N] alpha_pr;
}
transformed parameters {
  vector<lower=0, upper=1>[N] Arew = Phi_approx(mu_pr[1] + sigma[1] * Arew_pr);
  vector<lower=0, upper=1>[N] Apun = Phi_approx(mu_pr[2] + sigma[2] * Apun_pr);
  vector<lower=0, upper=5>[N] Kp = 5 * Phi_approx(mu_pr[3] + sigma[3] * Kp_pr);
  vector[N] betaF = mu_pr[4] + sigma[4] * betaF_pr;
  vector[N] betaPw = mu_pr[5] + sigma[5] * betaPw_pr;
  vector[N] betaPl = mu_pr[6] + sigma[6] * betaPl_pr;
  vector<lower=0, upper=2>[N] alpha = 2 * Phi_approx(mu_pr[7] + sigma[7] * alpha_pr);
}
model {
  mu_pr ~ normal(0, 1);
  sigma[{1, 2, 3, 7}] ~ normal(0, 0.2);
  sigma[{4, 5, 6}] ~ cauchy(0, 1);
  Arew_pr ~ std_normal();
  Apun_pr ~ std_normal();
  Kp_pr ~ std_normal();
  betaF_pr ~ std_normal();
  betaPw_pr ~ std_normal();
  betaPl_pr ~ std_normal();
  alpha_pr ~ std_normal();
  for (i in 1:N) {
    vector[Tsubj[i]] lp = orl_nova_trials(choice[i, 1:Tsubj[i]], win[i, 1:Tsubj[i]], loss[i, 1:Tsubj[i]], Tsubj[i], Arew[i], Apun[i], Kp[i], betaF[i], betaPw[i], betaPl[i], alpha[i], 0.0, 1.0);
    target += dot_product(lp, to_vector(train[i, 1:Tsubj[i]]));
  }
}
generated quantities {
  real mu_Arew = Phi_approx(mu_pr[1]);
  real mu_Apun = Phi_approx(mu_pr[2]);
  real mu_Kp = 5 * Phi_approx(mu_pr[3]);
  real mu_betaF = mu_pr[4];
  real mu_betaPw = mu_pr[5];
  real mu_betaPl = mu_pr[6];
  real mu_alpha = 2 * Phi_approx(mu_pr[7]);
  vector[N] log_lik;                               // each participant's log-likelihood of the training trials
  array[N, n_test] real lp_test;                   // log-probability of each held-out trial
  for (i in 1:N) {
    vector[Tsubj[i]] lp = orl_nova_trials(choice[i, 1:Tsubj[i]], win[i, 1:Tsubj[i]], loss[i, 1:Tsubj[i]], Tsubj[i], Arew[i], Apun[i], Kp[i], betaF[i], betaPw[i], betaPl[i], alpha[i], 0.0, 1.0);
    log_lik[i] = dot_product(lp, to_vector(train[i, 1:Tsubj[i]]));
    for (n in 1:n_test) lp_test[i, n] = lp[test_t[i, n]];
  }
}
