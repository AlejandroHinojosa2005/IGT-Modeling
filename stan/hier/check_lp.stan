/*
  Likelihood check for scripts/55a_stan_check.R: the log-probability of every trial for given parameters, from the
  same functions as the models in this folder. Run with the fixed_param sampler, one draw.
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
  int<lower=1> N;
  int<lower=1> T;
  array[N, T] int<lower=1, upper=4> choice;
  array[N, T] real win;
  array[N, T] real loss;
  int<lower=1, upper=3> family;                    // 1 = ORL family and NOVA, 2 = VSE family, 3 = EEF
  array[N] vector[9] par;                          // parameters in the order of the function's arguments
}
generated quantities {
  array[N, T] real lp;
  for (i in 1:N) {
    vector[T] x;
    if (family == 1) x = orl_nova_trials(choice[i], win[i], loss[i], T, par[i, 1], par[i, 2], par[i, 3], par[i, 4], par[i, 5], par[i, 6], par[i, 7], par[i, 8], par[i, 9]);
    if (family == 2) x = vse_trials(choice[i], win[i], loss[i], T, par[i, 1], par[i, 2], par[i, 3], par[i, 4], par[i, 5], par[i, 6]);
    if (family == 3) x = eef_trials(choice[i], win[i], loss[i], T, par[i, 1], par[i, 2], par[i, 3], par[i, 4]);
    lp[i] = to_array_1d(x);
  }
}
