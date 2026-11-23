"""Writes the hierarchical Stan models (stan/hier/) from one copy of the likelihood functions, so every model and the
likelihood check (stan/hier/check_lp.stan) use identical code. Edit this file, not the .stan files, then rerun it.
Run from the project folder: python3 scripts/55_make_stan.py ."""
import sys, os

FUNCTIONS = r"""
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
"""

# name, description, parameters [(name, lower, upper) or (name, None, None) for unbounded], call
MODELS = [
    ("orl", "ORL (Haines, Vassileva & Ahn, 2018)",
     [("Arew", 0, 1), ("Apun", 0, 1), ("Kp", 0, 5), ("betaF", None, None), ("betaP", None, None)],
     "orl_nova_trials({d}, Arew[i], Apun[i], Kp[i], betaF[i], betaP[i], betaP[i], 1.0, 0.0, 1.0)"),
    ("orl_split", "ORL + split: separate perseverance after gains and after losses",
     [("Arew", 0, 1), ("Apun", 0, 1), ("Kp", 0, 5), ("betaF", None, None), ("betaPw", None, None), ("betaPl", None, None)],
     "orl_nova_trials({d}, Arew[i], Apun[i], Kp[i], betaF[i], betaPw[i], betaPl[i], 1.0, 0.0, 1.0)"),
    ("orl_alpha", "ORL + alpha: outcome curvature",
     [("Arew", 0, 1), ("Apun", 0, 1), ("Kp", 0, 5), ("betaF", None, None), ("betaP", None, None), ("alpha", 0, 2)],
     "orl_nova_trials({d}, Arew[i], Apun[i], Kp[i], betaF[i], betaP[i], betaP[i], alpha[i], 0.0, 1.0)"),
    ("orl_split_alpha", "ORL + split + alpha",
     [("Arew", 0, 1), ("Apun", 0, 1), ("Kp", 0, 5), ("betaF", None, None), ("betaPw", None, None), ("betaPl", None, None), ("alpha", 0, 2)],
     "orl_nova_trials({d}, Arew[i], Apun[i], Kp[i], betaF[i], betaPw[i], betaPl[i], alpha[i], 0.0, 1.0)"),
    ("nova", "NOVA (Novelty, Outcome curvature, Value decay, Alternation): ORL + alpha with a starting value Q0 and decay D of the unchosen decks (script 50, version hybrid)",
     [("Arew", 0, 1), ("Apun", 0, 1), ("Kp", 0, 5), ("betaF", None, None), ("betaP", None, None), ("alpha", 0, 2), ("Q0", None, None), ("D", 0, 1)],
     "orl_nova_trials({d}, Arew[i], Apun[i], Kp[i], betaF[i], betaP[i], betaP[i], alpha[i], Q0[i], D[i])"),
    ("vse", "VSE (Ligneul, 2019)",
     [("theta", 0, 1), ("Delta", 0, 1), ("alpha", 0, 1), ("phi", None, None), ("beta", 0, 5)],
     "vse_trials({d}, theta[i], Delta[i], alpha[i], phi[i], beta[i], 0.0)"),
    ("vse_split", "VSE-split: after a net loss the chosen deck's exploration weight resets to psi instead of 0",
     [("theta", 0, 1), ("Delta", 0, 1), ("alpha", 0, 1), ("phi", None, None), ("beta", 0, 5), ("psi", None, None)],
     "vse_trials({d}, theta[i], Delta[i], alpha[i], phi[i], beta[i], psi[i])"),
    ("eef", "EEF (Yang, Xie & Wang, 2025)",
     [("theta", 0, 1), ("lambda", 0, 1), ("phi", None, None), ("beta", 0, 5)],
     "eef_trials({d}, theta[i], lambda[i], phi[i], beta[i])"),
]

HEADER = """/*
  {desc}: hierarchical Bayesian version, for the IGT-Modeling project (Hinojosa, 2026).
  Written by scripts/55_make_stan.py (do not edit by hand; edit the generator so every model keeps the same likelihood).
  Structure after hBayesDM's igt_orl.stan (CCS-Lab; Ahn, Haines & Zhang, 2017): group means and spreads on a probit
  scale, non-centered person-level parameters, bounded parameters through Phi_approx, unbounded ones directly.
  The likelihood is identical to the R functions used for every maximum-likelihood fit (checked by
  scripts/55a_stan_check.R).
  train[i, t] = 1 for the trials that count toward the fit (all of them for a full fit; for cross-validation, all but
  the held-out fold). lp_test holds each held-out trial's log-probability, for held-out scoring.
*/
"""

DATA = """data {
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
"""


def model_file(name, desc, params, call):
    K = len(params)
    bounded = [k for k, (_, lo, hi) in enumerate(params) if lo is not None]
    unbounded = [k for k, (_, lo, hi) in enumerate(params) if lo is None]
    out = [HEADER.format(desc=desc), "functions {" + FUNCTIONS + "}\n", DATA]
    out.append("parameters {\n  vector[%d] mu_pr;                                // group means (probit scale for bounded parameters)\n"
               "  vector<lower=0>[%d] sigma;                       // group spreads\n" % (K, K))
    for p, _, _ in params:
        out.append(f"  vector[N] {p}_pr;\n")
    out.append("}\n")
    tp = ["transformed parameters {\n"]
    for k, (p, lo, hi) in enumerate(params):
        if lo is None:
            tp.append(f"  vector[N] {p} = mu_pr[{k + 1}] + sigma[{k + 1}] * {p}_pr;\n")
        else:
            scale = "" if (lo, hi) == (0, 1) else f"{lo} + {hi - lo} * " if lo != 0 else f"{hi} * "
            tp.append(f"  vector<lower={lo}, upper={hi}>[N] {p} = {scale}Phi_approx(mu_pr[{k + 1}] + sigma[{k + 1}] * {p}_pr);\n")
    tp.append("}\n")
    out.append("".join(tp))
    m = ["model {\n  mu_pr ~ normal(0, 1);\n"]
    if bounded:
        m.append("  sigma[{%s}] ~ normal(0, 0.2);\n" % ", ".join(str(k + 1) for k in bounded))
    if unbounded:
        m.append("  sigma[{%s}] ~ cauchy(0, 1);\n" % ", ".join(str(k + 1) for k in unbounded))
    for p, _, _ in params:
        m.append(f"  {p}_pr ~ std_normal();\n")
    d = "choice[i, 1:Tsubj[i]], win[i, 1:Tsubj[i]], loss[i, 1:Tsubj[i]], Tsubj[i]"
    m.append("  for (i in 1:N) {\n")
    m.append(f"    vector[Tsubj[i]] lp = {call.format(d=d)};\n")
    m.append("    target += dot_product(lp, to_vector(train[i, 1:Tsubj[i]]));\n  }\n}\n")
    out.append("".join(m))
    g = ["generated quantities {\n"]
    for k, (p, lo, hi) in enumerate(params):
        if lo is None:
            g.append(f"  real mu_{p} = mu_pr[{k + 1}];\n")
        else:
            scale = "" if (lo, hi) == (0, 1) else f"{lo} + {hi - lo} * " if lo != 0 else f"{hi} * "
            g.append(f"  real mu_{p} = {scale}Phi_approx(mu_pr[{k + 1}]);\n")
    g.append("  vector[N] log_lik;                               // each participant's log-likelihood of the training trials\n")
    g.append("  array[N, n_test] real lp_test;                   // log-probability of each held-out trial\n")
    g.append("  for (i in 1:N) {\n")
    g.append(f"    vector[Tsubj[i]] lp = {call.format(d=d)};\n")
    g.append("    log_lik[i] = dot_product(lp, to_vector(train[i, 1:Tsubj[i]]));\n")
    g.append("    for (n in 1:n_test) lp_test[i, n] = lp[test_t[i, n]];\n  }\n}\n")
    out.append("".join(g))
    return "".join(out)


CHECK = """/*
  Likelihood check for scripts/55a_stan_check.R: the log-probability of every trial for given parameters, from the
  same functions as the models in this folder. Run with the fixed_param sampler, one draw.
*/
functions {""" + FUNCTIONS + """}
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
"""

if __name__ == "__main__":
    out_dir = os.path.join(sys.argv[1], "stan", "hier")
    os.makedirs(out_dir, exist_ok=True)
    for name, desc, params, call in MODELS:
        open(os.path.join(out_dir, f"igt_{name}.stan"), "w").write(model_file(name, desc, params, call))
    open(os.path.join(out_dir, "check_lp.stan"), "w").write(CHECK)
    print("wrote", len(MODELS), "models and check_lp.stan to", out_dir)
