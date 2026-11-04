/*
    Adapted from hBayesDM's igt_orl.stan (CCS-Lab).
    Modified by Alejandro Hinojosa (2026): added a training mask (train) and
    per-trial log-likelihood for cross-validation; removed y_pred;
    perseveration split into separate traces after gains and after losses (betaPw, betaPl).
    
    hBayesDM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    hBayesDM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with hBayesDM.  If not, see <http://www.gnu.org/licenses/>.
*/

data {
  int<lower=1> N;
  int<lower=1> T;
  array[N] int<lower=1, upper=T> Tsubj;
  array[N, T] int choice;
  array[N, T] real outcome;
  array[N, T] real sign_out;
  array[N, T] int<lower=0, upper=1> train;
}
transformed data {
  vector[4] initV;
  initV = rep_vector(0.0, 4);
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[6] mu_pr;
  vector<lower=0>[6] sigma;
  
  // Subject-level raw parameters (for Matt trick)
  vector[N] Arew_pr;
  vector[N] Apun_pr;
  vector[N] K_pr;
  vector[N] betaF_pr;
  vector[N] betaPw_pr;
  vector[N] betaPl_pr;
}
transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0, upper=1>[N] Arew;
  vector<lower=0, upper=1>[N] Apun;
  vector<lower=0, upper=5>[N] K;
  vector[N] betaF;
  vector[N] betaPw;
  vector[N] betaPl;
  
  for (i in 1 : N) {
    Arew[i] = Phi_approx(mu_pr[1] + sigma[1] * Arew_pr[i]);
    Apun[i] = Phi_approx(mu_pr[2] + sigma[2] * Apun_pr[i]);
    K[i] = Phi_approx(mu_pr[3] + sigma[3] * K_pr[i]) * 5;
  }
  betaF = mu_pr[4] + sigma[4] * betaF_pr;
  betaPw = mu_pr[5] + sigma[5] * betaPw_pr;
  betaPl = mu_pr[6] + sigma[6] * betaPl_pr;
}
model {
  // Hyperparameters
  mu_pr ~ normal(0, 1);
  sigma[1 : 3] ~ normal(0, 0.2);
  sigma[4 : 6] ~ cauchy(0, 1.0);
  
  // individual parameters
  Arew_pr ~ normal(0, 1.0);
  Apun_pr ~ normal(0, 1.0);
  K_pr ~ normal(0, 1.0);
  betaF_pr ~ normal(0, 1.0);
  betaPw_pr ~ normal(0, 1.0);
  betaPl_pr ~ normal(0, 1.0);
  
  for (i in 1 : N) {
    // Define values
    vector[4] ef;
    vector[4] ev;
    vector[4] PEfreq_fic;
    vector[4] PEval_fic;
    vector[4] persW; // perseveration after gains
    vector[4] persL; // perseveration after losses
    vector[4] util;
    
    real PEval;
    real PEfreq;
    real efChosen;
    real evChosen;
    real K_tr;
    
    // Initialize values
    ef = initV;
    ev = initV;
    persW = initV;
    persL = initV;
    util = initV;
    K_tr = pow(3, K[i]) - 1;
    
    for (t in 1 : Tsubj[i]) {
      // softmax choice
    if (train[i, t] == 1) {
      choice[i, t] ~ categorical_logit(util);
    }
      
      // Prediction error
      PEval = outcome[i, t] - ev[choice[i, t]];
      PEfreq = sign_out[i, t] - ef[choice[i, t]];
      PEfreq_fic = -sign_out[i, t] / 3 - ef;
      
      // store chosen deck ev
      efChosen = ef[choice[i, t]];
      evChosen = ev[choice[i, t]];
      
      if (outcome[i, t] >= 0) {
        // Update ev for all decks
        ef += Apun[i] * PEfreq_fic;
        // Update chosendeck with stored value
        ef[choice[i, t]] = efChosen + Arew[i] * PEfreq;
        ev[choice[i, t]] = evChosen + Arew[i] * PEval;
      } else {
        // Update ev for all decks
        ef += Arew[i] * PEfreq_fic;
        // Update chosendeck with stored value
        ef[choice[i, t]] = efChosen + Apun[i] * PEfreq;
        ev[choice[i, t]] = evChosen + Apun[i] * PEval;
      }
      
        // Perseverance updating
      if (outcome[i, t] >= 0) {
        persW[choice[i, t]] = 1;
        persL[choice[i, t]] = 0;
      } else {
        persL[choice[i, t]] = 1;
        persW[choice[i, t]] = 0;
      }
        persW /= 1 + K_tr;
        persL /= 1 + K_tr;  
      // Utility of expected value and perseverance
      util = ev + ef * betaF[i] + persW * betaPw[i] + persL * betaPl[i];
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_Arew;
  real<lower=0, upper=1> mu_Apun;
  real<lower=0, upper=5> mu_K;
  real mu_betaF;
  real mu_betaPw;
  real mu_betaPl;
  
  // For log likelihood calculation
  array[N, T] real log_lik;
  

  mu_Arew = Phi_approx(mu_pr[1]);
  mu_Apun = Phi_approx(mu_pr[2]);
  mu_K = Phi_approx(mu_pr[3]) * 5;
  mu_betaF = mu_pr[4];
  mu_betaPw = mu_pr[5];
  mu_betaPl = mu_pr[6];
  
  {
    // local section, this saves time and space
    for (i in 1 : N) {
      // Define values
      vector[4] ef;
      vector[4] ev;
      vector[4] PEfreq_fic;
      vector[4] PEval_fic;
      vector[4] persW; // perseveration after gains
      vector[4] persL; // perseveration after losses
      vector[4] util;
      
      real PEval;
      real PEfreq;
      real efChosen;
      real evChosen;
      real K_tr;
      
      // Initialize values
      for (t in 1 : T) log_lik[i, t] = 0;
      ef = initV;
      ev = initV;
      persW = initV;
      persL = initV;
      util = initV;
      K_tr = pow(3, K[i]) - 1;
      
      for (t in 1 : Tsubj[i]) {
        // softmax choice
        log_lik[i, t] = categorical_logit_lpmf(choice[i, t] | util);
        
        // Prediction error
        PEval = outcome[i, t] - ev[choice[i, t]];
        PEfreq = sign_out[i, t] - ef[choice[i, t]];
        PEfreq_fic = -sign_out[i, t] / 3 - ef;
        
        // store chosen deck ev
        efChosen = ef[choice[i, t]];
        evChosen = ev[choice[i, t]];
        
        if (outcome[i, t] >= 0) {
          // Update ev for all decks
          ef += Apun[i] * PEfreq_fic;
          // Update chosendeck with stored value
          ef[choice[i, t]] = efChosen + Arew[i] * PEfreq;
          ev[choice[i, t]] = evChosen + Arew[i] * PEval;
        } else {
          // Update ev for all decks
          ef += Arew[i] * PEfreq_fic;
          // Update chosendeck with stored value
          ef[choice[i, t]] = efChosen + Apun[i] * PEfreq;
          ev[choice[i, t]] = evChosen + Apun[i] * PEval;
        }
        
        // Perseverance updating
        if (outcome[i, t] >= 0) {
          persW[choice[i, t]] = 1;
          persL[choice[i, t]] = 0;
        } else {
          persL[choice[i, t]] = 1;
          persW[choice[i, t]] = 0;
        }
        persW /= 1 + K_tr;
        persL /= 1 + K_tr;
        
        // Utility of expected value and perseverance
        util = ev + ef * betaF[i] + persW * betaPw[i] + persL * betaPl[i];
      }
    }
  }
}

