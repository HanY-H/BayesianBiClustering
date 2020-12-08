#include "cluster_unc.h"

#include <cmath>

#include <vector>
#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"


void ClusterGPU::CalcCaches() {
  log_pi_S1 = log(pi_S);
  log_pi_S0 = log(1 - pi_S);

  // omega
  omega = vector<double>(G);
  log_omega1 = vector<double>(G);
  log_omega0 = vector<double>(G);
  for (int g = 0; g < G; ++g) {
    double n0 = a_w, n1 = a_w;
    for (int c = 0; c < C; ++c) {
      n0 += Z0[c][g];
      n1 += Z1[c][g];
    }
    omega[g] = n1 / (n0 + n1);
    log_omega1[g] = log(omega[g]);
    log_omega0[g] = log(1 - omega[g]);
  }

  // buff_log_l0
  buff_log_l0 = vector<double>(G);
  for (int g = 0; g < G; ++g) {
    int n0 = 0, n1 = 0;
    for (int c = 0; c < C; ++c) {
      n0 += Z0[c][g];
      n1 += Z1[c][g];
    }
    buff_log_l0[g] = log_pi_S0 + logBinomML(a_w, a_w, n1, n0);
  }

  // buff_log_binom_ML_theta
  buff_log_binom_ML_theta = matrix<double>(C + 1, C + 1);
  for (int n1 = 0; n1 <= C; ++n1)
    for (int n0 = 0; n0 <= C; ++n0)
      buff_log_binom_ML_theta[n1][n0] = logBinomML(a_t, a_t, n1, n0);

  // tmp variables
  pp = vector<double>(K, -10000000000);
  qq = vector<double>(K, -10000000000);
}

void ClusterGPU::InitMCMC() {
  // init MCMC states
  I = as<vector<int> >(static_cast<IntegerVector>(sample(K, C, true) - 1));
  S0 = vector<bool>(G, 0);
  S1 = vector<bool>(G, 1);
  theta = matrix<double>(G, K);
  log_theta1 = matrix<double>(G, K);
  log_theta0 = matrix<double>(G, K);
  for (int g = 0; g < G; ++g)
    for (int k = 0; k < K; ++k) {
      theta[g][k] = R::rbeta(a_t, a_t);
      log_theta1[g][k] = log(theta[g][k]);
      log_theta0[g][k] = log(1 - theta[g][k]);
    }

  pm_ll = -10000000000;
  t_I  = vector<vector<int> >(n_m);
  t_ll = vector<double>(n_b + n_m, 0);
}

void ClusterGPU::ClusterCells(int chain_id) {
  cout << "  Chain " << (chain_id + 1) << ":" << endl;

  for (m = 0; m < (n_b + n_m); ++m) {
    for (int c = 0; c < C; ++c) {
      if (verbose) cout << "ClusterCells::{m = " << m << " / " << (n_b + n_m) << ", c = " << c << endl;
      Update_I_c(c);
    }

    for (int g = 0; g < G; ++g)
      Update_S_g(g);

    for (int g = 0; g < G; ++g)
      Update_theta_g(g);

    SaveStates();
  }
  cout << "    Maximum posterior probability = " << pm_ll << ((pm_ll > ppm_ll) ? string(" (new max)") : string()) << endl;
}

// updating sub-population assignments under finite mixture model
void ClusterGPU::Update_I_c(int c) {
  // initialization of classifying probabilities
  for (int k = 0; k < K; ++k)
    pp[k] = log_f_Z_Ic(c, k);

  double MAX = -10000000000;

  for (int k = 0; k < K; ++k)
    if (pp[k] > MAX)
      MAX = pp[k];

  double sum = 0;
  for (int k = 0; k < K; ++k) {
    qq[k] = fast_exp(pp[k] - MAX);
    sum = sum + qq[k];
  }

  for (int k = 0; k < K; ++k)
    qq[k] = qq[k] / sum;

  for (int k = 1; k < K; ++k)
    qq[k] = qq[k] + qq[k - 1];

  double rand = R::runif(0, 1);
  int kk;
  for (kk = 0; kk < K; ++kk)
    if (rand <= qq[kk])
      break;

  if (kk < K) I[c] = kk;

//  if (!(kk >= 0 && kk < K)) {
//    cout << "Update_I_c::pp = ";
//    ctnutils::DispVector(pp);
//
//    for (int k = 0; k < K; ++k)
//      if (pp[k] > MAX)
//        MAX = pp[k];
//    cout << "Update_I_c::MAX = " << MAX << endl;
//
//    double sum = 0;
//    for (int k = 0; k < K; ++k) {
//      qq[k] = fast_exp(pp[k] - MAX);
//      sum = sum + qq[k];
//    }
//
//    cout << "Update_I_c::qq1 = ";
//    ctnutils::DispVector(qq);
//    for (int k = 0; k < K; ++k)
//      qq[k] = qq[k] / sum;
//
//    cout << "Update_I_c::qq2 = ";
//    ctnutils::DispVector(qq);
//
//    cout << endl;
//  }
}

void ClusterGPU::Update_S_g(int g) {
  double l0 = buff_log_l0[g];
  double l1 = 0;
  for (int k = 0; k < K; ++k) {
    double n0 = 0, n1 = 0;
    for (int c = 0; c < C; ++c) {
      n0 += (I[c] == k) * Z0[c][g];
      n1 += (I[c] == k) * Z1[c][g];
    }
    l1 += buff_log_binom_ML_theta[n1][n0];
  }
  double MAX = std::max(l0, l1);
  l0 = fast_exp(l0 - MAX);
  l1 = fast_exp(l1 - MAX);
  S1[g] = R::runif(0, 1) < (l1 / (l0 + l1));
  S0[g] = 1 - S1[g];
  // cout << "g = " << g << ",\tl1 - l0 = " << (l1 - l0) << ", S1[g] = " << S1[g] << endl;
}

void ClusterGPU::Update_theta_g(int g) {
  for (int k = 0; k < K; ++k) {
    double n0 = 0, n1 = 0;
    for (int c = 0; c < C; ++c) {
      n0 += (I[c] == k) * Z0[c][g];
      n1 += (I[c] == k) * Z1[c][g];
    }
    // theta[g][k] = (a_t + n1) / (2 * a_t + n0 + n1);
    theta[g][k] = R::rbeta(a_t + n1, a_t + n0);
    log_theta1[g][k] = log(theta[g][k]);
    log_theta0[g][k] = log(1 - theta[g][k]);
    // cout << "k = " << k << ", a_t = " << a_t << ", n1 = " << n1 << ", n0 = " << n0 << ", theta_g_k = " << theta[g][k] << endl;
  }
}

void ClusterGPU::SaveStates() {
  // show states
  int n_total = n_b + n_m;
  if (((m + 1) * 10) >= n_total && ((m + 1) * 10) % n_total == 0) {
    int pct = (m + 1) * 100 / n_total;
    if (pct % 20 == 0) {
      int nS = 0;
      for (int g = 0; g < G; ++g) nS += S1[g];
      int curr_ll = static_cast<int>(log_f_Z_I());
      cout << "    Finished " << numutils::EnoughSpaces(3, pct) << pct << "%, log-likelihood = " << curr_ll << ", # selected genes = " << nS << ", N(I) = {";
      for (int k = 0; k < K; ++k) {
        int n_k = 0;
        for (int c = 0; c < C; ++c) n_k += I[c] == k;
        cout << (k + 1) << ": " << n_k;
        if (k != K - 1) cout << ", ";
      }
      cout << "}" << endl;
    }
  }

  // calculate marginal likelihood
  t_ll[m] = log_f_Z_I();

  // save states
  if (m >= n_b) {
    int mm = m - n_b;
    t_I[mm] = I;
    if (t_ll[m] > pm_ll) {
      pm_ll = t_ll[m];
      pm_I  = I;
    }
  }
}

void ClusterGPU::GetPosteriorMode() {
  ll_traces.push_back(NumericVector(t_ll.begin(), t_ll.end()));
  if (pm_ll > ppm_ll) {
    ppm_I = pm_I;
    ppm_ll = pm_ll;
    ppm_t_I  = t_I;
    ppm_t_ll = t_ll;
  }
}

void ClusterGPU::PostProcess() {
  I = ppm_I;
  t_I = ppm_t_I;

  CalcPostProbs();
  CalcMargLikMethod1();

  cout << "Maximum posterior probability from multiple chains = " << ppm_ll << endl;
  cout << "Calculating marginal likelihood log{P(Z|M_K)} = " << log_p_Z_K << endl;
}

 double ClusterGPU::logBinomML(double a_p, double b_p, double a_c, double b_c) {
  return R::lbeta(a_p + a_c, b_p + b_c) - R::lbeta(a_p, b_p);
}

double ClusterGPU::log_f_Z_Ic(int c, int k) {
  double ll = 0;
  for (int g = 0; g < G; ++g) {
//    cout << "c = " << c << ", k = " << k << ", S0 = " << S0[g] << ", S1 = " << S1[g] << endl;
    ll += S0[g] * (Z0[c][g] * log_omega0[g] + Z1[c][g] * log_omega1[g]);
    ll += S1[g] * (Z0[c][g] * log_theta0[g][k] + Z1[c][g] * log_theta1[g][k]);
  }
  return ll;
}

double ClusterGPU::log_f_Z_I() {
  double ll = 0;
  for (int g = 0; g < G; ++g) {
    double l0 = buff_log_l0[g];
    double l1 = log_pi_S1;
    for (int k = 0; k < K; ++k) {
      int n1 = 0, n0 = 0;
      for (int c = 0; c < C; ++c) {
        n0 += (I[c] == k) * Z0[c][g];
        n1 += (I[c] == k) * Z1[c][g];
      }
      l1 += buff_log_binom_ML_theta[n1][n0];
      // cout << "k = " << k << ", n0 = " << N_k_g[0] << ", n1 = " << N_k_g[1] << ", ll = " << buff_log_binom_ML_theta[N_k_g[1]][N_k_g[0]] << endl;
    }
    // cout << "l0 = " << l0 << ", l1 = " << l1 << endl;
    ll += log_add(l0, l1);
  }
  return ll;
}

double ClusterGPU::log_f_Sg_I(int g) {
  double l0 = buff_log_l0[g];
  double l1 = log_pi_S1;
  for (int k = 0; k < K; ++k){
    double n0 = 0, n1 = 0;
    for (int c = 0; c < C; ++c) {
      n0 += (I[c] == k) * Z0[c][g];
      n1 += (I[c] == k) * Z1[c][g];
    }
    l1 += buff_log_binom_ML_theta[n1][n0];
  }
  return l1 - l0;
}

inline double ClusterGPU::fast_exp(double x) {
  if (x < -15) return 0;
  x = 1.0 + x / 32.0;
  x = x * x; x = x * x; x = x * x; x = x * x; x = x * x;
  return x;
}

inline double ClusterGPU::fast_log1p(double x) {
  return 0.9404 * x - 0.2537 * x * x;
}

inline double ClusterGPU::log_add(double a, double b) {
  if (a >= b) {
    return a + fast_log1p(fast_exp(b - a)); //log(1 + fast_exp(b - a));
  } else {
    return b + fast_log1p(fast_exp(a - b));
  }
}

void ClusterGPU::CalcPostProbs() {
  // p_S
  p_S = NumericVector(G);
  for (int g = 0; g < G; ++g)
    p_S[g] = calc_p_from_LLR(log_f_Sg_I(g));

  // p_I
  p_I = NumericMatrix(C, K);
  for (int m = 0; m < n_m; ++m)
    for (int c = 0; c < C; ++c)
      p_I(c, ppm_t_I[m][c]) += 1.0 / n_m;

  // p_I
  I_hat = IntegerVector(ppm_I.begin(), ppm_I.end());
}

void ClusterGPU::CalcMargLikMethod1() {
  // 1. logP(Z|M_K, I)
  I = ppm_I;
  double term1 = log_f_Z_I();
  // 2. logP(I|M_K)
  double term2 = C * log(1.0 / K);
  // 3. logP(I)
  double log_factorial = 0;
  for (int k = 1; k <= K; ++k)
    log_factorial += log(k);
  double n_samples = ppm_t_I.size();
  double prob_I = 0;
  for (int m = 0; m < n_samples; ++m)
    prob_I += (ppm_t_I[m] == I) / n_samples;
  double term3 = log(prob_I) - log_factorial;
  log_p_Z_I = term1;
  log_p_Z_K = term1 + term2 - term3;
  cout << "P(I_hat | Z) = " << prob_I << ", term1 = " << term1 << ", term2 = " << term2 << ", term3 = " << term3 << ", log_p_Z_K = " << log_p_Z_K << endl;
}

double ClusterGPU::calc_p_from_LLR(double LLR) {
  double mp = max(LLR, 0.0);
  double p0 = fast_exp(-mp);
  double p1 = fast_exp(LLR - mp);
  p1 = p1 / (p0 + p1);
  return p1;
}

pair<double, double> ClusterGPU::calc_log_pq_from_LLR(double LLR) {
  double p = calc_p_from_LLR(LLR);
  pair<double, double> log_pq;
  if (LLR >= 0) {
    log_pq.second = log(p);
    log_pq.first = log(p) - LLR;
  } else {
    log_pq.first = log(1 - p);
    log_pq.second = log(1 - p) + LLR;
  }
  return log_pq;
}
