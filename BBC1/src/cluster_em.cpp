#include "cluster_em.h"

#include <cmath>

#include <vector>
#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"


void ClusterEM::CalcCaches() {
  // compute omega hat and buff_l0
  omega = vector<double>(G);
  log_omega1 = vector<double>(G);
  log_omega0 = vector<double>(G);
  buff_l0 = vector<double>(G);
  for (int g = 0; g < G; ++g) {
    vector<double> nn(2, a_w);
    for (int c = 0; c < C; ++c) nn[Z[c][g]] ++;
    omega[g] = nn[1] / (nn[0] + nn[1]);
    log_omega1[g] = log(omega[g]);
    log_omega0[g] = log(1 - omega[g]);
    for (int c = 0; c < C; ++c)
      buff_l0[g] += log_omega1[g] * Z[c][g] + log_omega0[g] * (1 - Z[c][g]);
  }

  // tmp variables
  pp = vector<double>(K, -10000000000);
  qq = vector<double>(K, -10000000000);
}

void ClusterEM::InitEM() {
  // init MCMC states
  S = vector<bool>(G, true);
  p_I = matrix<double>(C, K);
  ll_trace = vector<double>(n_steps);

  // set pi = (1/K, ..., 1/K)
  pi = vector<double>(K);
  log_pi = vector<double>(K);
  for (int k = 0; k < K; ++k) {
    pi[k] = 1.0 / K;
    log_pi[k] = -log(K);
  }
  // random thetas from Beta(a_t, a_t)
  theta = matrix<double>(G, K);
  log_theta1 = matrix<double>(G, K);
  log_theta0 = matrix<double>(G, K);
  for (int g = 0; g < G; ++g)
    for (int k = 0; k < K; ++k) {
      theta[g][k] = R::rbeta(a_t, a_t);
      log_theta1[g][k] = log(theta[g][k]);
      log_theta0[g][k] = log(1 - theta[g][k]);
    }
}

void ClusterEM::RunEM(int round_id) {
  if (!parallel_print)
    cout << "  Round " << (round_id + 1) << ":" << endl;

  for (m = 0; m < n_steps; ++m) {
    // E-step
    for (int c = 0; c < C; ++c)
      Update_p_I(c);
    for (int g = 0; g < G; ++g)
      Sample_S(g);

    // M-step
    Update_pi();
    for (int g = 0; g < G; ++g)
      Update_theta(g);

    SaveStates();
    ShowStates();
  }
}

void ClusterEM::Update_p_I(int c) {
  pp = log_pi;
  for (int g = 0; g < G; ++g) {
    if (S[g] == 0) continue;
    bool Zcg = Z[c][g];
    if (Zcg)
      for (int k = 0; k < K; ++k) pp[k] += log_theta1[g][k];
    else
      for (int k = 0; k < K; ++k) pp[k] += log_theta0[g][k];
  }

  double MAX = -1e9, SUM = 0;
  for (int k = 0; k < K; ++k)
    MAX = std::max(MAX, pp[k]);
  for (int k = 0; k < K; ++k) {
    pp[k] = fast_exp(pp[k] - MAX);
    SUM += pp[k];
  }
  for (int k = 0; k < K; ++k)
    p_I[c][k] = pp[k] / SUM;
}

void ClusterEM::Sample_S(int g) {
  double l0 = buff_l0[g] + 3;
  double l1 = 0;
  for (int k = 0; k < K; ++k) {
    double n0 = 0, n1 = 0;
    for (int c = 0; c < C; ++c) {
      n0 += p_I[c][k] * (1 - Z[c][g]);
      n1 += p_I[c][k] * Z[c][g];
    }
    l1 += n0 * log_theta0[g][k] + n1 * log_theta1[g][k];
  }
  S[g] = l1 > l0;
//  double MAX = std::max(l0, l1);
//  l0 = fast_exp(l0 - MAX);
//  l1 = fast_exp(l1 - MAX);
//  S[g] = R::runif(0, 1) < (l1 / (l0 + l1));
}

void ClusterEM::Update_pi() {
  for (int k = 0; k < K; ++k) {
    double sum_p = 1;
    for (int c = 0; c < C; ++c)
      sum_p += p_I[c][k];
    pi[k] = sum_p / (C + K);
    log_pi[k] = log(pi[k]);
  }
}

void ClusterEM::Update_theta(int g) {
  for (int k = 0; k < K; ++k) {
    double n0 = a_t, n1 = a_t;
    for (int c = 0; c < C; ++c) {
      n0 += p_I[c][k] * (1 - Z[c][g]);
      n1 += p_I[c][k] * (Z[c][g]);
    }
    theta[g][k] = n1 / (n0 + n1);
    log_theta1[g][k] = log(theta[g][k]);
    log_theta0[g][k] = log(1 - theta[g][k]);
  }
}

void ClusterEM::SaveStates() {
  // calculate the observed likelihood given current S
  ll_trace[m] = log_f_obs_S();
}

void ClusterEM::ShowStates() {
  int pct = (m + 1) * 100 / n_steps;
  if (((m + 1) * 10) % n_steps != 0 || pct % 20 != 0) return;

  if (!parallel_print) {
    int curr_nS = 0;
    for (int g = 0; g < G; ++g) curr_nS += S[g];
    cout << "    Finished " << numutils::EnoughSpaces(3, pct) << pct
         << "%, log-likelihood = " << static_cast<int>(ll_trace[m])
         << ", # selected genes = " << curr_nS << ", N(I) = {";
    for (int k = 0; k < K; ++k) {
      double n_k = 0;
      for (int c = 0; c < C; ++c) n_k += p_I[c][k];
      cout << (k + 1) << ": " << round(n_k);
      if (k != K - 1) cout << ", ";
    }
    cout << "}" << endl;
  }
}

void ClusterEM::GetMLE() {
  ll_traces.push_back(NumericVector(ll_trace.begin(), ll_trace.end()));
  if (ll_trace.back() > opt_ll) {
    opt_ll = ll_trace.back();
    opt_S = S;
    opt_p_I = p_I;
    opt_pi = pi;
    opt_theta = theta;
  }
}

void ClusterEM::PostProcess() {
  CalcBIC();
  Cpp2Robjects();
  Reorder_I();

  if (parallel_print) cout << "[K = " << K << "] ";
  cout << "Maximum log-likelihood = " << opt_ll << ", BIC = " << BIC << endl;
}

double ClusterEM::log_f_obs_S() {
  double ll = 0;
  for (int c = 0; c < C; ++c) {
    for (int g = 0; g < G; ++g)
      if (!S[g])
        ll += Z[c][g] * log_omega1[g] + (1 - Z[c][g]) * log_omega0[g];
    vector<double> ll_k = log_pi;
    for (int k = 0; k < K; ++k) {
      for (int g = 0; g < G; ++g)
        if (S[g]) {
          if (Z[c][g])
            ll_k[k] += log_theta1[g][k];
          else
            ll_k[k] += log_theta0[g][k];
        }
    }
    ll += log_sum(ll_k);
  }
  return ll;
}

inline double ClusterEM::fast_exp(double x) {
  if (x < -15) return 0;
  x = 1.0 + x / 32.0;
  x = x * x; x = x * x; x = x * x; x = x * x; x = x * x;
  return x;
}

inline double ClusterEM::fast_log1p(double x) {
  return 0.9404 * x - 0.2537 * x * x;
}

inline double ClusterEM::log_add(double a, double b) {
  if (a >= b) {
    return a + fast_log1p(fast_exp(b - a)); //log(1 + fast_exp(b - a));
  } else {
    return b + fast_log1p(fast_exp(a - b));
  }
}

double ClusterEM::log_sum(const vector<double> & v) {
  double max_val = -10000000000, max_id = -1;

  for (unsigned i = 0; i < v.size(); ++i) {
    if (v[i] > max_val) {
      max_val = v[i];
      max_id  = i;
    }
  }

  double sum = 1;
  for (unsigned i = 0; i < v.size(); i++) {
    if (i == max_id)
      continue;
    sum += fast_exp(v[i] - max_val);
  }

  return max_val + log(sum);
}

void ClusterEM::Cpp2Robjects() {
  // S
  R_S = LogicalVector(opt_S.begin(), opt_S.end());

  // p_I
  R_p_I = NumericMatrix(C, K);
  for (int c = 0; c < C; ++c)
    for (int k = 0; k < K; ++k)
      R_p_I(c, k) = opt_p_I[c][k];

  // I_hat
  R_I_hat = IntegerVector(C);
  for (int c = 0; c < C; ++c) {
    double MAX = 0;
    for (int k = 0; k < K; ++k) {
      if (opt_p_I[c][k] > MAX) {
        MAX = opt_p_I[c][k];
        R_I_hat(c) = k;
      }
    }
  }

  // pi and theta
  R_pi = NumericVector(opt_pi.begin(), opt_pi.end());
  R_theta = NumericMatrix(G, K);
  for (int g = 0; g < G; ++g)
    for (int k = 0; k < K; ++k)
      R_theta(g, k) = opt_theta[g][k];
}

void ClusterEM::CalcBIC() {
  BIC = 0;
}

void ClusterEM::Reorder_I() {
  const double FRAC_THRESH = 0.05;
  const int MIN_N = max(3, static_cast<int>(std::ceil(FRAC_THRESH * C)));
  vector<float> eff_n(K, 0.1);
  vector<float> avg_c(K);

  for (int k = 0; k < K; ++k) {
    for (int c = 0; c < C; ++c) {
      eff_n[k] += R_p_I(c, k);
      avg_c[k] += R_p_I(c, k) * c;
    }
    avg_c[k] /= eff_n[k];
    if (eff_n[k] < MIN_N) avg_c[k] += 1e6;
  }

  vector<int> reorder(K);
  for (int k1 = 0; k1 < K; ++k1) reorder[k1] = k1;
  for (int k1 = 0; k1 < K; ++k1) {
    int rk1 = reorder[k1];
    for (int k2 = k1 + 1; k2 < K; ++k2) {
      int rk2 = reorder[k2];
      double avgc_1 = avg_c[rk1];
      double avgc_2 = avg_c[rk2];
      if (avgc_1 > avgc_2) { // swap rk1 and rk2
        reorder[k2] = rk1;
        rk1 = rk2;
      }
    }
    reorder[k1] = rk1;
  }

  IntegerVector I_hat_copy(clone(R_I_hat));
  NumericMatrix p_I_copy(clone(R_p_I));
  for (int c = 0; c < C; ++c)
    for (int k = 0; k < K; ++k) {
      if (I_hat_copy(c) == reorder[k]) R_I_hat(c) = k;
      R_p_I(c, k) = p_I_copy(c, reorder[k]);
    }
}
