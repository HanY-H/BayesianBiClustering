#include "cluster.h"

#include <cmath>

#include <vector>
#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"


void Cluster::CalcCaches() {
  buff_log_l0 = vector<double>(G);
  buff_log_l1 = vector<double>(G);
  for (int g = 0; g < G; ++g) {
    int n1 = 0;
    for (int c = 0; c < C; ++c) n1 += Z[c][g];
    buff_log_l0[g] = log_1_pi_S + logBinomML(a_w, a_w, n1, C - n1);
  }
  buff_log_binom_ML_theta = matrix<double>(C + 1, C + 1);
  for (int n1 = 0; n1 <= C; ++n1)
    for (int n0 = 0; n0 <= C; ++n0)
      buff_log_binom_ML_theta[n1][n0] = logBinomML(a_t, a_t, n1, n0);

  buff_log_binom_ML_theta_gain = matrix<vector<double> >(C + 1, C + 1, vector<double>(2));
  buff_log_binom_ML_theta_lose = matrix<vector<double> >(C + 1, C + 1, vector<double>(2));
  for (int n1 = 0; n1 <= C; ++n1) {
    for (int n0 = 0; n0 <= C; ++n0) {
      if (n1 != C)
        buff_log_binom_ML_theta_gain[n1][n0][1] = buff_log_binom_ML_theta[n1 + 1][n0] - buff_log_binom_ML_theta[n1][n0];
      if (n0 != C)
        buff_log_binom_ML_theta_gain[n1][n0][0] = buff_log_binom_ML_theta[n1][n0 + 1] - buff_log_binom_ML_theta[n1][n0];
      if (n1 != 0)
        buff_log_binom_ML_theta_lose[n1][n0][1] = buff_log_binom_ML_theta[n1 - 1][n0] - buff_log_binom_ML_theta[n1][n0];
      if (n0 != 0)
        buff_log_binom_ML_theta_lose[n1][n0][0] = buff_log_binom_ML_theta[n1][n0 - 1] - buff_log_binom_ML_theta[n1][n0];
    }
  }
}

void Cluster::Compute_Nk() {
  N_k = matrix<vector<int> >(K, G, vector<int>(2));
  for (int g = 0; g < G; ++g)
    for (int c = 0; c < C; ++c)
      N_k[I[c]][g][Z[c][g]] ++;
}

void Cluster::Update_Nk(int c, int old_k, int new_k) {
  if (old_k == new_k) return;
  auto& N_k_old = N_k[old_k];
  auto& N_k_new = N_k[new_k];
  auto& Z_c = Z[c];
  for (int g = 0; g < G; ++g) {
    bool Z_c_g = Z_c[g];
    N_k_old[g][Z_c_g] --;
    N_k_new[g][Z_c_g] ++;
  }
}

void Cluster::InitMCMC() {
  // init MCMC states
  I = as<vector<int> >(static_cast<IntegerVector>(sample(K, C, true) - 1));
  if (I_start.size() > 0) {
    int min_K = 10000, max_K = -1;
    for (int c = 0; c < C; ++c) {
      min_K = min(I_start[c], min_K);
      max_K = max(I_start[c], max_K);
    }
    if (min_K >= 0 && max_K < K) I = I_start;
  }

  Compute_Nk();
  Rebuild_Buff_Log_L1();

  // init MCMC traces
  pm_ll = -10000000000;
  t_I  = vector<vector<int> >(n_m);
  t_ll = vector<double>(n_b + n_m, 0);

  // tmp variables
  pp  = vector<double>(K, -10000000000);
  qq  = vector<double>(K, -10000000000);
  update_I_freq = vector<double>(C, 1.0);
}

void Cluster::ClusterCells(int chain_id) {
  if (!parallel_print)
    cout << "  Chain " << (chain_id + 1) << ":" << endl;

  for (m = 0; m < (n_b + n_m); ++m) {
    for (int c = 0; c < C; ++c) {
      if (verbose) cout << "ClusterCells::{m = " << m << " / " << (n_b + n_m) << ", c = " << c << ", update_I_freq=" << update_I_freq[c] << "}" << endl;
      Update_Ic(c);
    }

    ShowStates();

    SaveStates();
  }
  if (parallel_print)
    cout << "[K = " << K << "] Finished Chain " << chain_id << ", maximum posterior probability = " << pm_ll << ", # selected genes = " << curr_nS << ((pm_ll > ppm_ll) ? string("  (new max)") : string()) << endl;
  else
    cout << "    Maximum posterior probability = " << pm_ll << ((pm_ll > ppm_ll) ? string(" (new max)") : string()) << endl;
}

void Cluster::Rebuild_Buff_Log_L1() {
  for (int g = 0; g < G; ++g) {
    double l1 = log_pi_S;
    for (int k = 0; k < K; ++k) {
      auto& N_k_g = N_k[k][g];
      l1 += buff_log_binom_ML_theta[N_k_g[1]][N_k_g[0]];
    }
    buff_log_l1[g] = l1;
  }
}

// updating sub-population assignments under finite mixture model
void Cluster::Update_Ic(int c) {
//  if (R::runif(0, 1) > update_I_freq[c])
//    return;

  int old_k = I[c];

  // initialization of classifying probabilities
  for (int k = 0; k < K; ++k)
    pp[k] = log_f_Z_Ic(c, k);

  double MAX = -10000000000;

  // compute the annealing temperature
  double T = 1;
  if (m < n_b * 0.6)
    T = T0 + m / (n_b * 0.6) * (1 - T0);

  for (int k = 0; k < K; ++k)
    if (pp[k] > MAX)
      MAX = pp[k];

  double sum = 0;
  for (int k = 0; k < K; ++k) {
    qq[k] = fast_exp(pp[k] - MAX);
    qq[k] = pow(qq[k], 1 / T);	// tempering
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

  if (!(kk >= 0 && kk < K)) {
    cout << "Update_Ic::pp = ";
    ctnutils::DispVector(pp);

    for (int k = 0; k < K; ++k)
      if (pp[k] > MAX)
        MAX = pp[k];
    cout << "Update_Ic::MAX = " << MAX << endl;

    double sum = 0;
    for (int k = 0; k < K; ++k) {
      qq[k] = fast_exp(pp[k] - MAX);
      qq[k] = pow(qq[k], 1 / T);	// tempering
      sum = sum + qq[k];
    }

    cout << "Update_Ic::qq1 = ";
    ctnutils::DispVector(qq);
    for (int k = 0; k < K; ++k)
      qq[k] = qq[k] / sum;

    cout << "Update_Ic::qq2 = ";
    ctnutils::DispVector(qq);

    cout << endl;
  }

  if (kk < K) {
    if (I[c] != kk)
    for (int g = 0; g < G; ++g) {
      bool Z_c_g = Z[c][g];
//      buff_log_l1[g] = buff_log_l1[g] + (buff_log_binom_ML_theta[N_k[old_k][g][1] - Z_c_g][N_k[old_k][g][0] - 1 + Z_c_g] - buff_log_binom_ML_theta[N_k[old_k][g][1]][N_k[old_k][g][0]]) + (buff_log_binom_ML_theta[N_k[kk][g][1] + Z_c_g][N_k[kk][g][0] + 1 - Z_c_g] - buff_log_binom_ML_theta[N_k[kk][g][1]][N_k[kk][g][0]]);
      buff_log_l1[g] += buff_log_binom_ML_theta_lose[N_k[old_k][g][1]][N_k[old_k][g][0]][Z_c_g] + buff_log_binom_ML_theta_gain[N_k[kk][g][1]][N_k[kk][g][0]][Z_c_g];
    }
    I[c] = kk;
    Update_Nk(c, old_k, kk);
  }

  if (old_k != I[c])
    update_I_freq[c] = min(update_I_freq[c] + 0.2, 1.0);
  else
    update_I_freq[c] = max(update_I_freq[c] - 0.1, MIN_UPDATE_FREQ);
}

void Cluster::SaveStates() {
  // calculate marginal likelihood
  double ll = 0;
  for (int g = 0; g < G; ++g) {
    double l0 = buff_log_l0[g];
    double l1 = buff_log_l1[g];
    ll += log_add(l0, l1);
  }
  t_ll[m] = ll;

  // save states
  if (m >= n_b) {
    int mm = m - n_b;
    t_I[mm] = I;
    if (ll > pm_ll) {
      pm_ll = ll;
      pm_I  = I;
    }
  }
}

void Cluster::ShowStates() {
  int n_total = n_b + n_m;
  if (((m + 1) * 10) >= n_total && ((m + 1) * 10) % n_total == 0) {
    int pct = (m + 1) * 100 / n_total;
    if (pct % 20 == 0) {
      curr_nS = 0;
      for (int g = 0; g < G; ++g) curr_nS += (log_f_Sg_I(g) > 0);
      curr_ll = static_cast<int>(log_f_Z_I());
      if (!parallel_print) {
        cout << "    Finished " << numutils::EnoughSpaces(3, pct) << pct << "%, log-likelihood = " << curr_ll << ", # selected genes = " << curr_nS << ", N(I) = {";
        for (int k = 0; k < K; ++k) {
          int n_k = 0;
          for (int c = 0; c < C; ++c) n_k += I[c] == k;
          cout << (k + 1) << ": " << n_k;
          if (k != K - 1) cout << ", ";
        }
        cout << "}" << endl;
      }
    }
  }
}

void Cluster::GetPosteriorMode() {
  ll_traces.push_back(NumericVector(t_ll.begin(), t_ll.end()));
  if (pm_ll > ppm_ll) {
    ppm_I = pm_I;
    ppm_ll = pm_ll;
    ppm_t_I  = t_I;
    ppm_t_ll = t_ll;
  }
}

void Cluster::PostProcess() {
  I = ppm_I;
  t_I = ppm_t_I;

  CalcPostProbs();
  CalcMargLikMethod1();

  if (parallel_print) {
    cout << "[K = " << K << "] Maximum posterior probability from multiple chains = " << ppm_ll << endl;
    // cout << "[K = " << K << "] Cache hit rate = " << (static_cast<double>(cache_hit) / (cache_hit + cache_miss)) << endl;
    cout << "[K = " << K << "] Marginal likelihood log{P(Z|M_K)} = " << log_p_Z_K << endl;
  } else {
    cout << "Maximum posterior probability from multiple chains = " << ppm_ll << endl;
    cout << "Cache hit rate = " << (static_cast<double>(cache_hit) / (cache_hit + cache_miss)) << endl;
    cout << "Calculating marginal likelihood log{P(Z|M_K)} = " << log_p_Z_K << endl;
  }

  Reorder_I();
}

 double Cluster::logBinomML(double a_p, double b_p, double a_c, double b_c) {
  return R::lbeta(a_p + a_c, b_p + b_c) - R::lbeta(a_p, b_p);
}

double Cluster::log_f_Z_Ic(int c, int new_k) {
  int old_k = I[c];
  I[c] = new_k;

  auto iter = ll_db.find(I);
  if (iter != ll_db.end()) {
    cache_hit++;
    I[c] = old_k;
    return iter->second;
  } else {
    cache_miss++;
  }

  auto& N_k_old_k = N_k[old_k];
  auto& N_k_new_k = N_k[new_k];
  // calculate ll
  double ll = 0;
  for (int g = 0; g < G; ++g) {
    double l0 = buff_log_l0[g];
    double l1 = buff_log_l1[g];
    if (new_k != old_k) {
      bool Z_c_g = Z[c][g];
      auto& n_old = N_k_old_k[g];
      auto& n_new = N_k_new_k[g];
      l1 += buff_log_binom_ML_theta_lose[n_old[1]][n_old[0]][Z_c_g] + buff_log_binom_ML_theta_gain[n_new[1]][n_new[0]][Z_c_g];
    }
    ll += (l1 - l0 > 10) ? l1 : log_add(l0, l1);
  }

  if (ll_db.size() > max_cache_size) ll_db.clear();
  ll_db[I] = ll;

  I[c] = old_k;
  return ll;
}

double Cluster::log_f_Z_I() {
  return log_f_Z_I(I);
}

double Cluster::log_f_Z_I(const vector<int> & Ip) {
  double ll = 0;
  for (int g = 0; g < G; ++g) {
    double l0 = buff_log_l0[g];
    double l1 = log_pi_S;
    for (int k = 0; k < K; ++k) {
      auto& N_k_g = N_k[k][g];
      l1 += buff_log_binom_ML_theta[N_k_g[1]][N_k_g[0]];
      // cout << "k = " << k << ", n0 = " << N_k_g[0] << ", n1 = " << N_k_g[1] << ", ll = " << buff_log_binom_ML_theta[N_k_g[1]][N_k_g[0]] << endl;
    }
    // cout << "l0 = " << l0 << ", l1 = " << l1 << endl;
    ll += log_add(l0, l1);
  }
  return ll;
}

double Cluster::log_f_Sg_I(int g) {
  double l0 = buff_log_l0[g];
  double l1 = log_pi_S;
  for (int k = 0; k < K; ++k) l1 += buff_log_binom_ML_theta[N_k[k][g][1]][N_k[k][g][0]];
  return l1 - l0;
}

inline double Cluster::fast_exp(double x) {
  if (x < -15) return 0;
  x = 1.0 + x / 32.0;
  x = x * x; x = x * x; x = x * x; x = x * x; x = x * x;
  return x;
}

inline double Cluster::fast_log1p(double x) {
  return 0.9404 * x - 0.2537 * x * x;
}

double Cluster::log_beta_pdf(double x, double a, double b) {
  return log(R::dbeta(x, a, b, 0));
}

inline double Cluster::log_add(double a, double b) {
//  if (a >= b) {
//    return a + fast_log1p(fast_exp(b - a)); //log(1 + fast_exp(b - a));
//  } else {
//    return b + fast_log1p(fast_exp(a - b));
//  }
  if (a >= b) {
    return a + log(1 + exp(b - a)); //log(1 + fast_exp(b - a));
  } else {
    return b + log(1 + exp(a - b));
  }
}

double Cluster::log_sum(const vector<double> & v) {
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

void Cluster::CalcPostProbs() {
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

vector<double> Cluster::calc_n_g(const vector<int>& II, int g) {
  vector<double> n_g(2);
  for (int c = 0; c < C; ++c) n_g[Z[c][g]]++;
  return n_g;
}

matrix<double> Cluster::calc_n_k(const vector<int>& II, int g) {
  matrix<double> n_k(K, 2);
  for (int c = 0; c < C; ++c) n_k[II[c]][Z[c][g]] ++;
  return n_k;
}

void Cluster::EstimateDeltaK() {
  vector<matrix<double> > t_logS(n_m, matrix<double>(G, 2));
  // run Gibbs sampling
  for (int m = 0; m < n_m; ++m) {
    for (int c = 0; c < C; ++c)
      Update_Ic(c);
    t_I[m] = I;
    for (int g = 0; g < G; ++g) {
      pair<double, double> log_pq = calc_log_pq_from_LLR(log_f_Sg_I(g));
      t_logS[m][g][0] = log_pq.first;
      t_logS[m][g][1] = log_pq.second;
    }
  }

  // calculate omega*
  vector<double> omega_s(G);
  for (int g = 0; g < G; ++g) {
    for (int m = 0; m < n_m; ++m) {
      vector<double> n_g = calc_n_g(t_I[m], g);
      omega_s[g] += (a_w + n_g[1]) / (2 * a_w + n_g[0] + n_g[1]) / double(n_m);
    }
  }

  // calculate theta*
  matrix<double> theta_s(G, K);
  for (int g = 0; g < G; ++g) {
    for (int m = 0; m < n_m; ++m) {
      matrix<double> n_k = calc_n_k(t_I[m], g);
      for (int k = 0; k < K; ++k)
        theta_s[g][k] += (a_t + n_k[k][1]) / (2 * a_t + n_k[k][0] + n_k[k][1]) / double(n_m);
    }
  }

  // 1st term as product of observed likelihoods
  double term_1 = 0;
  for (int c = 0; c < C; ++c) {
    for (int g = 0; g < G; ++g) {
      double lS0 = log_1_pi_S + Z[c][g] * log(omega_s[g]) + (1 - Z[c][g]) * log(1 - omega_s[g]);
      vector<double> log_f_k(K, 0);
      for (int k = 0; k < K; ++k)
        log_f_k[k] = Z[c][g] * log(theta_s[g][k]) + (1 - Z[c][g]) * log(1 - theta_s[g][k]);
      double lS1 = log_pi_S - log(K) + log_sum(log_f_k);
      term_1 += log_add(lS0, lS1);
    }
  }

  // 2nd term as product of beta priors
  double term_2 = 0;
  for (int g = 0; g < G; ++g)
    term_2 += log_beta_pdf(omega_s[g], a_w, a_w);
  for (int g = 0; g < G; ++g)
    for (int k = 0; k < K; ++k)
      term_2 += log_beta_pdf(theta_s[g][k], a_t, a_t);

  // approximate the 3rd term
  // without permutation of I - note that it is wrong!
  double term_3 = 0;
  vector<double> log_ps(n_m, 0);
  for (int m = 0; m < n_m; ++m) {
    for (int g = 0; g < G; ++g) {
      double lS0 = t_logS[m][g][0];
      vector<double> n_g = calc_n_g(t_I[m], g);
      lS0 += log_beta_pdf(omega_s[g], a_w + n_g[1], a_w + n_g[0]);
      for (int k = 0; k < K; ++k)
        lS0 += log_beta_pdf(theta_s[g][k], a_t, a_t);

      double lS1 = t_logS[m][g][1];
      lS1 += log_beta_pdf(omega_s[g], a_w, a_w);
      matrix<double> n_k = calc_n_k(t_I[m], g);
      for (int k = 0; k < K; ++k)
        lS1 += log_beta_pdf(theta_s[g][k], a_t + n_k[k][1], a_t + n_k[k][0]);

      log_ps[m] += log_add(lS0, lS1);
    }
  }
  term_3 = log_sum(log_ps) - log(n_m);
  log_p_Z_K = term_1 + term_2 - term_3;
}

void Cluster::CalcMargLikMethod1() {
  // 1. logP(Z|M_K, I)
  I = ppm_I;
  Compute_Nk();
  double term1 = log_f_Z_I(I);
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

void Cluster::Reorder_I() {
  const double FRAC_THRESH = 0.05;
  const int MIN_N = max(3, static_cast<int>(std::ceil(FRAC_THRESH * C)));
  vector<float> eff_n(K, 0.1);
  vector<float> avg_c(K);

  for (int k = 0; k < K; ++k) {
    for (int c = 0; c < C; ++c) {
      eff_n[k] += p_I(c, k);
      avg_c[k] += p_I(c, k) * c;
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

  IntegerVector I_hat_copy(clone(I_hat));
  NumericMatrix p_I_copy(clone(p_I));
  for (int c = 0; c < C; ++c)
    for (int k = 0; k < K; ++k) {
      if (I_hat_copy(c) == reorder[k]) I_hat(c) = k;
      p_I(c, k) = p_I_copy(c, reorder[k]);
    }
}

double Cluster::calc_p_from_LLR(double LLR) {
  double mp = max(LLR, 0.0);
  double p0 = fast_exp(-mp);
  double p1 = fast_exp(LLR - mp);
  p1 = p1 / (p0 + p1);
  return p1;
}

pair<double, double> Cluster::calc_log_pq_from_LLR(double LLR) {
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
