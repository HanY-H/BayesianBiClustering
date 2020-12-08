#ifndef CLUSTER_Z_H_INCLUDED
#define CLUSTER_Z_H_INCLUDED

#include <map>
#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"


struct Config {
  int C;
  int G;

  matrix<bool> Z;
  vector<string> cell_names;
  vector<string> gene_names;

  vector<int> I_start;

  int K;

  int n_chain;
  int n_burn;
  int n_mcmc;
  double T0;

  double a_w;
  double a_t;

  double pi_S;

  int max_cache_size = 250000;

  bool verbose;
  bool parallel_print;
};

class Cluster {
 public:
  Cluster(const Config& cfg) :
      Z{cfg.Z},
      cell_names{cfg.cell_names},
      gene_names{cfg.gene_names},
      C{cfg.C},
      G{cfg.G},
      K{cfg.K},
      n_b{cfg.n_burn},
      n_m{cfg.n_mcmc},
      T0{cfg.T0},
      a_w{cfg.a_w},
      a_t{cfg.a_t},
      pi_S{cfg.pi_S},
      I_start{cfg.I_start},
      max_cache_size{cfg.max_cache_size},
      verbose{cfg.verbose},
      parallel_print{cfg.parallel_print} {
    // log caches
    log_pi_S = log(pi_S);
    log_1_pi_S = log(1 - pi_S);
    if (pi_S < 0.000001) {
      log_pi_S = -1e6;
      log_1_pi_S = 0;
    }
    if (pi_S > 0.999999) {
      log_pi_S = 0;
      log_1_pi_S = -1e6;
    }

    ll_db = map<vector<int>, double>();
    ppm_ll = -10000000000;
    CalcCaches();
  }

  ~Cluster() {}

  // initialize the Gibbs sampler
  void InitMCMC();

  // Gibbs sampling
  void ClusterCells(int chain_id);

  // find posterior mode
  void GetPosteriorMode();

  // posterior processing
  void PostProcess();

  // output: posterior probabilities
  NumericVector p_S;
  NumericMatrix p_I;
  IntegerVector I_hat;
  List ll_traces;

  // output: marginal likelihood log{P(Z|I)} and log{P(Z|M_K)}
  double log_p_Z_I;
  double log_p_Z_K;

 private:
  // calculate caches for likelihood functions
  void CalcCaches();

  // compute N_k
  void Compute_Nk();

  // update N_k
  void Update_Nk(int c, int old_k, int new_k);

  // rebuild buff_log_l1
  void Rebuild_Buff_Log_L1();

  // update Ic
  void Update_Ic(int c);

  // monitor and save MCMC chain
  void SaveStates();

  // show the states of current MCMC step
  void ShowStates();

  // calculate the posterior probabilities
  void CalcPostProbs();

  // estimate logP(Z | M_K)
  void EstimateDeltaK();
  void CalcMargLikMethod1();

  // reorder clusters so that small i have small I
  void Reorder_I();

  /**
  *** calculations of probability distributions
  **/
  // log Binomial-Beta marginal likelihood
  double logBinomML(double a_p, double b_p, double a_c, double b_c);

  double log_f_Z_Ic(int c, int k);

  double log_f_Z_I();

  double log_f_Z_I(const vector<int>& Ip);

  double log_f_Sg_I(int g);

  // calculate n_g and n_g_k based on I
  vector<double> calc_n_g(const vector<int>& II, int g);
  matrix<double> calc_n_k(const vector<int>& II, int g);

  /**
  *** utility functions
  **/
  // fast exp(x); very accurate for x < 0
  inline double fast_exp(double x);

  // fast log(1 + x); very accurate for 0 < x < 1
  inline double fast_log1p(double x);

  inline double log_add(double a, double b);

  double log_beta_pdf(double x, double a, double b);

  double log_sum(const vector<double> & v);

  double calc_p_from_LLR(double LLR);

  pair<double, double> calc_log_pq_from_LLR(double LLR);

  /** Member variables */
  // input data
  matrix<bool> Z;
  vector<string> cell_names;
  vector<string> gene_names;

  // num of cells and genes
  int C;
  int G;

  // running parameters
  int n_b;
  int n_m;
  double T0;

  // model hyperparameters
  double a_t;  // alpha theta
  double a_w;  // alpha omega
  double pi_S;

  // computed log caches
  double log_pi_S;
  double log_1_pi_S;

  // log-likelihood function caches
  vector<double> buff_log_l0;
  vector<double> buff_log_l1;
  matrix<double> buff_log_binom_ML_theta;
  matrix<vector<double> > buff_log_binom_ML_theta_gain;
  matrix<vector<double> > buff_log_binom_ML_theta_lose;
  map<vector<int>, double> ll_db;
  int cache_hit = 0;
  int cache_miss = 0;

  // MCMC states
  int m;
  int K;
  vector<int> I;

  // MCMC traces
  vector<double> t_ll;
  vector<vector<int> > t_I;

  // posterior mode for current MCMC chain
  double pm_ll;
  vector<int> pm_I;

  // posterior mode over all MCMC chains
  double ppm_ll = -10000000000;
  vector<int> ppm_I;

  // MCMC trace corresponding to the mode
  vector<double> ppm_t_ll;
  vector<vector<int> > ppm_t_I;

  // tmp variables
  vector<double> pp;
  vector<double> qq;
  matrix<vector<int> > N_k;
  vector<double> update_I_freq;
  int curr_ll;
  int curr_nS;

  // constants
  const double MIN_UPDATE_FREQ = 0.1;

  vector<int> I_start;
  int max_cache_size;
  bool verbose;
  bool parallel_print;
  RNGScope rngScope;
};


#endif // CLUSTER_Z_H_INCLUDED


