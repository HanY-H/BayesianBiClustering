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

  int K;

  int n_chain;
  int n_burn;
  int n_mcmc;

  double a_w;
  double a_t;

  double pi_S;

  bool verbose;
};

class ClusterGPU {
 public:
  ClusterGPU(const Config& cfg) :
      Z1{cfg.Z},
      cell_names{cfg.cell_names},
      gene_names{cfg.gene_names},
      C{cfg.C},
      G{cfg.G},
      K{cfg.K},
      n_b{cfg.n_burn},
      n_m{cfg.n_mcmc},
      a_w{cfg.a_w},
      a_t{cfg.a_t},
      pi_S{cfg.pi_S},
      verbose{cfg.verbose} {
    // Z0
    Z0 = matrix<bool>(C, G);
    for (int c = 0; c < C; ++c)
      for (int g = 0; g < G; ++g)
        Z0[c][g] = 1 - Z1[c][g];

    ppm_ll = -10000000000;
    CalcCaches();
  }

  ~ClusterGPU() {}

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

  // update I_c
  void Update_I_c(int c);

  // update S_g
  void Update_S_g(int g);

  // update theta_g
  void Update_theta_g(int g);

  // monitor and save MCMC chain
  void SaveStates();

  // calculate the posterior probabilities
  void CalcPostProbs();

  // estimate logP(Z | M_K)
  void CalcMargLikMethod1();

  /**
  *** calculations of probability distributions
  **/
  // log Binomial-Beta marginal likelihood
  double logBinomML(double a_p, double b_p, double a_c, double b_c);

  double log_f_Z_Ic(int c, int k);

  double log_f_Z_I();

  double log_f_Sg_I(int g);

  /**
  *** utility functions
  **/
  // fast exp(x); very accurate for x < 0
  inline double fast_exp(double x);

  // fast log(1 + x); very accurate for 0 < x < 1
  inline double fast_log1p(double x);

  inline double log_add(double a, double b);

  double calc_p_from_LLR(double LLR);

  pair<double, double> calc_log_pq_from_LLR(double LLR);

  /** Member variables */
  // input data
  matrix<bool> Z0;
  matrix<bool> Z1;
  vector<string> cell_names;
  vector<string> gene_names;

  // num of cells and genes
  int C;
  int G;

  // running parameters
  int n_b;
  int n_m;

  // model hyperparameters
  double a_t;  // alpha theta
  double a_w;  // alpha omega
  double pi_S;

  // caches
  double log_pi_S1;
  double log_pi_S0;
  vector<double> omega;
  vector<double> log_omega1;
  vector<double> log_omega0;
  matrix<double> log_theta1;
  matrix<double> log_theta0;
  vector<double> buff_log_l0;
  matrix<double> buff_log_binom_ML_theta;

  // MCMC states
  int m;
  int K;
  vector<int> I;
  vector<bool> S0;
  vector<bool> S1;
  matrix<double> theta;

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

  bool verbose;
  RNGScope rngScope;
};


#endif // CLUSTER_Z_H_INCLUDED


