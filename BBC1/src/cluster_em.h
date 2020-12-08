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

  int n_rounds;
  int n_steps;

  double a_w;
  double a_t;

  bool verbose;
  bool parallel_print;
};

class ClusterEM {
 public:
  ClusterEM(const Config& cfg) :
      Z{cfg.Z},
      cell_names{cfg.cell_names},
      gene_names{cfg.gene_names},
      C{cfg.C},
      G{cfg.G},
      K{cfg.K},
      n_steps{cfg.n_steps},
      a_w{cfg.a_w},
      a_t{cfg.a_t},
      verbose{cfg.verbose},
      parallel_print{cfg.parallel_print} {
    CalcCaches();
  }

  ~ClusterEM() {}

  // initialize the Gibbs sampler
  void InitEM();

  // Gibbs sampling
  void RunEM(int round_id);

  // find posterior mode
  void GetMLE();

  // posterior processing
  void PostProcess();

  // outputs to Rs
  LogicalVector R_S;
  NumericMatrix R_p_I;
  IntegerVector R_I_hat;
  NumericVector R_pi;
  NumericMatrix R_theta;
  List ll_traces;
  double BIC;

 private:
  // calculate caches for likelihood functions
  void CalcCaches();

  // update p_I
  void Update_p_I(int c);

  // sample S
  void Sample_S(int g);

  // update pi
  void Update_pi();

  // update theta
  void Update_theta(int g);

  // monitor and save EM chain
  void SaveStates();

  // show the states of current EM step
  void ShowStates();

  // calculate posterior probabilities
  void Cpp2Robjects();

  // calculate BIC
  void CalcBIC();

  // reorder clusters so that small i have small I
  void Reorder_I();

  /**
  *** calculations of probability distributions
  **/
  double log_f_obs_S();

  /**
  *** utility functions
  **/
  // fast exp(x); very accurate for x < 0
  inline double fast_exp(double x);

  // fast log(1 + x); very accurate for 0 < x < 1
  inline double fast_log1p(double x);

  inline double log_add(double a, double b);

  double log_sum(const vector<double> & v);

  double calc_p_from_LLR(double LLR);

  /** Member variables */
  // input data
  matrix<bool> Z;
  vector<string> cell_names;
  vector<string> gene_names;

  // num of cells and genes
  int C;
  int G;

  // running parameters
  int m;
  int K;
  int n_steps;

  // model hyperparameters
  double a_t;  // alpha theta
  double a_w;  // alpha omega

  // caches
  vector<double> omega;
  vector<double> log_omega1;
  vector<double> log_omega0;
  vector<double> buff_l0;

  // EM states
  vector<bool> S;
  matrix<double> p_I;
  vector<double> pi;
  matrix<double> theta;
  vector<double> ll_trace;
  // log caches
  vector<double> log_pi;
  matrix<double> log_theta1;
  matrix<double> log_theta0;

  // MLE over all EM chains
  double opt_ll = -10000000000;
  vector<bool> opt_S;
  matrix<double> opt_p_I;
  vector<double> opt_pi;
  matrix<double> opt_theta;

  // tmp variables
  vector<double> pp;
  vector<double> qq;

  bool verbose;
  bool parallel_print;
  RNGScope rngScope;
};


#endif // CLUSTER_Z_H_INCLUDED


