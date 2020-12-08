#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;

#include "cluster_unc.h"
#include "utils.h"


Config R2CppConfig(const List& r_cfg) {
  Config cfg;
  IntegerMatrix Z = as<IntegerMatrix>(r_cfg["Z"]);
  cfg.C = Z.nrow();
  cfg.G = Z.ncol();
  cfg.Z = matrix<bool>(cfg.C, cfg.G);
  for (int c = 0; c < cfg.C; ++c)
    for (int g = 0; g < cfg.G; ++g)
      cfg.Z[c][g] = Z(c, g) > 0;
  cfg.cell_names = as<vector<string> >(rownames(Z));
  cfg.gene_names = as<vector<string> >(colnames(Z));
  cfg.K = as<int>(r_cfg["K"]);
  cfg.n_chain = as<int>(r_cfg["n_chain"]);
  cfg.n_burn = as<int>(r_cfg["n_burn"]);
  cfg.n_mcmc = as<int>(r_cfg["n_mcmc"]);
  cfg.a_w = as<double>(r_cfg["a_w"]);
  cfg.a_t = as<double>(r_cfg["a_t"]);
  cfg.pi_S = as<double>(r_cfg["pi_S"]);
  cfg.verbose = r_cfg.containsElementNamed("verbose") && as<bool>(r_cfg["verbose"]);
  return cfg;
}

void showParameters(const Config& cfg) {
  cout << "# Cells = " << cfg.C << ", # Genes = " << cfg.G << endl;
  cout << "# MCMC chain = " << cfg.n_chain << ", burn-in = " << cfg.n_burn << ", mcmc-steps = " << cfg.n_mcmc << endl;
  cout << "Parameters: K = " << cfg.K << ", a_w = " << cfg.a_w << ", a_t = " << cfg.a_t << ", pi_S = " << cfg.pi_S << endl;
}

// [[Rcpp::export]]
List unc_cluster_K(const List& r_cfg) {
  Config cfg = R2CppConfig(r_cfg);
  // showParameters(cfg);

  ClusterGPU obj(cfg);

  cout << "Start Gibbs sampling......" << endl;
  for (int i = 0; i < cfg.n_chain; ++i) {
    // initialize the MCMC sampling
    obj.InitMCMC();

    // Gibbs sampling
    obj.ClusterCells(i);

    // get posterior mode
    obj.GetPosteriorMode();
  }

  // posterior processing
  obj.PostProcess();

  return List::create(Named("p_S") = obj.p_S,
                      Named("p_I") = obj.p_I,
                      Named("I_hat") = obj.I_hat,
                      Named("ll_traces") = obj.ll_traces,
                      Named("log_p_Z_I") = obj.log_p_Z_I,
                      Named("log_p_Z_K") = obj.log_p_Z_K);
}
