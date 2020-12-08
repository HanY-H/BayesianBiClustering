#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;

#include "cluster_em.h"
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
  cfg.n_rounds = as<int>(r_cfg["n_rounds"]);
  cfg.n_steps = as<int>(r_cfg["n_steps"]);
  cfg.a_w = as<double>(r_cfg["a_w"]);
  cfg.a_t = as<double>(r_cfg["a_t"]);
  cfg.verbose = r_cfg.containsElementNamed("verbose") && as<bool>(r_cfg["verbose"]);
  cfg.parallel_print = r_cfg.containsElementNamed("parallel_print") && as<bool>(r_cfg["parallel_print"]);
  return cfg;
}

void showParameters(const Config& cfg) {
  cout << "# Cells = " << cfg.C << ", # Genes = " << cfg.G << endl;
  cout << "# EM rounds = " << cfg.n_rounds << ", steps = " << cfg.n_steps << ", K = " << cfg.K << ", a_w = " << cfg.a_w << ", a_t = " << cfg.a_t << endl;
}

// [[Rcpp::export]]
List cpp_run_EM_K(const List& r_cfg) {
  Config cfg = R2CppConfig(r_cfg);
  // showParameters(cfg);

  ClusterEM obj(cfg);

  if (cfg.parallel_print)
    cout << "[K = " << cfg.K << "] Starting EM......" << endl;
  else
    cout << "Starting EM......" << endl;
  for (int i = 0; i < cfg.n_rounds; ++i) {
    obj.InitEM();

    obj.RunEM(i);

    obj.GetMLE();
  }
  obj.PostProcess();

  return List::create(Named("S") = obj.R_S,
                      Named("p_I") = obj.R_p_I,
                      Named("I_hat") = obj.R_I_hat,
                      Named("pi") = obj.R_pi,
                      Named("theta") = obj.R_theta,
                      Named("ll_traces") = obj.ll_traces,
                      Named("BIC") = obj.BIC);
}
