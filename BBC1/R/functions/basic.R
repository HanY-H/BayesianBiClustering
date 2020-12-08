library(Rcpp)

basic = function(Z, K_from=2, K_to=3, a_w=1, a_t=1, T0=1, pi_S=0.1, n_chain=1, n_burn=100, n_mcmc=100, verbose=F, unc=F, gpu=F, I_start=NULL) {
  set.seed(0)

  cell_names = rownames(Z)
  gene_names = colnames(Z)
  
  res = list()
  cfg = list(Z=Z, K=1, a_w=a_w, a_t=a_t, T0=T0, pi_S=pi_S, n_chain=n_chain, n_burn=n_burn, n_mcmc=n_mcmc, verbose=verbose)
  if (!is.null(I_start)) cfg$I_start = I_start

  cpp_func = NULL
  if (unc) {
    cpp_func = unc_cluster_K;
  } else if (gpu) {
    cpp_func = gpu_cluster_K;
  } else {
    cpp_func = cpp_cluster_K;
  }
  
  cat('K = 1\n')
  res_1 = cpp_func(cfg)
  for(K in K_from:K_to) {
    cat(paste0('\nK = ', K, ':\n'))
    cfg$K = K
    res_K = cpp_func(cfg)
    res_K$K = K
    res_K$cfg = cfg
    res_K$L_K = res_K$log_p_Z_I - res_1$log_p_Z_I
    res_K$D_K = res_K$log_p_Z_K - res_1$log_p_Z_K #log_p_Z_K Marginal Log-Likelihood P(Z|M_K)
    names(res_K$p_S) = gene_names
    rownames(res_K$p_I) = cell_names
    colnames(res_K$p_I) = paste('Cluster', 1:K)
    res[[paste0('K', K)]] = res_K
    cat(paste0('D', K, ' = ', format(res_K$D_K, 1), "\n"))
  }
  
  return(res)
}

basic_EM = function(Z, K_from=2, K_to=3, a_w=1, a_t=1, n_rounds=1, n_steps=100, verbose=F) {
  set.seed(0)
  
  cell_names = rownames(Z)
  gene_names = colnames(Z)

  res = list()
  cfg = list(Z=Z, a_w=a_w, a_t=a_t, n_rounds=n_rounds, n_steps=n_steps, verbose=verbose)

  for(K in K_from:K_to) {
    cat(paste0('\nK = ', K, ':\n'))
    cfg$K = K
    res_K = cpp_run_EM_K(cfg)
    res_K$K = K
    res_K$cfg = cfg
    names(res_K$S) = gene_names
    rownames(res_K$p_I) = cell_names
    colnames(res_K$p_I) = paste('Cluster', 1:K)
    names(res_K$pi) = paste('Cluster', 1:K)
    names(res_K$I_hat) = cell_names
    rownames(res_K$theta) = gene_names
    colnames(res_K$theta) = paste('Cluster', 1:K)
    res[[paste0('K', K)]] = res_K
    cat(paste0('D', K, ' = ', format(res_K$D_K, 1), "\n"))
  }
  
  return(res)
}

select_D_K = function(res) {
  best_K = 0
  best_D_K = -1e9
  for(res_K in res) {
    K = res_K$cfg$K
    if (res_K$D_K > best_D_K) {
      best_K = K
      best_D_K = res_K$D_K
    }
  }
  for(res_K in res) {
    K = res_K$cfg$K
    cat(paste0('D(', K, ') = ', round(res_K$D_K)))
    if (K == best_K) cat(' *')
    cat('\n')
  }
}


select_BIC = function(res, gamma=0.5) {
  best_K = 0
  best_BIC = -1e9
  Ks = numeric()
  BICs = numeric()
  for(res_K in res) {
    n = length(res_K$I_hat)
    K = res_K$cfg$K
    L_K = res_K$L_K
    BIC = L_K - gamma * log(n) * K
    if (BIC > best_BIC) {
      best_K = K
      best_BIC = BIC
    }
    Ks = c(Ks, K)
    BICs = c(BICs, BIC)
  }
  for(i in 1:length(Ks)) {
    cat(paste0('BIC(', Ks[i], ') = ', round(BICs[i])))
    if (Ks[i] == best_K) cat(' *')
    cat('\n')
  }
}

select_K = function(res, gamma=0.5) {
  select_D_K(res)
  cat('\n')
  select_BIC(res, gamma)
}


log_f_Z_I = function(Z, I, pi_S, a_w, a_t) {
  
}

cppFunction('
            double log_f_Z_I(IntegerMatrix Z, IntegerVector I, double pi_S, double a_w, double a_t) {
            double log_pi_S0 = log(1 - pi_S);
            double log_pi_S1 = log(pi_S);
            int C = Z.nrow();
            int G = Z.ncol();
            int max_I = max(I);
            double ll = 0;
            for (int g = 0; g < G; ++g) {
            // l0
            int n1 = 0;
            for (int c = 0; c < C; ++c) n1 += Z(c, g);
            double l0 = log_pi_S0 + R::lbeta(a_w + n1, a_w + C - n1) - R::lbeta(a_w, a_w);
            // l1
            double l1 = log_pi_S1;
            for (int k = 0; k <= max_I; ++k) {
            int n0 = 0, n1 = 0;
            for (int c = 0; c < C; ++c)
            if (I(c) == k) {
            if (Z(c, g)) n1 ++;
            else n0 ++;
            }
            l1 += R::lbeta(a_t + n1, a_t + n0) - R::lbeta(a_t, a_t);
            // Rcout << "k = " << k << ", n0 = " << n0 << ", n1 = " << n1 << ", ll = " << R::lbeta(a_t + n1, a_t + n0) - R::lbeta(a_t, a_t) << std::endl;
            }
            // Rcout << "l0 = " << l0 << ", l1 = " << l1 << std::endl;
            // log_add
            if (l0 > l1)
            ll += l0 + log(1 + exp(l1 - l0));
            else
            ll += l1 + log(1 + exp(l0 - l1));
            }
            return ll;
            }
            ')

