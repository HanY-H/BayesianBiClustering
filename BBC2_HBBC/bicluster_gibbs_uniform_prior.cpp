//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
vec rdirichlet_cpp(vec alpha, int L){
  vec gamma_draw(L);
  double sum_draw = 0;
  for(int i=0;i<L;i++){
    double draw = R::rgamma(alpha[i], 1.0);
    //printf("%f,", draw);
    gamma_draw[i] = draw;
    sum_draw = sum_draw + draw;
  }
  vec res(L);
  for(int i=0;i<L;i++){
    res[i] = gamma_draw[i]/sum_draw;
  }
  return res;
}



// [[Rcpp::export]]
double sum_cpp(vec x, int n){
  double tmp = 0;
  //int n = x.size();
  for(int j = 0; j<n; j++){
    tmp = tmp + x[j];
  }
  return tmp;
}

// [[Rcpp::export]]
vec compute_prob(vec a, int L){
  double b = max(a);
  vec a0(L);
  for(int i=0;i<L;i++){
    a0[i] = exp(a[i] - b);
  }
  double sum_a0 = sum_cpp(a0, L);
  vec a1(L);
  for(int i=0;i<L;i++){
    a1[i] = a0[i]/sum_a0;
  }
  return a1;
}

// [[Rcpp::export]]
double log_sum_exp(vec a, int L){
  double b = max(a);
  vec a0(L);
  for(int i=0;i<L;i++){
    a0[i] = exp(a[i] - b);
  }
  double sum_a0 = sum_cpp(a0, L);
  return log(sum_a0) + b;
}
  
  

// [[Rcpp::export]]
int sum_int_cpp(ivec x, int n){
  int tmp = 0;
  //int n = x.size();
  for(int j = 0; j<n; j++){
    tmp = tmp + x[j];
  }
  return tmp;
}

// [[Rcpp::export]]  
vec lgamma_vec_cpp(vec x, int n){
  //int n = x.size();
  vec res(n);
  for(int j = 0; j < n; j++){
    res[j] = lgamma(x[j]);
  }
  return res;
}  

// [[Rcpp::export]]
double lgamma_cpp(double x){
  double res;
  res = lgamma(x);
  return res;
}

// [[Rcpp::export]]
double LogB_cpp(vec x, int n){
  vec lgamma_x_vec = lgamma_vec_cpp(x, n);
  double up = sum_cpp(lgamma_x_vec, n);
  double sum_x = sum_cpp(x, n);
  double down = lgamma_cpp(sum_x);
  return up - down;
}

// [[Rcpp::export]]
double likelihood_cpp(imat S, int sum_S, ivec C, ivec Z, icube H1, imat H0, imat dat, int J, int M, int K, 
                      double LogB_GAMMA_vec, vec GAMMA_vec){

  double lfac1 = -(sum_S + J)*LogB_GAMMA_vec;
  double lfac2 = 0;
  double lfac3 = 0;
  for(int j=0; j<J; j++){
    ivec tmp_sum(M, fill::zeros);
    for(int k=0; k<K; k++){
      if(S(k,j) == 1){
        ivec tmp_count = H1.slice(j).col(k);
        lfac2 = lfac2 + LogB_cpp(GAMMA_vec + tmp_count, M);
        tmp_sum = tmp_sum + tmp_count;
      }
    }
    ivec H0j = H0.col(j);
    lfac3 = lfac3 + LogB_cpp(GAMMA_vec + H0j, M);
  }
  return lfac1 + lfac2 + lfac3;
}
// [[Rcpp::export]]
ivec column_sum(imat S){
  ivec a = arma::conv_to<arma::ivec>::from(sum(S, 0));
  return a;
}


// [[Rcpp::export]]
double prior_SZ_cpp(imat S, int sum_S, int K, int I, int J, double ALPHA, double BETA){
  double fac1 = -I*log(K) + (K-1)*log(ALPHA) - lgamma(K);
  double fac2 = sum_S*log(BETA) + (K*J-sum_S)*log(1-BETA);
  ivec S_col_sum = column_sum(S);
  int S_col_sum_0 = 0;
  for(int j = 0; j<J; j++){
    if(S_col_sum[j] == K){
      S_col_sum_0++;
    }
  }
  double fac3 = S_col_sum_0*(log(K-(K-1)*BETA) - log(BETA));
  return fac1 + fac2 + fac3;
}

// [[Rcpp::export]]
double prior_Z_cpp(int K, int I, double ALPHA){
  double fac1 = -I*log(K)  + (K-1)*log(ALPHA) - lgamma(K);
  return fac1;
}

// [[Rcpp::export]]
double posterior_bruteforce(imat S, ivec Z, imat dat, int I, int J, int M, int K, 
                            double LogB_GAMMA_vec, vec GAMMA_vec, double ALPHA, double BETA){
  icube H1(M, K, J, fill::zeros);
  imat H0(M, J, fill::zeros);
  ivec C(I, fill::zeros);
  for(int i=0;i<I;i++){
    int k = Z[i];
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      H1(obs, k, j) = H1(obs, k, j) + 1;
      if(S(k, j) == 0){
        H0(obs, j) = H0(obs, j) + 1;
      }
    }
    C[k] = C[k] + 1;
  }
  int sum_S = 0;
  for(int k=0; k<K; k++){
    for(int j=0; j<J; j++){
      sum_S = sum_S + S(k,j);
    }
  }
  double lik = likelihood_cpp(S, sum_S, C, Z, H1, H0, dat, J, M, K, LogB_GAMMA_vec, GAMMA_vec);
  double prior = prior_SZ_cpp(S, sum_S, K, I, J, ALPHA, BETA);
  return lik + prior;
}

// [[Rcpp::export]]
double likelihood_bruteforce(imat S, ivec Z, imat dat, int I, int J, int M, int K, 
                            double LogB_GAMMA_vec, vec GAMMA_vec){
  icube H1(M, K, J, fill::zeros);
  imat H0(M, J, fill::zeros);
  ivec C(I, fill::zeros);
  for(int i=0;i<I;i++){
    int k = Z[i];
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      H1(obs, k, j) = H1(obs, k, j) + 1;
      if(S(k, j) == 0){
        H0(obs, j) = H0(obs, j) + 1;
      }
    }
    C[k] = C[k] + 1;
  }
  int sum_S = 0;
  for(int k=0; k<K; k++){
    for(int j=0; j<J; j++){
      sum_S = sum_S + S(k,j);
    }
  }
  double lik = likelihood_cpp(S, sum_S, C, Z, H1, H0, dat, J, M, K, LogB_GAMMA_vec, GAMMA_vec);
  return lik;
}

// [[Rcpp::export]]
double posterior_one_cluster(imat dat, int I, int J, int M,  
                             double LogB_GAMMA_vec, vec GAMMA_vec){
  imat Hjs(M, J, fill::zeros);
  for(int i=0;i<I;i++){
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      Hjs(obs, j) = Hjs(obs, j) + 1;
    }
  }
  
  double fac = 0;
  for(int j=0;j<J;j++){
    fac = fac + LogB_cpp(Hjs.col(j) + GAMMA_vec, M);
  }
  double lik = -J*LogB_GAMMA_vec + fac;
  double prior = 0.0;
  return lik + prior;
}


// [[Rcpp::export]]
int sample_discrete(vec log_prob, int L){
  double max_log_prob = max(log_prob);
  vec prob_adj = vec(L);
  for(int ll=0; ll<L; ll++){
    prob_adj[ll] = exp(log_prob[ll] - max_log_prob);
  }
  ivec targets = Range(0, L-1);
  ivec sim_vec = Rcpp::RcppArmadillo::sample(targets, 1, true, prob_adj);
  int sim = sim_vec[0];
  return sim;
}

// [[Rcpp::export]]
int sample_discrete_equal_weight(int L){
  ivec targets = Range(0, L-1);
  ivec sim_vec = Rcpp::RcppArmadillo::sample(targets, 1, false);
  int sim = sim_vec[0];
  return sim;
}

// [[Rcpp::export]]
imat modify_S(imat S, int K, int J){
  imat Stmp(S);
  for(int j=0;j<J;j++){
    int col_j_sum = sum_int_cpp(S.col(j), K);
    if(col_j_sum == K-1){
      ivec S1(K, fill::ones);
      Stmp.col(j) = S1;
    }
  }
  return Stmp;
}

// [[Rcpp::export]]
Rcpp::List generate_H(imat dat, ivec Ztmp, imat Stmp, int I, int J, int M, int K){
  icube H1tmp(M, K, J, fill::zeros);
  imat H0tmp(M, J, fill::zeros);
  ivec Ctmp(I, fill::zeros);
  ivec C0tmp(J, fill::zeros);
  imat Hjs(M, J, fill::zeros);
  for(int i=0;i<I;i++){
    int k = Ztmp[i];
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      Hjs(obs, j) = Hjs(obs, j) + 1;
      H1tmp(obs, k, j) = H1tmp(obs, k, j) + 1;
      if(Stmp(k, j) == 0){
        H0tmp(obs, j) = H0tmp(obs, j) + 1;
        C0tmp[j] = C0tmp[j] + 1;
      }
    }
    Ctmp[k] = Ctmp[k] + 1;
  }
  return Rcpp::List::create(Rcpp::Named("H1") = H1tmp, Rcpp::Named("H0") = H0tmp,
                            Rcpp::Named("C") = Ctmp, Rcpp::Named("C0") = C0tmp,
                            Rcpp::Named("Hjs") = Hjs);

}



// [[Rcpp::export]]
ivec update_Z_2(ivec Zinit, imat dat, int iter_num, int I, int J, int M, vec GAMMA_vec){
  ivec Ztmp(Zinit);
  int K = 2;
  icube H1tmp(M, K, J, fill::zeros);
  ivec Ctmp(I, fill::zeros);
  imat Hjs(M, J, fill::zeros);
  for(int i=0;i<I;i++){
    int k = Ztmp[i];
    
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      Hjs(obs, j) = Hjs(obs, j) + 1;
      H1tmp(obs, k, j) = H1tmp(obs, k, j) + 1;
    }
    Ctmp[k] = Ctmp[k] + 1;
  }
  
  double sum_GAMMA_vec = sum_cpp(GAMMA_vec, M);
  
  
  for(int n = 0; n < iter_num; n++){
    
    for(int i = 0; i < I; i++){
      int k1 = Ztmp[i];
      if(Ctmp[k1] > 1){
        vec weight_z = vec(K);
        weight_z[k1] = 0;
        for(int k2=0; k2<K; k2++){
          if(k2 != k1){
            double fac1 = 0;
            double fac2 = 0;
            double fac5 = 0;
            
            for(int j=0; j<J; j++){
              int obs = dat(i,j);
              
              fac1 = fac1 + log(sum_GAMMA_vec + Ctmp[k1] - 1) -
                log(GAMMA_vec[obs] + H1tmp(obs, k1, j) - 1);
              fac2 = fac2 - log(sum_GAMMA_vec + Ctmp[k2]) +
                log(GAMMA_vec[obs] + H1tmp(obs, k2, j));
            }
            fac5 = log(Ctmp[k2]) - log(Ctmp[k1]-1);
            weight_z[k2] = fac1+fac2+fac5;
          }
        }
        
        int sim_k = sample_discrete(weight_z, K);
        
        if(sim_k != k1){
          Ztmp[i] = sim_k;
          Ctmp[k1] = Ctmp[k1] - 1;
          Ctmp[sim_k] = Ctmp[sim_k] + 1;
          for(int j = 0; j<J; j++){
            int obs = dat(i,j);
            H1tmp(obs, k1, j) = H1tmp(obs, k1, j)-1;
            H1tmp(obs, sim_k, j) = H1tmp(obs, sim_k, j)+1;
          }
        }
      }
    }
  }
  return Ztmp;
}





// [[Rcpp::export]]
ivec update_S_cpp(icube H1, imat Hjs, int J, int K, int M, double LogB_GAMMA_vec, vec GAMMA_vec, double temperature,
                  imat all_S_config, vec prior_all_s){
  ivec sim_m_all(J);
  for(int j=0; j<J; j++){
    int all_S_num = pow(2, K) - K;
    vec weight_s(all_S_num);
    for(int m=0; m < all_S_num; m++){
      ivec s = all_S_config.col(m);
      double tmp_w1 = 0;
      ivec tmp_sum(M, fill::zeros);
      for(int k=0; k < K; k++){
        if(s[k] == 1){
          ivec H1_k_j = H1.slice(j).col(k);
          tmp_w1 = tmp_w1 + LogB_cpp(H1_k_j + GAMMA_vec, M) - LogB_GAMMA_vec;
          tmp_sum = tmp_sum + H1_k_j;
        }
      }
      weight_s[m] = tmp_w1 + LogB_cpp(Hjs.col(j) - tmp_sum + GAMMA_vec, M);
    }
    weight_s = weight_s + prior_all_s;
    
    // if(j == 100){
    //   double pb_j = posterior_bruteforce(Stmp, Ztmp, dat, I, J, M, K, LogB_ALPHA, ALPHA_theta, log_ALPHA_crp, BETA);
    //   vec weight_s_2 = vec(all_S_num);
    //   for(int m=0; m<all_S_num; m++){
    //     imat Stmp_2(Stmp);
    //     Stmp_2.col(j) = all_S_config.col(m);
    //     double pb_j_m = posterior_bruteforce(Stmp_2, Ztmp, dat, I, J, M, K, LogB_ALPHA, ALPHA_theta, log_ALPHA_crp, BETA);
    //     weight_s_2[m] = pb_j_m - pb_j;
    //   }
    //   for(int m=0; m<all_S_num; m++){
    //     printf("%f,",weight_s[m]-weight_s_2[m]);
    //   }
    //   printf("\n");
    // }
    
    int sim_m = sample_discrete(weight_s * temperature, all_S_num);
    sim_m_all[j] = sim_m;
  }
  return sim_m_all;
}


// [[Rcpp::export]]
Rcpp::List gibbs_iter(imat Stmp, ivec Ztmp, imat dat,
                      int iter_num, int I, int J, int K, int M, imat all_S_config, vec prior_all_s,
                      double LogB_GAMMA_vec, vec GAMMA_vec, double ALPHA, double BETA, bool update_S = true,
                      bool resample_Z = true, int resample_Z_until = 100, int resample_Z_iter = 5, double temperature = 1.0){
  
  icube H1tmp(M, K, J, fill::zeros);
  imat H0tmp(M, J, fill::zeros);
  ivec Ctmp(I, fill::zeros);
  ivec C0tmp(J, fill::zeros);
  imat Hjs(M, J, fill::zeros);
  for(int i=0;i<I;i++){
    int k = Ztmp[i];
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      Hjs(obs, j) = Hjs(obs, j) + 1;
      H1tmp(obs, k, j) = H1tmp(obs, k, j) + 1;
      if(Stmp(k, j) == 0){
        H0tmp(obs, j) = H0tmp(obs, j) + 1;
        C0tmp[j] = C0tmp[j] + 1;
      }
    }
    Ctmp[k] = Ctmp[k] + 1;
  }
  
  
  imat simZ(I, iter_num);
  icube simS(K,J,iter_num);
  vec sim_lik(iter_num);
  vec sim_post(iter_num);
  double sum_GAMMA_vec = sum_cpp(GAMMA_vec, M);
  
  //iter_num = 1;
  
  int resample_index = -2;
  
  for(int n = 0; n < iter_num; n++){
    
    //Sample Z
    //printf("temp = %f ", temperature);
    //printf("n = %i: ", n);
    //printf("row_simS = %i: ", simS.n_rows);
    
    // bool zero_cluster = false;
    // for(int kk = 0;kk < K;kk++){
    //   printf("%i, ", Ctmp[kk]);
    //   if(Ctmp[kk] == 0){
    //     zero_cluster = true;
    //   }
    // }
    
    //if(resample_Z && n <= resample_Z_until && (zero_cluster || (n > 0 && (n % 10) == 0))){
    if(resample_Z && n <= resample_Z_until && ((n > 0 && (n % 10) == 0))){
      
      //printf("\nresampling Z:");
      
      ivec Ztmp_new(Ztmp);
      //imat Stmp_new(K+1, J);
      //for(int kk=0;kk<K;kk++){
      //  Stmp_new.row(kk) = Stmp.row(kk);
      //}
      int kk = sample_discrete_equal_weight(K);
      while(Ctmp[kk] <= 1){
        kk = sample_discrete_equal_weight(K);
      }

      //printf("split:%i\n", kk);
      uvec ids_kk = find(Ztmp == kk);
      //printf("C_k2:%i,len_ids_k2:%i\n", Ctmp[k2], ids_k2.size());
      
      for(int kkk = 0; kkk<K;kkk++){
        uvec ids_kkk = find(Ztmp == kkk);
        //printf("C_k:%i,len_ids_k:%i\n", Ctmp[kkk], ids_kkk.size());
      }
      
      imat dat_kk = dat.rows(ids_kk);
      ivec Ztmp_kk(Ctmp[kk], fill::zeros);
      
      ivec range_kk = Range(0, Ctmp[kk]-1);
      ivec shuffle_kk = shuffle(range_kk);
      
      for(int ii = Ctmp[kk]/2; ii < Ctmp[kk];ii++){
        Ztmp_kk[shuffle_kk[ii]] = 1;
      }
      
      //printf("good here");
      
      ivec Ztmp_kk_update = update_Z_2(Ztmp_kk, dat_kk, resample_Z_iter, Ctmp[kk], J, M, GAMMA_vec);
      //printf("good here2");
      
      for(int ii = 0; ii < Ctmp[kk];ii++){
        if(Ztmp_kk_update[ii] == 1){
          Ztmp_new[ids_kk[ii]] = K;
        }
      }
      //for(int j=0; j<J; j++){
      //  Stmp_new(kk,j) = 1;
      //  Stmp_new(K,j) = 1;
      //}
      
      imat merge_id_table(2, (K*(K+1))/2);
      vec merge_id_weight((K*(K+1))/2);
      int id_tmp = 0;
      for(int k1=0;k1<K;k1++){
        for(int k2=k1+1;k2<K+1;k2++){
          ivec Ztmp_new_copy(Ztmp_new);
          merge_id_table(0,id_tmp) = k1;
          merge_id_table(1,id_tmp) = k2;
          
          if(k2 != K){
            uvec ids_k2 = find(Ztmp_new == k2);
            uvec ids_K = find(Ztmp_new == K);
            Ztmp_new_copy.elem(ids_k2).fill(k1);
            Ztmp_new_copy.elem(ids_K).fill(k2);
          }else{
            uvec ids_k2 = find(Ztmp_new == k2);
            Ztmp_new_copy.elem(ids_k2).fill(k1);
          }
          
          H1tmp.fill(0);
          H0tmp.fill(0);
          Ctmp.fill(0);
          C0tmp.fill(0);
          for(int ii=0;ii<I;ii++){
            int kkk = Ztmp_new_copy[ii];
            for(int j=0;j<J;j++){
              int obs = dat(ii,j);
              H1tmp(obs, kkk, j) = H1tmp(obs, kkk, j) + 1;
            }
            Ctmp[kkk] = Ctmp[kkk] + 1;
          }
          
          imat Stmp_new(K,J,fill::ones);
          ivec sim_m_all = update_S_cpp(H1tmp, Hjs, J, K, M, LogB_GAMMA_vec, GAMMA_vec, temperature, all_S_config, prior_all_s);
          for(int j=0;j<J;j++){
            int sim_m = sim_m_all[j];
            ivec sim_s = all_S_config.col(sim_m);
            Stmp_new.col(j) = sim_s;
            ivec H0_tmp_sum(M, fill::zeros);
            for(int k=0; k < K; k++){
              if(sim_s[k] == 0){
                H0_tmp_sum = H0_tmp_sum + H1tmp.slice(j).col(k);
              }
            }
            H0tmp.col(j) = H0_tmp_sum;
            C0tmp[j] = sum_int_cpp(H0_tmp_sum, M);
          }
          
          int sum_Stmp_new = 0;
          for(int k=0; k<K; k++){
            for(int j=0; j<J; j++){
              sum_Stmp_new = sum_Stmp_new + Stmp_new(k,j);
            }
          }
          double sim_lik_tmp = likelihood_cpp(Stmp_new, sum_Stmp_new, Ctmp, Ztmp_new_copy, H1tmp, H0tmp, dat, J, M, K, 
                                              LogB_GAMMA_vec, GAMMA_vec);
          double sim_prior_tmp;
          if(update_S){
            sim_prior_tmp = prior_SZ_cpp(Stmp_new, sum_Stmp_new, K, I, J, ALPHA, BETA);
          }else{
            sim_prior_tmp = prior_Z_cpp(K, I, ALPHA);
          }
          merge_id_weight[id_tmp] = sim_lik_tmp + sim_prior_tmp;
          id_tmp++;
          }
        }
      int sample_id = sample_discrete(merge_id_weight, (K*(K+1))/2);
      
      int k1 = merge_id_table(0, sample_id);
      int k2 = merge_id_table(1, sample_id);
      
      ivec Ztmp_new_copy(Ztmp_new);
      
      if(k2 != K){
        uvec ids_k2 = find(Ztmp_new == k2);
        uvec ids_K = find(Ztmp_new == K);
        Ztmp_new_copy.elem(ids_k2).fill(k1);
        Ztmp_new_copy.elem(ids_K).fill(k2);
      }else{
        uvec ids_k2 = find(Ztmp_new == k2);
        Ztmp_new_copy.elem(ids_k2).fill(k1);
      }
      
      Ztmp = Ztmp_new_copy;
      
      H1tmp.fill(0);
      H0tmp.fill(0);
      Ctmp.fill(0);
      C0tmp.fill(0);
      for(int ii=0;ii<I;ii++){
        int kk = Ztmp[ii];
        for(int j=0;j<J;j++){
          int obs = dat(ii,j);
          H1tmp(obs, kk, j) = H1tmp(obs, kk, j) + 1;
        }
        Ctmp[kk] = Ctmp[kk] + 1;
      }
      
      imat Stmp_new(K,J,fill::ones);
      ivec sim_m_all = update_S_cpp(H1tmp, Hjs, J, K, M, LogB_GAMMA_vec, GAMMA_vec, temperature, all_S_config, prior_all_s);
      for(int j=0;j<J;j++){
        int sim_m = sim_m_all[j];
        ivec sim_s = all_S_config.col(sim_m);
        Stmp_new.col(j) = sim_s;
        ivec H0_tmp_sum(M, fill::zeros);
        for(int k=0; k < K; k++){
          if(sim_s[k] == 0){
            H0_tmp_sum = H0_tmp_sum + H1tmp.slice(j).col(k);
          }
        }
        H0tmp.col(j) = H0_tmp_sum;
        C0tmp[j] = sum_int_cpp(H0_tmp_sum, M);
      }
        
      Stmp = Stmp_new;
      
      //printf("After resampling:");
      for(int kkk = 0;kkk<K;kkk++){
        uvec ids_kkk = find(Ztmp == kkk);
        //printf("C_k:%i,len_ids_k:%i\n", Ctmp[kkk], ids_kkk.size());
      }
      //printf("\n");
    }
    
    for(int i = 0; i < I; i++){
      //printf("i:%i\n",i);
      for(int kk = 0; kk<K;kk++){
        ivec ids_kk = arma::conv_to<ivec>::from(find(Ztmp == kk));
        if(i == resample_index + 1 && int(ids_kk.size()) != Ctmp[kk]){
          //printf("error:i:%i,kk:%i,C_k2:%i,len_ids_k2:%i\n", i, kk, Ctmp[kk], ids_kk.size());
        }
      }
      //for(int i = 0; i < 1; i++){
      int k1 = Ztmp[i];
      
        
        
        
      
      vec weight_z = vec(K);
      weight_z[k1] = 0;
      for(int k2=0; k2<K; k2++){
        if(k2 != k1){
          double fac1 = 0;
          double fac2 = 0;
          double fac3 = 0;
          double fac4 = 0;
          double fac5 = 0;
          
          for(int j=0; j<J; j++){
            int obs = dat(i,j);
            if(Stmp(k1,j) == 1){
              fac1 = fac1 + log(sum_GAMMA_vec + Ctmp[k1] - 1) -
                log(GAMMA_vec[obs] + H1tmp(obs, k1, j) - 1);
            }
            if(Stmp(k2,j) == 1){
              fac2 = fac2 - log(sum_GAMMA_vec + Ctmp[k2]) +
                log(GAMMA_vec[obs] + H1tmp(obs, k2, j));
            }
            if(Stmp(k1, j) == 1 && Stmp(k2, j) == 0){
              fac3 = fac3 + log(H0tmp(obs,j) + GAMMA_vec[obs])-
                log(C0tmp[j] + sum_GAMMA_vec);
            }
            if(Stmp(k1, j) == 0 && Stmp(k2, j) == 1){
              fac4 = fac4 - log(H0tmp(obs,j) + GAMMA_vec[obs] - 1)+
                log(C0tmp[j] + sum_GAMMA_vec - 1);
            }
          }
          
          weight_z[k2] = fac1+fac2+fac3+fac4+fac5;
        }
      }
      
      // if(i == 100){
      //   double pb_i = posterior_bruteforce(Stmp, Ztmp, dat, I, J, M, K, LogB_ALPHA, ALPHA_theta, log_ALPHA_crp, BETA);
      //   vec weight_z_2 = vec(K);
      //   for(int k2=0; k2<K; k2++){
      //     ivec Ztmp_2(Ztmp);
      //     Ztmp_2[i] = k2;
      //     double pb_i_k2 = posterior_bruteforce(Stmp, Ztmp_2, dat, I, J, M, K, LogB_ALPHA, ALPHA_theta, log_ALPHA_crp, BETA);
      //     weight_z_2[k2] = pb_i_k2 - pb_i;
      //   }
      //   for(int k2=0; k2<K; k2++){
      //     printf("%f,",weight_z[k2]);
      //     printf("%f,",weight_z_2[k2]);
      //   }
      //   printf("\n");
      // }
      
      int sim_k = sample_discrete(weight_z * temperature, K);
      // int sim_k;
      // if(k1 == 2){
      //   sim_k = 1;
      // }else{
      //   sim_k = 2;
      // }
      
      if(sim_k != k1){
        Ztmp[i] = sim_k;
        Ctmp[k1] = Ctmp[k1] - 1;
        Ctmp[sim_k] = Ctmp[sim_k] + 1;
        for(int j = 0; j<J; j++){
          int obs = dat(i,j);
          H1tmp(obs, k1, j) = H1tmp(obs, k1, j)-1;
          H1tmp(obs, sim_k, j) = H1tmp(obs, sim_k, j)+1;
          if(Stmp(k1,j) == 0){
            H0tmp(obs, j) = H0tmp(obs, j) - 1;
            C0tmp[j] = C0tmp[j] - 1;
          }
          if(Stmp(sim_k,j) == 0){
            H0tmp(obs, j) = H0tmp(obs, j) + 1;
            C0tmp[j] = C0tmp[j] + 1;
          }
        }
      }
    }
    
    
    simZ.col(n) = Ztmp;
    

    
    //Sample S
    if(update_S){
      ivec sim_m_all = update_S_cpp(H1tmp, Hjs, J, K, M, LogB_GAMMA_vec, GAMMA_vec, temperature, all_S_config, prior_all_s);
      for(int j=0;j<J;j++){
        int sim_m = sim_m_all[j];
        ivec sim_s = all_S_config.col(sim_m);
        Stmp.col(j) = sim_s;
        simS.slice(n).col(j) = sim_s;
        ivec H0_tmp_sum(M, fill::zeros);
        for(int k=0; k < K; k++){
          if(sim_s[k] == 0){
            H0_tmp_sum = H0_tmp_sum + H1tmp.slice(j).col(k);
          }
        }
        H0tmp.col(j) = H0_tmp_sum;
        C0tmp[j] = sum_int_cpp(H0_tmp_sum, M);
      }
    }
    
    
    int sum_Stmp = 0;
    for(int k=0; k<K; k++){
      for(int j=0; j<J; j++){
        sum_Stmp = sum_Stmp + Stmp(k,j);
      }
    }
    double sim_lik_tmp = likelihood_cpp(Stmp, sum_Stmp, Ctmp, Ztmp, H1tmp, H0tmp, dat, J, M, K, LogB_GAMMA_vec, GAMMA_vec);
    double sim_prior_tmp;
    if(update_S){
      sim_prior_tmp = prior_SZ_cpp(Stmp, sum_Stmp, K, I, J, ALPHA, BETA);
    }else{
      sim_prior_tmp = prior_Z_cpp(K, I, ALPHA);
    }
    sim_lik[n] = sim_lik_tmp;
    sim_post[n] = sim_lik_tmp + sim_prior_tmp;
    
    //printf("posterior = %f\n", sim_post[n]);
  }
  return Rcpp::List::create(Rcpp::Named("sim_lik") = sim_lik,
                            Rcpp::Named("sim_post") = sim_post,
                            Rcpp::Named("simS") = simS,
                            Rcpp::Named("simZ") = simZ,
                            Rcpp::Named("H1tmp") = H1tmp, Rcpp::Named("H0tmp") = H0tmp,
                            Rcpp::Named("Ctmp") = Ctmp, Rcpp::Named("C0tmp") = C0tmp,
                            Rcpp::Named("Ztmp") = Ztmp, Rcpp::Named("Stmp") = Stmp,
                            Rcpp::Named("Hjs") = Hjs);
}

// [[Rcpp::export]]
double log_vec_power(vec a, vec b, int L){
  double res = 0;
  for(int i=0;i<L;i++){
    res = res + b[i]*log(a[i]);
  }
  return res;
}

// [[Rcpp::export]]
bool all_equal(ivec a, ivec b, int L){
  int res = 0;
  for(int i=0;i<L;i++){
    if(a[i] != b[i]){
      res = 1;
      break;
    };
  }
  return res == 0;
}


// [[Rcpp::export]]
bool equivalent_clustering(ivec a, ivec b, int K){
  int res = 0;
  for(int k=0;k<K;k++){
    uvec ids = find(a == k);
    int ii = ids.size();
    if(ii > 1){
      ivec b_k = b.elem(ids);
      int b_k0 = b_k[0];
      for(int j=1;j<ii;j++){
        if(b_k[j] != b_k0){
          res = 1;
          break;
        }
      }
    }
    if(res == 1){
      break;
    }
  }
  return res == 0;
}

// [[Rcpp::export]]
ivec get_id_order(ivec Z, int K, int I){
  ivec order_Z(K);
  order_Z.fill(-1);
  int id_k = 0;
  ivec val_summary(K, fill::zeros);
  for(int i=0;i<I;i++){
    int zi = Z[i];
    int if_equal = 0;
    for(int k=0;k<K;k++){
      if(zi==order_Z[k]){
        if_equal = 1;
      }
    }
    if(if_equal == 0){
      order_Z[id_k] = zi;
      val_summary[zi] = 1;
      id_k++;
    }
  }
  for(int k=0;k<K;k++){
    if(val_summary[k]==0){
      order_Z[id_k] = k;
      id_k++;
    }
  }
  return order_Z;
}

// [[Rcpp::export]]
Rcpp::List adjust_SZ(imat simZ, icube simS, int I, int J, int K, int sim_num){
  imat simZ1(simZ);
  icube simS1(simS);
  for(int ii=0;ii < sim_num;ii++){
    //printf("%i", ii);
    ivec Ztmp(simZ.col(ii));
    imat Stmp(simS.slice(ii));
    ivec Ztmp_copy(Ztmp);
    imat Stmp_copy(Stmp);
    ivec order_Ztmp = get_id_order(Ztmp, K, I);
    for(int k=0;k<K;k++){
      uvec id_k = find(Ztmp_copy == order_Ztmp[k]);
      Ztmp.elem(id_k).fill(k);
      Stmp.row(k) = Stmp_copy.row(order_Ztmp[k]);
    }
    simZ1.col(ii) = Ztmp;
    simS1.slice(ii) = Stmp;
  }
  return Rcpp::List::create(Rcpp::Named("simZ") = simZ1, Rcpp::Named("simS") = simS1);
}


std::pair<imat, ivec> count_table(imat a, int I, int sim_num){
  imat vals(I, sim_num, fill::zeros);
  int pattern = 0;
  ivec conts(sim_num, fill::zeros);
  for(int k=0;k<sim_num;k++){
    ivec ak = a.col(k);
    int find = 0;
    int find_id = 0;
    for(int j=0;j<pattern;j++){
      if(all_equal(ak, vals.col(j), I)){
        find = 1;
        find_id = j;
      }
    }
    if(find == 0){
      vals.col(pattern) = ak;
      conts[pattern] = 1;
      pattern = pattern + 1;
    }else{
      conts[find_id] = conts[find_id] + 1;
    }
  }
  imat vals0 = vals.cols(0, pattern-1);
  ivec conts0 = conts.subvec(0, pattern-1);
  return std::make_pair(vals0, conts0);
}


// [[Rcpp::export]]
Rcpp::List evaluate_posterior_K_method_1(int K, ivec Zhat, imat Shat, imat simZ, icube simS, int iter_num,
                                imat dat, int I, int J, int M, double LogB_GAMMA_vec, vec GAMMA_vec, 
                                double ALPHA, double BETA, imat all_S_config, vec prior_all_s){
  
  ivec order_Zhat = get_id_order(Zhat, K, I);
  //vec res(iter_num);
  int all_s_count = pow(2, K)-K;
  ivec Shat_match(J);
  for(int j=0;j<J;j++){
    ivec sj_hat = Shat.col(j);
    int m_match = -1;
    for(int m=0;m < all_s_count;m++){
      ivec sm = all_S_config.col(m);
      if(all_equal(sm, sj_hat, K)){
        m_match = m;
        break;
      }
      
    }
    Shat_match[j] = m_match;
  }
  
  // for(int ii=0;ii<iter_num;ii++){
  //   ivec Ztmp(simZ.col(ii));
  //   ivec Ztmp_copy(Ztmp);
  //   ivec order_Ztmp = get_id_order(Ztmp, K, I);
  //   for(int k=0;k<K;k++){
  //     uvec id_k = find(Ztmp_copy == order_Ztmp[k]);
  //     Ztmp.elem(id_k).fill(order_Zhat[k]);
  //   }
  //   simZ.col(ii) = Ztmp;
  // }
  
  std::pair<imat, ivec> Z_counts = count_table(simZ, I, iter_num);
  imat Zpattern = Z_counts.first;
  ivec Zcount = Z_counts.second;
  
  int sum_Shat = 0;
  for(int k=0;k<K;k++){
    for(int j=0;j<J;j++){
      sum_Shat = sum_Shat + Shat(k,j);
    }
  }

  int total_pattern = (int)Zpattern.n_cols;
  vec res(total_pattern);

  for(int ss = 0;ss<total_pattern;ss++){
    ivec Ztmp = Zpattern.col(ss);
    double fac1 = -posterior_bruteforce(Shat, Ztmp, dat, I, J, M, K, LogB_GAMMA_vec, GAMMA_vec, ALPHA, BETA);
    double fac2 = -lgamma(K+1);

    icube H1tmp(M, K, J, fill::zeros);
    imat H0tmp(M, J, fill::zeros);
    ivec Ctmp(I, fill::zeros);
    imat Hjs(M, J, fill::zeros);
    for(int i=0;i<I;i++){
      int k = Ztmp[i];
      for(int j=0;j<J;j++){
        int obs = dat(i,j);
        H1tmp(obs, k, j) = H1tmp(obs, k, j) + 1;
        Hjs(obs, j) = Hjs(obs, j) + 1;
        if(Shat(k, j) == 0){
          H0tmp(obs, j) = H0tmp(obs, j) + 1;
        }
      }
      Ctmp[k] = Ctmp[k] + 1;
    }

    mat prob_S(all_s_count, J);

    for(int j=0; j<J; j++){
      for(int m=0; m < all_s_count; m++){
        ivec s = all_S_config.col(m);
        double tmp_w1 = 0;
        ivec tmp_sum(M, fill::zeros);
        for(int k=0; k < K; k++){
          if(s[k] == 1){
            ivec H1_k_j = H1tmp.slice(j).col(k);
            tmp_w1 = tmp_w1 + LogB_cpp(H1_k_j + GAMMA_vec, M) - LogB_GAMMA_vec;
            tmp_sum = tmp_sum + H1_k_j;
          }
        }
        double weight_s_m = tmp_w1 + LogB_cpp(Hjs.col(j) - tmp_sum + GAMMA_vec, M) + prior_all_s[m];
        prob_S(m, j) = weight_s_m;
      }
    }

    for(int j=0;j<J;j++){
      prob_S.col(j) = compute_prob(prob_S.col(j), all_s_count);
    }


    double prob_Shat = 0;
    for(int j=0;j<J;j++){
      prob_Shat = prob_Shat + log(prob_S(Shat_match[j], j));
    }
    //prob_Shat = log P(Shat|Zhat, theta, dat)

    double fac3 = prob_Shat;
    //printf("%f,%f,%f,\n", fac1, fac2, fac3);
    res[ss] = fac1 + fac2 + fac3;

  }
  return Rcpp::List::create(Rcpp::Named("val") = res, Rcpp::Named("conts") = Zcount);
}


// [[Rcpp::export]]
Rcpp::List evaluate_posterior_K_method_2(int K, ivec Zhat, imat Shat, imat simZ, icube simS, int iter_num,
                            imat dat, int I, int J, int M, double LogB_GAMMA_vec, vec GAMMA_vec, 
                            double ALPHA, double BETA, imat all_S_config, vec prior_all_s){
  //if(K == 1){
  //  return posterior_one_cluster(dat, I, J, M, LogB_GAMMA_vec, GAMMA_vec);
  //}else{
  int equi_count = 0;
  for(int ii = 0;ii<iter_num;ii++){
    if(equivalent_clustering(Zhat, simZ.col(ii), K)){
      equi_count++;
    }
  }
  double equiv_portion = (double)equi_count/((double)iter_num);
  double fac4 = -log(equiv_portion);
  
  
  
  int sum_Shat = 0;
  for(int k=0;k<K;k++){
    for(int j=0;j<J;j++){
      sum_Shat = sum_Shat + Shat(k,j);
    }
  }
  double fac1 = posterior_bruteforce(Shat, Zhat, dat, I, J, M, K, LogB_GAMMA_vec, GAMMA_vec, ALPHA, BETA);
  double fac2 = lgamma(K+1);
  
  icube H1hat(M, K, J, fill::zeros);
  imat H0hat(M, J, fill::zeros);
  ivec Chat(I, fill::zeros);
  imat Hjs(M, J, fill::zeros);
  for(int i=0;i<I;i++){
    int k = Zhat[i];
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      H1hat(obs, k, j) = H1hat(obs, k, j) + 1;
      Hjs(obs, j) = Hjs(obs, j) + 1;
      if(Shat(k, j) == 0){
        H0hat(obs, j) = H0hat(obs, j) + 1;
      }
    }
    Chat[k] = Chat[k] + 1;
  }
  
  
    
  int all_s_count = pow(2, K)-K;
  mat prob_S(all_s_count, J);
  
  for(int j=0; j<J; j++){
    for(int m=0; m < all_s_count; m++){
      ivec s = all_S_config.col(m);
      double tmp_w1 = 0;
      ivec tmp_sum(M, fill::zeros);
      for(int k=0; k < K; k++){
        if(s[k] == 1){
          ivec H1_k_j = H1hat.slice(j).col(k);
          tmp_w1 = tmp_w1 + LogB_cpp(H1_k_j + GAMMA_vec, M) - LogB_GAMMA_vec;
          tmp_sum = tmp_sum + H1_k_j;
        }
      }
      double weight_s_m = tmp_w1 + LogB_cpp(Hjs.col(j) - tmp_sum + GAMMA_vec, M) + prior_all_s[m];
      prob_S(m, j) = weight_s_m;
    }
  }
    
  for(int j=0;j<J;j++){
    prob_S.col(j) = compute_prob(prob_S.col(j), all_s_count);
  }
  
  
  double prob_Shat = 0;
  for(int j=0;j<J;j++){
    ivec sj_hat = Shat.col(j);
    int m_match = -1;
    for(int m=0;m<all_s_count;m++){
      ivec sm = all_S_config.col(m);
      if(all_equal(sm, sj_hat, K)){
        m_match = m;
        break;
      }
      
    }
    prob_Shat = prob_Shat + log(prob_S(m_match, j));
  }
  //prob_Shat = log P(Shat|Zhat, theta, dat)
  
  double fac3 = -prob_Shat;
  //printf("%f,\n", prob_Shat);
    
  return Rcpp::List::create(Rcpp::Named("post_K") = fac1 + fac2 + fac3 + fac4,
                            Rcpp::Named("equiv_count") = equi_count,
                            Rcpp::Named("equiv_portion") = equiv_portion);
}


// [[Rcpp::export]]
double evaluate_posterior_K_method_3(int K, ivec Zhat, imat Shat, imat simZ, icube simS, int iter_num,
                            imat dat, int I, int J, int M, double LogB_GAMMA_vec, vec GAMMA_vec,
                            double ALPHA, double BETA, imat all_S_config, vec prior_all_s){
  //if(K == 1){
  //  return posterior_one_cluster(dat, I, J, M, LogB_GAMMA_vec, GAMMA_vec);
  //}else{
  int sum_Shat = 0;
  for(int k=0;k<K;k++){
    for(int j=0;j<J;j++){
      sum_Shat = sum_Shat + Shat(k,j);
    }
  }
  double fac1 = posterior_bruteforce(Shat, Zhat, dat, I, J, M, K, LogB_GAMMA_vec, GAMMA_vec, ALPHA, BETA);
  double fac2 = lgamma(K+1);
  
  icube H1hat(M, K, J, fill::zeros);
  imat H0hat(M, J, fill::zeros);
  ivec Chat(I, fill::zeros);
  imat Hjs(M, J, fill::zeros);
  for(int i=0;i<I;i++){
    int k = Zhat[i];
    for(int j=0;j<J;j++){
      int obs = dat(i,j);
      H1hat(obs, k, j) = H1hat(obs, k, j) + 1;
      Hjs(obs, j) = Hjs(obs, j) + 1;
      if(Shat(k, j) == 0){
        H0hat(obs, j) = H0hat(obs, j) + 1;
      }
    }
    Chat[k] = Chat[k] + 1;
  }
  
  
  
  int all_s_count = pow(2, K)-K;
  mat prob_S(all_s_count, J);
  
  for(int j=0; j<J; j++){
    for(int m=0; m < all_s_count; m++){
      ivec s = all_S_config.col(m);
      double tmp_w1 = 0;
      ivec tmp_sum(M, fill::zeros);
      for(int k=0; k < K; k++){
        if(s[k] == 1){
          ivec H1_k_j = H1hat.slice(j).col(k);
          tmp_w1 = tmp_w1 + LogB_cpp(H1_k_j + GAMMA_vec, M) - LogB_GAMMA_vec;
          tmp_sum = tmp_sum + H1_k_j;
        }
      }
      double weight_s_m = tmp_w1 + LogB_cpp(Hjs.col(j) - tmp_sum + GAMMA_vec, M) + prior_all_s[m];
      prob_S(m, j) = weight_s_m;
    }
  }
  
  for(int j=0;j<J;j++){
    prob_S.col(j) = compute_prob(prob_S.col(j), all_s_count);
  }
  
  
  double prob_Shat = 0;
  for(int j=0;j<J;j++){
    ivec sj_hat = Shat.col(j);
    int m_match = -1;
    for(int m=0;m<all_s_count;m++){
      ivec sm = all_S_config.col(m);
      if(all_equal(sm, sj_hat, K)){
        m_match = m;
        break;
      }
      
    }
    prob_Shat = prob_Shat + log(prob_S(m_match, j));
  }
  //prob_Shat = log P(Shat|Zhat, theta, dat)
  
  double fac4 = -prob_Shat;
  


  vec fac3_vec(iter_num);

  //mat Theta0tmp(M, J);
  //cube Thetatmp(M, K, J);

  for(int iter = 0; iter<iter_num; iter++){
    mat Theta0(M, J);
    cube Theta(M, K, J);
    ivec Ztmp = simZ.col(iter);
    imat Stmp = simS.slice(iter);

    icube H1tmp(M, K, J, fill::zeros);
    imat H0tmp(M, J, fill::zeros);
    //ivec Ctmp(I, fill::zeros);
    for(int i=0;i<I;i++){
      int k = Ztmp[i];
      for(int j=0;j<J;j++){
        int obs = dat(i,j);
        H1tmp(obs, k, j) = H1tmp(obs, k, j) + 1;
        if(Shat(k, j) == 0){
          H0tmp(obs, j) = H0tmp(obs, j) + 1;
        }
      }
      //Ctmp[k] = Chat[k] + 1;
    }

    for(int j=0;j<J;j++){
      for(int k=0;k<K;k++){
        vec theta_kj(M);
        if(Stmp(k,j) == 1){
          theta_kj = rdirichlet_cpp(H1tmp.slice(j).col(k) + GAMMA_vec, M);
        }else{
          theta_kj = rdirichlet_cpp(GAMMA_vec, M);
        }
        Theta.slice(j).col(k) = theta_kj;
      }
      vec theta_0j(M);
      theta_0j = rdirichlet_cpp(H0tmp.col(j) + GAMMA_vec, M);
      Theta0.col(j) = theta_0j;
    }

    mat prob_Z(K, I);
    for(int i=0;i<I;i++){
      for(int k=0;k<K;k++){
        double prob_zik = 0;
        for(int j=0;j<J;j++){
          int obs = dat(i,j);
          if(Stmp(k,j)==1){
            prob_zik = prob_zik + log(Theta.slice(j)(obs,k));
          }else{
            prob_zik = prob_zik + log(Theta0(obs, j));
          }
        }
        prob_Z(k,i) = prob_zik;
      }
    }
    for(int i=0;i<I;i++){
      prob_Z.col(i) = compute_prob(prob_Z.col(i), K);
    }
    double prob_Zhat = 0;
    for(int i=0;i<I;i++){
      prob_Zhat = prob_Zhat + log(prob_Z(Zhat[i], i));
    }
    //prob_Zhat = log P(Zhat|S, theta, dat)


    //prob_Shat = log P(Shat|Zhat, theta, dat)

    fac3_vec[iter] = prob_Zhat;
    //printf("%i,%f\n", iter, prob_Zhat);

  }
  double fac3 = log(iter_num) - log_sum_exp(fac3_vec, iter_num);
  //printf("%f,%f,%f,", fac1, fac2, fac3);
  return fac1 + fac2 + fac3 + fac4;
   // return Rcpp::List::create(Rcpp::Named("post_K") = fac1 + fac2 + fac3,
   //                           Rcpp::Named("prob_S") = probS);
}

// [[Rcpp::export]]
imat post_mode_Z(imat simZ, int K, int I, int sim_num){
  mat res(I, I, fill::zeros);
  for(int jj = 0;jj<sim_num; jj++){
    ivec Z(simZ.col(jj));
    for(int i=0;i<I;i++){
      for(int j=0;j<I;j++){
        if(Z[i] == Z[j]){
          res(i,j) = res(i,j) + 1;
        }
      }
    }
  }
  double base = ((double)sim_num)/((double)K);
  imat res0(I, I, fill::zeros);
  for(int i=0;i<I;i++){
    for(int j=0;j<I;j++){
      if(res(i,j) > base){
        res0(i,j) = 1;
      }
    }
  }
  return(res0);
}

