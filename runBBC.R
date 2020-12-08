library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)
source("BBC1/R/functions/basic.R")
sourceCpp('BBC1/src/api_cluster.cpp', cacheDir='cppCache')
source('BBC2_HBBC/bicluster_gibbs.R')
sourceCpp('BBC2_HBBC/bicluster_gibbs_uniform_prior.cpp')

dat = as.matrix(read.table('simuData.txt'))

#run BBC1
colnames(dat)=as.character(seq(1,ncol(dat)))
row.names(dat)=as.character(seq(1,nrow(dat)))
res_I_1 = basic(dat, K_from=4, K_to=4, n_chain=2, n_burn=100, n_mcmc=200, a_t=1)
##access row clustering result
#res_I_1$K4$I_hat

#run BBC2
colnames(dat) = c()
M = max(dat)+1
GAMMA_vec = rep(1, M)
ALPHA = 0.05
BETA = 0.05
LogB_GAMMA_vec = LogB(GAMMA_vec)

K=4
result_cpp = gibbs_bicluster_cpp(dat, K = K, M = M, iter_num = 400,burnin = 200, resample_Z = TRUE,resample_Z_until = 200)
##access row clustering result
#result_cpp$Zhat

#run HBBC
result_tree = gibbs_bicluster_tree(dat, M, iter_num = 400, burnin = 200, min_count = 10, max_K = 50, resample_Z = TRUE)
##hierarchical clustering result saved in result_tree$classification
