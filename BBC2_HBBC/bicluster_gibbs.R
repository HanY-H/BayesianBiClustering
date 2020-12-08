#### dat: 0 to M-1
#### Zinit: 0 to K-1

#### generate all possible configs of S.
generate_all_S = function(K){
  if(K == 1){
    return(matrix(1, 1, 1))
  }
  res = matrix(0, 1, K)
  for(i in 1:K){
    if(i != K-1){ #equivalent configurations, use (1,1...,1)
      com = combn(1:K, i) #generate index of 1s for k number of 1s out K.
      for(j in 1:ncol(com)){
        x = rep(0, K)
        x[com[,j]] = 1
        res = rbind(res, x)
      }
    }
  }
  rownames(res) = NULL
  res
}

#### returns the row ids for cluster k.
Rk <- function(k, Z){
  which(Z == k)
}

#### the beta function (log)
LogB <- function(alpha) {
  # to prevent data explosion, we use log_gamma instead of gamma function
  up <- sum(lgamma(alpha))
  # numerator: sum of log_gammas
  down <- lgamma(sum(alpha))
  # denominator: log_gamma of sum
  res <- up - down
  return(res)
}

#### the count table for dat[Rk_Value, j]
H <- function(Rk_Value, j, dat) {
  subset <- dat[Rk_Value, j] 
  tmp.vec <- table(subset)
  fnl.vec <- count0
  fnl.vec[names(tmp.vec)] <- tmp.vec
  return(fnl.vec)
}

#### likelihood function
likelihood <- function(S, C, Z, H1, H0, Hjs, dat){
  J = ncol(dat)
  lfac1 <- -(sum(S)+J)*LogB_ALPHA
  lfac2 <- 0
  lfac3 <- 0
  for(j in 1:J){
    tmp.sum <- rep(0,M)
    for(k in which(S[,j]==1)){
      tmp.count <- H1[,k,j] #H(R[[k]],j, dat = dat)
      lfac2 <- lfac2 + LogB(ALPHA_theta + tmp.count)
      tmp.sum <- tmp.sum + tmp.count
    }
    lfac3 <- lfac3 + LogB(ALPHA_theta + H0[,j])
  }
  return(lfac1+lfac2+lfac3)
}

prior_SZ = function(S, C, K, J){
  fac1 = K*log(ALPHA_crp) + sum(lgamma(C))
  fac2 = sum(S)*log(BETA) + (K*J - sum(S))*log(1-BETA)
  S_sum = apply(S, 2, sum)
  fac3 = sum(S_sum == K)*(log(K-(K-1)*BETA) - log(BETA)) 
  fac1 + fac2 + fac3
}

prior_S_j = function(S_j, K){
  fac2 = sum(S_j)*log(BETA) + (K - sum(S_j))*log(1-BETA)
  S_j_sum = sum(S_j)
  fac3 = (S_j_sum == K)*(log(K-(K-1)*BETA) - log(BETA))
  fac2 + fac3
}

###########################################
##################gibbs function
gibbs_bicluster = function(dat, K, M, Zinit = c(), Sinit = c(), iter_num = 1000){
  I = nrow(dat)
  J = ncol(dat)
  
  all_S_config = generate_all_S(K)
  ps = apply(all_S_config, 1, prior_S_j, K = K)
  
  simZ <- matrix(NA,I,iter_num)
  simS <- array(NA,c(K,J,iter_num))
  sim_lik <- rep(NA,iter_num)
  sim_post <- rep(NA,iter_num)
  
  Hjs <- matrix(0,M,J)
  row.names(Hjs) <- elements
  for(j in 1:J){
    Hjs[,j] <- H(1:I, j=j, dat = dat)
  }
  
  if(is.null(Zinit)){
    Ztmp <- sample(1:K, I, replace = TRUE)
  }else{
    Ztmp = Zinit
  }
  
  if(is.null(Sinit)){
    Stmp <- matrix(1,K,J)
  }else{
    Stmp = Sinit
  }
  
  simZ[,1] <- Ztmp #c(1:K,rep(0,(I-K)))
  simS[,,1] <- Stmp
  
  Rtmp = lapply(1:K, Rk, Z = Ztmp)
  Ctmp = sapply(Rtmp, length)
  
  H1tmp = array(NA, c(M, K, J))
  H0tmp = matrix(0, M, J)
  for(j in 1:J){
    for(k in 1:K){
      H1tmp[,k,j] = H(Rtmp[[k]], j, dat)
      if(Stmp[k, j] == 0){
        H0tmp[,j] = H0tmp[,j] + H1tmp[,k,j]
      }
    }
  }
  
  C0tmp = apply(H0tmp, 2, sum)
  
  for(n in 1:(iter_num-1)){
    print(n)
    ##########sample Z or R
    
    for(i in 1:I){
      if(Ctmp[Ztmp[i]] > 1){
        k1 <- Ztmp[i]
        weight.Z <- rep(NA,K)
        weight.Z[k1] <- 0
        for(k2 in 1:K){
          if(k2!=k1){
            
            fac1 = 0
            Sk1 = which(Stmp[k1,] == 1)
            Sk2 = which(Stmp[k2,] == 1)
            
            for(j in Sk1){
              obs = dat[i,j]
              fac1 = fac1 + log(sum(ALPHA_theta) + Ctmp[k1] - 1) -
                log(ALPHA_theta[obs] + H1tmp[obs, k1, j] - 1)
            }
            
            fac2 = 0
            for(j in Sk2){
              obs = dat[i,j]
              fac2 = fac2 - log(sum(ALPHA_theta) + Ctmp[k2]) +
                log(ALPHA_theta[obs] + H1tmp[obs, k2, j])
            }
            
            fac3 = 0
            for(j in setdiff(Sk1,Sk2)){
              obs = dat[i,j]
              fac3 = fac3 + log(H0tmp[obs,j] + ALPHA_theta[obs])-
                log(C0tmp[j] + sum(ALPHA_theta))
            }
            
            fac4 = 0
            for(j in setdiff(Sk2,Sk1)){
              obs = dat[i,j]
              fac4 = fac4 - log(H0tmp[obs,j] + ALPHA_theta[obs] - 1)+
                log(C0tmp[j] + sum(ALPHA_theta) - 1)
            }
            
            fac5 = log(Ctmp[k2]) - log(Ctmp[k1]-1)
            
            weight.Z[k2] <- fac1+fac2+fac3+fac4+fac5
          }
        }
        weight.Z.final <- weight.Z - max(weight.Z)
        sim.k <- sample(1:K,1,prob=exp(weight.Z.final))
        if(sim.k != k1){
          Ztmp[i] = sim.k
          Ctmp[k1] = Ctmp[k1]-1
          Ctmp[sim.k] = Ctmp[sim.k]+1
          for(j in 1:J){
            obs = dat[i,j]
            H1tmp[obs,k1,j] = H1tmp[obs,k1,j] - 1
            H1tmp[obs,sim.k,j] = H1tmp[obs,sim.k,j] + 1
            if(Stmp[k1,j] == 0){
              H0tmp[obs,j] = H0tmp[obs,j] - 1
              C0tmp[j] = C0tmp[j] - 1
            }
            if(Stmp[sim.k,j] == 0){
              H0tmp[obs,j] = H0tmp[obs,j] + 1
              C0tmp[j] = C0tmp[j] + 1
            }
          }
        }
      }
    }
    simZ[,n+1] <- Ztmp
    
    #########sample S
    for(j in 1:J){
      weight.S <- rep(NA,2^K-K)
      
      for(m in 1:(2^K-K)){
        s <- all_S_config[m,]
        tmp.w1 <- 0
        tmp.sum <- rep(0,M)
        for(k in which(s == 1)){
          tmp.w1 <- tmp.w1 + LogB(H1tmp[,k,j] + ALPHA_theta) - LogB_ALPHA
          tmp.sum = tmp.sum + H1tmp[,k,j]
        }
        weight.S[m] <- tmp.w1 + LogB(Hjs[,j] - tmp.sum + ALPHA_theta)
      }
      weight.S.final <- weight.S - max(weight.S) + ps
      sim.m <- sample(1:(2^K-K),1,prob=exp(weight.S.final))
      sim.s = all_S_config[sim.m,]
      Stmp[,j] = sim.s
      simS[,j,n+1] <- sim.s
      k0 = which(sim.s == 0)
      H0tmp[,j] = apply(H1tmp[,,j][,k0,drop = FALSE], 1, sum)
      C0tmp[j] = sum(H0tmp[,j])
    }
    sim_lik[n+1] <- likelihood(Stmp, Ctmp, Ztmp, H1tmp, H0tmp, Hjs, dat)
    sim_post[n+1] = sim_lik[n+1] + prior_SZ(Stmp, Ctmp, K, J)
  }
  list(sim_lik = sim_lik, sim_post = sim_post, simS = simS, simZ = simZ, H1tmp = H1tmp, H0tmp = H0tmp, C0tmp = C0tmp)
}


find_groups = function(A){
  z = rep(NA, nrow(A))
  node_left = 1:nrow(A)
  res = list()
  group_id = 1
  while(length(node_left) > 0){
    i = node_left[1]
    id_i = which(A[i,] == 1)
    res[[group_id]] = id_i
    z[id_i] = group_id - 1
    group_id = group_id + 1
    node_left = setdiff(node_left, id_i)
  }
  z
}


gibbs_bicluster_cpp = function(dat, K, M, Zinit = c(), Sinit = c(), iter_num = 1000, burnin = 300, post_K = TRUE, ...){
  I = nrow(dat)
  J = ncol(dat)
  
  all_S_config = generate_all_S(K)
  prior_all_s = apply(all_S_config, 1, prior_S_j, K = K)
  
  if(is.null(Zinit)){
    Ztmp <- sample(0:(K-1), I, replace = TRUE)
  }else{
    Ztmp = Zinit
  }
  
  if(is.null(Sinit)){
    Stmp <- matrix(1,K,J)
  }else{
    Stmp = Sinit
  }
  
  all_S_config = t(all_S_config)
  
  A = gibbs_iter(Stmp, Ztmp, dat,
                 iter_num, I, J, K, M, all_S_config, prior_all_s, 
                 LogB_GAMMA_vec, GAMMA_vec, ALPHA, BETA, ...)
  
  simZ_all = A$simZ
  simS_all = A$simS
  
  B = adjust_SZ(simZ_all, simS_all, I, J, K, iter_num)
  A$simZ = B$simZ
  A$simS = B$simS
  
  sim_post = A$sim_post[(burnin+1):iter_num]
  simZ = A$simZ[,(burnin+1):iter_num]
  simS = A$simS[,,(burnin+1):iter_num]
  
  max_post = max(sim_post)
  max_id = which(sim_post == max_post)[1]
  Zhat = simZ[,max_id]
  Shat = simS[,,max_id]
  
  Zmean = find_groups(post_mode_Z(simZ, K, I, ncol(simZ)))
  Smean = apply(simS, c(1,2), mean)
  
  if(post_K){
  
    post_K0 = evaluate_posterior_K_method_1(K, Zhat, Shat, simZ, simS, length(sim_post),
                                           dat, I, J, M, LogB_GAMMA_vec, GAMMA_vec, 
                                           ALPHA, BETA, all_S_config, prior_all_s)
    
    post_K = evaluate_posterior_K_method_2(K, Zhat, Shat, simZ, simS, length(sim_post),
                         dat, I, J, M, LogB_GAMMA_vec, GAMMA_vec, 
                         ALPHA, BETA, all_S_config, prior_all_s)
    
    post_K1 = evaluate_posterior_K_method_3(K, Zhat, Shat, simZ, simS, length(sim_post),
                                           dat, I, J, M, LogB_GAMMA_vec, GAMMA_vec, 
                                           ALPHA, BETA, all_S_config, prior_all_s)
    A$post_K0 = post_K0
    A$post_K = post_K
    A$post_K1 = post_K1
  }
  A$Zhat = Zhat
  A$Shat = Shat
  A$Zmean = Zmean
  A$Smean = Smean
  
  return(A)
}

gibbs_bicluster_sa_cpp = function(dat, K, M, temps, Zinit = c(), Sinit = c(), sa_step = 30, iter_num = 1000, ...){
  I = nrow(dat)
  J = ncol(dat)
  sim_post = c()
  simZ = matrix(NA, I, iter_num)
  simS = array(NA, c(K, J, iter_num))
  index0 = 1
  for(i in 1:length(temps)){
    if(temps[i] < 1){
      gibbs_step = sa_step
    }else{
      gibbs_step = iter_num - sa_step * (length(temps)-1)
    }
    A = gibbs_bicluster_cpp(dat, K, M, Zinit, Sinit, iter_num = gibbs_step, temperature = temps[i], ...)

    sim_post = c(sim_post, A$sim_post)
    simZ[, seq(index0, length.out = gibbs_step)] = A$simZ
    simS[,,seq(index0, length.out = gibbs_step)] = A$simS
    Zinit = A$simZ[, gibbs_step]
    Sinit = A$simS[,,gibbs_step]
    index0 = index0 + gibbs_step
  }
  list(simZ = simZ, simS = simS, sim_post = sim_post)
}

gibbs_bicluster_tree = function(dat, M, iter_num = 1000, burnin = 100, min_count = 10, max_K = 10, resample_Z = TRUE, ...){
  signal_cols = list()
  signal_cols2 = list() 
  classification = list()
  sim_lik_tree = list()
  initial_post_seq = c()
  divide_post_seq = c()
  divide_order = c()
  J = ncol(dat)
  post_prob_K = list()
  result = list()
  result[[1]] = list()
  result[[1]]$R = 1:nrow(dat)
  result[[1]]$I = nrow(dat)
  bt=1 #P_split=p/(TLevel^bt)
  result[[1]]$TLevel=1 
  
  invalid_Z_mean = c()
  result_list = list()
  diff_post_list = list()
  
  result_list_EndOfEachStep=list() 
  p_split=0.5 #0.05
  for(k in 1:max_K){
    kk = length(result)
    post_prob_K[[k]] = list(init = rep(NA, kk), post_1 = rep(NA, kk), post_2 = rep(NA, kk), post_3 = rep(NA, kk))
    
    diff_post = rep(-Inf, kk) #w_r for the current leaves
    for(i in 1:kk){
      if(result[[i]]$I >= 2*min_count){
        if(is.null(result[[i]]$init_post)){
          result[[i]]$init_post =  posterior_one_cluster(dat[result[[i]]$R,], result[[i]]$I, J, M,  
                                                         LogB_GAMMA_vec, GAMMA_vec)
          A = gibbs_bicluster_cpp(dat[result[[i]]$R,], 2, M, iter_num = iter_num, burnin = burnin, 
                                  resample_Z = resample_Z, resample_Z_until = burnin, update_S = TRUE)
          #max_post = max(A$sim_post)
          #max_id = which(A$sim_post == max_post)[1]
          Z.hat = A$Zhat #MCMC run with max posterior
          S.hat = A$Shat #MCMC run with max posterior
          result[[i]]$Z.hat = Z.hat
          result[[i]]$S.hat = S.hat[1,]
          result[[i]]$Z.mean = A$Zmean
          result[[i]]$S.mean = A$Smean[1,] #average over MCMC runs
          result[[i]]$post_1 = get_post_K0(A)[[1]][1]
          result[[i]]$post_2 = A$post_K$post_K
          result[[i]]$post_3 = A$post_K1
          result[[i]]$divide_post = A$post_K1
          result[[i]]$sim_post = A$sim_post
        }
        diff_post[i] = result[[i]]$divide_post - result[[i]]$init_post+log(p_split/ALPHA/(1-p_split))

        post_prob_K[[k]]$init[i] = result[[i]]$init_post
        post_prob_K[[k]]$post_1[i] = result[[i]]$post_1
        post_prob_K[[k]]$post_2[i] = result[[i]]$post_2
        post_prob_K[[k]]$post_3[i] = result[[i]]$post_3
      }
    }
    
    diff_post_list[[k]] = diff_post
    
    result_list[[k]] = result #estimated results if each of the current node are splitted.
    
    if(k==1){ 
      result_list_EndOfEachStep[[1]]=result
    }
    
    if(max(diff_post) <= log(k)){
      break
    }else{
      ii=0 
      descreasing_order=order(diff_post,decreasing = TRUE) #position of Largest to Smallest element
      for(choose_r in 1: length(descreasing_order)){
        iii=descreasing_order[choose_r]
        if(diff_post[iii]>log(k) && sum(result[[iii]]$Z.hat==0)>=min_count && sum(result[[iii]]$Z.hat==1)>=min_count){
          ii=iii
          break
        }
      }
      
      if(ii>0){
      divide_order = c(divide_order, ii)
      signal_cols[[k]] = result[[ii]]$S.mean #### S.hat ###selected S
      signal_cols2[[k]] = result[[ii]]$S.hat
      classification[[k]] = list()
      if(max(result[[ii]]$Z.mean) != 1){
        classification[[k]][[1]] = result[[ii]]$R[which(result[[ii]]$Z.hat == 0)]
        classification[[k]][[2]] = result[[ii]]$R[which(result[[ii]]$Z.hat == 1)]
        invalid_Z_mean[k] = 1
      }else{
        classification[[k]][[1]] = result[[ii]]$R[which(result[[ii]]$Z.mean == 0)]
        classification[[k]][[2]] = result[[ii]]$R[which(result[[ii]]$Z.mean == 1)]
        invalid_Z_mean[k] = 0
      }
      sim_lik_tree[[k]] = result[[ii]]$sim_post
      
      initial_post_seq = c(initial_post_seq, result[[ii]]$init_post)
      divide_post_seq = c(divide_post_seq, result[[ii]]$divide_post)
      
      result[[ii]]$init_post = NULL ##ii and ii+1 will store the newly created two nodes, replacing the currrect parent node ii. 
      #current (ii+1):end will be moved to (ii+2):(end+1)
      result[[ii]]$divide_post = NULL
      result[[ii]]$R = classification[[k]][[1]]
      result[[ii]]$I = length(classification[[k]][[1]])
      result[[ii]]$TLevel=result[[ii]]$TLevel+1 
      #S,Z,Post information associated with ii is inherited from parent node, and should NOT be used.
      
      result[[kk+1]] = list()
      result[seq(ii+2, length.out = kk-ii)] = result[seq(ii+1, length.out = kk-ii)]
      
      result[[ii+1]] = list()
      result[[ii+1]]$init_post = NULL
      result[[ii+1]]$divide_post = NULL
      result[[ii+1]]$R = classification[[k]][[2]]
      result[[ii+1]]$I = length(classification[[k]][[2]])
      result[[ii+1]]$TLevel=result[[ii]]$TLevel
      #S,Z,Post information will only be availabe in the next step.
      
      result_list_EndOfEachStep[[k+1]]=result 
      }else{
        break
      }
    }
  }
  
  list(signal_cols = signal_cols,signal_cols2=signal_cols2, classification = classification, sim_lik_tree = sim_lik_tree,
       initial_post_seq = initial_post_seq, divide_post_seq = divide_post_seq, divide_order = divide_order,
       result = result, post_prob_K = post_prob_K, diff_post_list = diff_post_list, result_list = result_list,
       invalid_Z_mean = invalid_Z_mean,result_list_New=result_list_EndOfEachStep)
}

get_post_K0 = function(result_cpp, min_p = 0.1){
  a0 = result_cpp$post_K0
  val = a0$val[,1]
  conts = a0$conts[,1]
  
  prop = conts/sum(conts)
  kk = which(prop > min_p)
  
  if(length(kk) == 0){
    res1 = NA
  }else{
    res1 = log(length(kk)) - log_sum_exp(val[kk]+log(prop[kk]), length(kk))
  }
  
  max_count = max(conts)[1]
  max_id = which(conts == max_count)
  res2 = -val[max_id] - log(max_count/sum(conts))
  list(c(res1, res2), res1-res2)
}