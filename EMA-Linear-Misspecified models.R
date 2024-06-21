#------------------------------------------------------
# Title: The EMA Algorithm - Example 1: a linear regression model.

# ---- Remove all variables ----
if(sys.nframe() == 0L){rm(list=ls()); gc()}

# address = ''
# setwd('')

# ---- Necessary packages ----
pacman::p_load(expm,mvtnorm,MASS,quadprog,mice,snowfall,parallel)

# ---- Functions required ----
source("MIS_covariate.R")
source("MIS.R")
source("density function_full.R")
source("density covariate fun.R")


# -------- Perform the EM algorithm to obtain estimates of \mu, \Sigma and \Pi --------
EM = function(data,          # dataset, including covariates and responses
              mu,            # initial value of mu (conditional mean vector of continuous variables given the value of discrete variables)
              Pi,            # initial value of Pi (the probability with which discrete variables take a certain value)
              Sigma.matrix   # initial value of Sigma.matrix (conditional covariance matrix of continuous variables given the value of discrete ones)
){
  
  # Initial parameter setting
  data_EM = data[-which(complete.cases(data[,c(2:5)])),]
  data_ori = data[which(complete.cases(data[,c(2:5)])),-1]
  mu.current = mu
  Pi.current = as.matrix(Pi)
  Sigma.matrix.current = Sigma.matrix
  mu.sum = matrix(c(0,0,1,1,0,1,0,1),nrow = 2, ncol = 4,byrow = TRUE)
  iter = 0
  
  while (iter < 50){
    data_gen = list()
    for (z in 1:dim(data_EM)[1]){
      data_gen[[z]] = MIS(data_EM[z,],mu = mu.current, Pi = Pi.current,Sigma.matrix = Sigma.matrix.current,
                          beta = c(1000,1000,1000,1000),sigma = 1000,model_index = 15)
    }
    
    mu.previous = mu.current
    Pi.previous = Pi.current
    Sigma.matrix.previous = Sigma.matrix.current
    
    # Calculate the estimates of mu_r and pi_r by the EM algorithm
    for (i in 1:4){
      E_f_mu = c(0,0)
      P_f_mu = 0
      for (j in 1:length(data_gen)){
        data_tem = matrix(unlist(data_gen[[j]]),nrow = 50,ncol = 4, byrow = FALSE)
        r.index = intersect(which(data_tem[,3] == mu.sum[1,i]),which(data_tem[,4] == mu.sum[2,i]))
        P_f_mu = P_f_mu + length(r.index)/50
        if (length(r.index) == 0){
          E_f_mu = E_f_mu
        }else if (length(r.index) == 1){
          E_f_mu = E_f_mu + data_tem[r.index,c(1,2)]/50
        }else{  
          E_f_mu = E_f_mu + colSums(data_tem[r.index,c(1,2)])/50
        }
      }
      r.index.ori = intersect(which(data_ori[,3] == mu.sum[1,i]),which(data_ori[,4] == mu.sum[2,i]))
      P_f_mu = P_f_mu + length(r.index.ori)
      if (length(r.index.ori) == 0){
        E_f_mu = E_f_mu
      }else if(length(r.index.ori) == 1){
        E_f_mu = E_f_mu + data_ori[r.index.ori,c(1,2)]
      }else{
        E_f_mu = E_f_mu + colSums(data_ori[r.index.ori,c(1,2)])
      }
      mu.current[,i] = E_f_mu/P_f_mu
      Pi.current[i,1] = P_f_mu/n
    }
    
    # Calculate the estimates of Omega by the EM algorithm
    Sigma.matrix.current = matrix(c(0,0,0,0),ncol = 2, nrow = 2)
    for (i in 1:4){
      for (j in 1:length(data_gen)){
        data_tem = matrix(unlist(data_gen[[j]]),nrow = 50,ncol = 4, byrow = FALSE)
        r.index = intersect(which(data_tem[,3] == mu.sum[1,i]),which(data_tem[,4] == mu.sum[2,i]))
        Sigma.matrix.current_tem = matrix(c(0,0,0,0),ncol = 2, nrow = 2)
        if (length(r.index) == 0){
          Sigma.matrix.current_tem = Sigma.matrix.current_tem
        }else{
          for (l in 1:length(r.index)){
            Sigma.matrix.current_tem = Sigma.matrix.current_tem + (data_tem[r.index[l],c(1,2)] - mu.current[,i])%*%t(data_tem[r.index[l],c(1,2)] - mu.current[,i])
          }
        }
        Sigma.matrix.current = Sigma.matrix.current + Sigma.matrix.current_tem/50
      }
      r.index.ori = intersect(which(data_ori[,3] == mu.sum[1,i]),which(data_ori[,4] == mu.sum[2,i]))
      if (length(r.index.ori) == 0){
        Sigma.matrix.current = Sigma.matrix.current
      }else{
        for (k in 1:length(r.index.ori)){
          Sigma.matrix.current = Sigma.matrix.current + (data_ori[r.index.ori[k],c(1,2)] - mu.current[,i])%*%t(data_ori[r.index.ori[k],c(1,2)] - mu.current[,i])
        }
      }
    }
    Sigma.matrix.current = Sigma.matrix.current/n
    
    # Stop criterion
    if ((norm(mu.previous - mu.current,"2")<0.01)&&(norm(Pi.previous - Pi.current,"2")<0.01)&&(norm(Sigma.matrix.previous - Sigma.matrix.current,"2")<0.01)){
      break
    }else{}
    
    iter = iter + 1
  }
  return(list(mu.current,Pi.current,Sigma.matrix.current))
}


# ----------------------------- E-MS algorithm -----------------------------
EMS = function(data,          # data set, responses of all subjects
               mu,            # estimate of mu produced by the EM algorithm
               Pi,            # estimate of Pi produced by the EM algorithm
               Sigma.matrix,  # estimate of covariance matrix produced by the EM algorithm
               beta,          # initial estimate of the coefficient
               sigma,         # initial estimate of the variance
               CandidateModel.list
){
  
  # Record computing time
  start_time_EMS = Sys.time()
  
  iter_EMS = 0
  model_index = 15 # Assume the initial model for the EMS algorithm is the full model
  number.candidatemodels = 15
  
  data_miss = data[-which(complete.cases(data)),]
  data_obs = data[which(complete.cases(data)),]
  
  while (iter_EMS < 50){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current model M_c
    data_gen_EMS = list()
    for (z in 1:dim(data_miss)[1]){
      data_gen_EMS[[z]] = MIS(data_miss[z,],mu = mu, Pi = Pi,Sigma.matrix = Sigma.matrix,
                              beta = beta,sigma = sigma,model_index = model_index)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma under different candidate models using the obtained full sample
    parameter.est = list()
    loss.summary = matrix(0,nrow = number.candidatemodels,ncol = 1)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 4){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 5)&&(k <= 10)){
        num_covariate = 2
        model_cur_loc = k - 4
      }else if ((k >= 11)&&(k <= 14)){
        num_covariate = 3
        model_cur_loc = k - 10
      }else{
        num_covariate = 4
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate,ncol = num_covariate)
      for (i in 1:length(data_gen_EMS)){
        data_gen_tem = matrix(unlist(data_gen_EMS[[i]]),nrow = 50, ncol = 5, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate,ncol = num_covariate)
        for (j in 1:50){
          S_1_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_gen_tem[j,1]
          S_2_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
            t(data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(n-length(data_gen_EMS))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
          t(data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
      }
      
      ## 2.3 Estimation of beta, sigma^2 and loss function
      beta_est = solve(S_2)%*%S_1
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/n
      
      RSS_est = 0
      for (i in 1:length(data_gen_EMS)){
        data_gen_tem = matrix(unlist(data_gen_EMS[[i]]),nrow = 50, ncol = 5, byrow = FALSE)
        RSS_est_tem = 0
        for (j in 1:50){
          RSS_est_tem = RSS_est_tem + (data_gen_tem[j,1] - t(data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc])+1])%*%beta_est)^2
        }
        RSS_est = RSS_est + RSS_est_tem/50
      }
      for (i in 1:(n-length(data_gen_EMS))){
        RSS_est = RSS_est + (data_obs[i,1] - t(data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])%*%beta_est)^2
      }
      loss_est = -n*log(sigma_square_est)/2 - 1/(2*sigma_square_est)*RSS_est
      parameter.est[[k]] = list(beta_est,sigma_square_est)
      loss.summary[k,1] = loss_est
    }
    
    # 3. Find the optimal model based on (E_c{Q(M)|y_0}+\lambda_n|M|)
    loss.summary = loss.summary - log(n)/2*c(rep(1,4),rep(2,6),rep(3,4),4)
    model_index = which.max(loss.summary)
    if (model_index <= 4){
      num_covariate = 1
      model_cur_loc = model_index
    }else if ((model_index >= 5)&&(model_index <= 10)){
      num_covariate = 2
      model_cur_loc = model_index - 4
    }else if ((model_index >= 11)&&(model_index <= 14)){
      num_covariate = 3
      model_cur_loc = model_index - 10
    }else{
      num_covariate = 4
      model_cur_loc = 1
    }
    beta = c(0,0,0,0)
    beta[c(CandidateModel.list[[num_covariate]][,model_cur_loc])] = c(unlist(parameter.est[[model_index]][1]))
    sigma = unlist(parameter.est[[model_index]][2])
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.01)&&(norm(beta.previous - beta,"2")<0.01)){
      break
    }else{}
    
    iter_EMS = iter_EMS + 1
  }
  
  # loss_beta = norm(beta - c(1,0,1,0),"2")
  loss_mu = norm(data_obs[,1] - data_obs[,c(2:5)]%*%beta,"2")
  
  # Record computing time
  end_time_EMS = Sys.time()
  running_time_EMS = end_time_EMS - start_time_EMS
  
  return(list(beta,sigma,min(loss.summary),model_index,loss_mu,running_time_EMS))
}


# ----------------------------- E-MA algorithm with \lambda=2 -----------------------------
EMA_A = function(data,            # data set, responses of all subjects
               mu,                # estimate of mu produced by the EM algorithm
               Pi,                # estimate of Pi produced by the EM algorithm
               Sigma.matrix,      # estimate of covariance matrix produced by the EM algorithm
               beta,              # initial estimate of the coefficient
               sigma,             # initial estimate of the variance
               CandidateModel.list# candidate model set
){
  
  # Record computing time
  start_time_EMA = Sys.time()
  
  iter_EMA = 0
  model_index = 15 
  number.candidatemodels = 15
  
  data_miss = data[-which(complete.cases(data)),]
  data_obs = data[which(complete.cases(data)),]
  
  while (iter_EMA < 100){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EMA = list()
    data_gen_EMA_matrix = matrix(0,nrow = dim(data_miss)[1]*50,ncol = 5)
    for (z in 1:dim(data_miss)[1]){
      data_gen_EMA[[z]] = MIS(data_miss[z,],mu = mu, Pi = Pi,Sigma.matrix = Sigma.matrix,
                              beta = beta,sigma = sigma,model_index = model_index)
      data_gen_EMA_matrix[c((50*z-49):(50*z)),] = matrix(unlist(data_gen_EMA[[z]]),nrow = 50,ncol = 5,byrow = FALSE)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma of different candidate models using the obtained full sample
    parameter_est_EMA = list()
    parameter_est_EMA_matrix = matrix(0,nrow = 4,ncol = 15)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 4){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 5)&&(k <= 10)){
        num_covariate = 2
        model_cur_loc = k - 4
      }else if ((k >= 11)&&(k <= 14)){
        num_covariate = 3
        model_cur_loc = k - 10
      }else{
        num_covariate = 4
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate,ncol = num_covariate)
      for (i in 1:length(data_gen_EMA)){
        data_gen_tem = matrix(unlist(data_gen_EMA[[i]]),nrow = 50, ncol = 5, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate,ncol = num_covariate)
        for (j in 1:50){
          S_1_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_gen_tem[j,1]
          S_2_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
            t(data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(n-length(data_gen_EMA))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
          t(data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
      }
      
      ## 2.3 Estimation of beta and sigma^2
      beta_est_tem = solve(S_2)%*%S_1
      beta_est_full = c(0,0,0,0)
      beta_est_full[c(CandidateModel.list[[num_covariate]][,model_cur_loc])] = beta_est_tem
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/n
      parameter_est_EMA_matrix[,k] = beta_est_full
      parameter_est_EMA[[k]] = list(beta_est_full,sigma_square_est)
    }
    
    # 3. Find the optimal model weight based on the weight choice criterion by using solve.QP function
    ## Initialize parameter setting
    sigma_square_EMA = unlist(parameter_est_EMA[[15]][2])
    A_1 = data_gen_EMA_matrix[,1]
    B_1 = data_gen_EMA_matrix[,c(2:5)]%*%parameter_est_EMA_matrix
    A_2 = data_obs[,1]
    B_2 = data_obs[,c(2:5)]%*%parameter_est_EMA_matrix
    nu = c(rep(1,4),rep(2,6),rep(3,4),4)
    
    ## Calculate the optimal model weight
    Dmat = (1/50*t(B_1)%*%B_1 + t(B_2)%*%B_2)/sigma_square_EMA + (1e-5)*diag(number.candidatemodels)
    dvec = (1/50*t(B_1)%*%A_1 + t(B_2)%*%A_2)/sigma_square_EMA - nu
    Amat = t(rbind(matrix(1,1,number.candidatemodels),diag(number.candidatemodels)))
    bvec = matrix(c(1,rep(1e-10,number.candidatemodels)))
    W_MA = solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)$solution
    
    ## Calculate the optimal model averaging estimator for beta
    beta_MA = parameter_est_EMA_matrix%*%W_MA
    
    # 4. Replace beta and sigma
    beta = beta_MA
    sigma = sigma_square_EMA
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.01)&&(norm(beta.previous - beta,"2")<0.01)){
      break
    }else{}
    
    iter_EMA = iter_EMA + 1
  }
  
  # loss_beta = norm(beta - c(1,0,1,0),"2")
  loss_mu = norm(data_obs[,1] - data_obs[,c(2:5)]%*%beta,"2")
  
  # Record computing time
  end_time_EMA = Sys.time()
  running_time_EMA = end_time_EMA - start_time_EMA
  
  return(list(beta,sigma,W_MA,loss_mu,running_time_EMA,parameter_est_EMA_matrix))
}


# ----------------------------- E-MA algorithm with \lambda=\log n -----------------------------
EMA_B = function(data,            # data set, responses of all subjects
                 mu,                # estimate of mu produced by the EM algorithm
                 Pi,                # estimate of Pi produced by the EM algorithm
                 Sigma.matrix,      # estimate of covariance matrix produced by the EM algorithm
                 beta,              # initial estimate of the coefficient
                 sigma,             # initial estimate of the variance
                 CandidateModel.list# candidate model set
){
  
  # Record computing time
  start_time_EMA = Sys.time()
  
  iter_EMA = 0
  model_index = 15 
  number.candidatemodels = 15
  
  data_miss = data[-which(complete.cases(data)),]
  data_obs = data[which(complete.cases(data)),]
  
  while (iter_EMA < 100){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EMA = list()
    data_gen_EMA_matrix = matrix(0,nrow = dim(data_miss)[1]*50,ncol = 5)
    for (z in 1:dim(data_miss)[1]){
      data_gen_EMA[[z]] = MIS(data_miss[z,],mu = mu, Pi = Pi,Sigma.matrix = Sigma.matrix,
                              beta = beta,sigma = sigma,model_index = model_index)
      data_gen_EMA_matrix[c((50*z-49):(50*z)),] = matrix(unlist(data_gen_EMA[[z]]),nrow = 50,ncol = 5,byrow = FALSE)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma of different candidate models using the obtained full sample
    parameter_est_EMA = list()
    parameter_est_EMA_matrix = matrix(0,nrow = 4,ncol = 15)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 4){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 5)&&(k <= 10)){
        num_covariate = 2
        model_cur_loc = k - 4
      }else if ((k >= 11)&&(k <= 14)){
        num_covariate = 3
        model_cur_loc = k - 10
      }else{
        num_covariate = 4
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate,ncol = num_covariate)
      for (i in 1:length(data_gen_EMA)){
        data_gen_tem = matrix(unlist(data_gen_EMA[[i]]),nrow = 50, ncol = 5, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate,ncol = num_covariate)
        for (j in 1:50){
          S_1_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_gen_tem[j,1]
          S_2_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
            t(data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(n-length(data_gen_EMA))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
          t(data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
      }
      
      ## 2.3 Estimation of beta and sigma^2
      beta_est_tem = solve(S_2)%*%S_1
      beta_est_full = c(0,0,0,0)
      beta_est_full[c(CandidateModel.list[[num_covariate]][,model_cur_loc])] = beta_est_tem
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/n
      parameter_est_EMA_matrix[,k] = beta_est_full
      parameter_est_EMA[[k]] = list(beta_est_full,sigma_square_est)
    }
    
    # 3. Find the optimal model weight based on the weight choice criterion by using solve.QP function
    ## Initialize parameter setting
    sigma_square_EMA = unlist(parameter_est_EMA[[15]][2])
    A_1 = data_gen_EMA_matrix[,1]
    B_1 = data_gen_EMA_matrix[,c(2:5)]%*%parameter_est_EMA_matrix
    A_2 = data_obs[,1]
    B_2 = data_obs[,c(2:5)]%*%parameter_est_EMA_matrix
    nu = c(rep(1,4),rep(2,6),rep(3,4),4)
    
    ## Calculate the optimal model weight
    Dmat = (1/50*t(B_1)%*%B_1 + t(B_2)%*%B_2)/sigma_square_EMA + (1e-5)*diag(number.candidatemodels)
    dvec = (1/50*t(B_1)%*%A_1 + t(B_2)%*%A_2)/sigma_square_EMA - log(n)/2*nu
    Amat = t(rbind(matrix(1,1,number.candidatemodels),diag(number.candidatemodels)))
    bvec = matrix(c(1,rep(1e-10,number.candidatemodels)))
    W_MA = solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)$solution
    
    ## Calculate the optimal model averaging estimator for beta
    beta_MA = parameter_est_EMA_matrix%*%W_MA
    
    # 4. Replace beta and sigma
    beta = beta_MA
    sigma = sigma_square_EMA
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.01)&&(norm(beta.previous - beta,"2")<0.01)){
      break
    }else{}
    
    iter_EMA = iter_EMA + 1
  }
  
  # loss_beta = norm(beta - c(1,0,1,0),"2")
  loss_mu = norm(data_obs[,1] - data_obs[,c(2:5)]%*%beta,"2")
  
  # Record computing time
  end_time_EMA = Sys.time()
  running_time_EMA = end_time_EMA - start_time_EMA
  
  return(list(beta,sigma,W_MA,loss_mu,running_time_EMA,parameter_est_EMA_matrix))
}


# ----------------------------- E-MA algorithm with equal weighting -----------------------------
EMA_EW = function(data,            # data set, responses of all subjects
                  mu,                # estimate of mu produced by the EM algorithm
                  Pi,                # estimate of Pi produced by the EM algorithm
                  Sigma.matrix,      # estimate of covariance matrix produced by the EM algorithm
                  beta,              # initial estimate of the coefficient
                  sigma,             # initial estimate of the variance
                  CandidateModel.list# candidate model set
){
  
  # Record computing time
  start_time_EMAEW = Sys.time()
  
  iter_EMA = 0
  model_index = 15 
  number.candidatemodels = 15
  
  data_miss = data[-which(complete.cases(data)),]
  data_obs = data[which(complete.cases(data)),]
  
  while (iter_EMA < 50){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EMA = list()
    data_gen_EMA_matrix = matrix(0,nrow = dim(data_miss)[1]*50,ncol = 5)
    for (z in 1:dim(data_miss)[1]){
      data_gen_EMA[[z]] = MIS(data_miss[z,],mu = mu, Pi = Pi,Sigma.matrix = Sigma.matrix,
                              beta = beta,sigma = sigma,model_index = model_index)
      data_gen_EMA_matrix[c((50*z-49):(50*z)),] = matrix(unlist(data_gen_EMA[[z]]),nrow = 50,ncol = 5,byrow = FALSE)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma of different candidate models using the obtained full sample
    parameter_est_EMA = list()
    parameter_est_EMA_matrix = matrix(0,nrow = 4,ncol = 15)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 4){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 5)&&(k <= 10)){
        num_covariate = 2
        model_cur_loc = k - 4
      }else if ((k >= 11)&&(k <= 14)){
        num_covariate = 3
        model_cur_loc = k - 10
      }else{
        num_covariate = 4
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate,ncol = num_covariate)
      for (i in 1:length(data_gen_EMA)){
        data_gen_tem = matrix(unlist(data_gen_EMA[[i]]),nrow = 50, ncol = 5, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate,ncol = num_covariate)
        for (j in 1:50){
          S_1_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_gen_tem[j,1]
          S_2_tem = data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
            t(data_gen_tem[j,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(n-length(data_gen_EMA))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1]%*%
          t(data_obs[i,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1])
      }
      
      ## 2.3 Estimation of beta and sigma^2
      beta_est_tem = solve(S_2)%*%S_1
      beta_est_full = c(0,0,0,0)
      beta_est_full[c(CandidateModel.list[[num_covariate]][,model_cur_loc])] = beta_est_tem
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/n
      parameter_est_EMA_matrix[,k] = beta_est_full
      parameter_est_EMA[[k]] = list(beta_est_full,sigma_square_est)
    }
    
    # 3. Obtain the model averaging estimator with equal weights
    ## Calculate the EW estimator for beta
    beta_MA = parameter_est_EMA_matrix%*%rep(1/number.candidatemodels,number.candidatemodels)
    
    # 4. Replace beta and sigma
    beta = beta_MA
    sigma = unlist(parameter_est_EMA[[15]][2])
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.01)&&(norm(beta.previous - beta,"2")<0.01)){
      break
    }else{}
    
    iter_EMA = iter_EMA + 1
  }
  
  # loss_beta = norm(beta - c(1,0,1,0),"2")
  loss_mu = norm(data_obs[,1] - data_obs[,c(2:5)]%*%beta,"2")
  
  # Record computing time
  end_time_EMAEW = Sys.time()
  running_time_EMAEW = end_time_EMAEW - start_time_EMAEW
  
  return(list(beta,sigma,loss_mu,running_time_EMAEW))
}


# ----------------------------- EM algorithm using the full model -----------------------------
EM_F = function(data,          # data set, responses of all subjects
                mu,            # estimate of mu produced by the EM algorithm
                Pi,            # estimate of Pi produced by the EM algorithm
                Sigma.matrix,  # estimate of covariance matrix produced by the EM algorithm
                beta,          # initial estimate of the coefficient
                sigma          # initial estimate of the variance
){
  
  # Record computing time
  start_time_EMF = Sys.time()
  
  iter_EM_F = 0
  model_index = 15 
  
  data_miss = data[-which(complete.cases(data)),]
  data_obs = data[which(complete.cases(data)),]
  
  while (iter_EM_F < 50){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EM = list()
    for (z in 1:dim(data_miss)[1]){
      data_gen_EM[[z]] = MIS(data_miss[z,],mu = mu, Pi = Pi, Sigma.matrix = Sigma.matrix,
                             beta = beta, sigma = sigma, model_index = model_index)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_f{Q(M)|y_0} and estimate paramters beta and sigma under the full model
    ## 2.1 Compute (1) S_0 = \sum_{i=1}^n E_f(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_f(X_{i,c}Y_i|y_0,x_0)
    ##             (3) S_2 = \sum_{i=1}^n E_f(X_{i,c}X_{i,c}^\prime|y_0,x_0)
    S_0 = 0
    num_covariate = 4
    S_1 = matrix(0,nrow = num_covariate,ncol = 1)
    S_2 = matrix(0,nrow = num_covariate,ncol = num_covariate)
    for (i in 1:length(data_gen_EM)){
      data_gen_tem = matrix(unlist(data_gen_EM[[i]]),nrow = 50, ncol = 5, byrow = FALSE)
      S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
      S_1_tem = matrix(0,nrow = num_covariate,ncol = 1)
      S_2_tem = matrix(0,nrow = num_covariate,ncol = num_covariate)
      for (j in 1:50){
        S_1_tem = data_gen_tem[j,c(2:5)]*data_gen_tem[j,1]
        S_2_tem = data_gen_tem[j,c(2:5)]%*%t(data_gen_tem[j,c(2:5)])
      }
      S_1 = S_1 + S_1_tem/50
      S_2 = S_2 + S_2_tem/50
    }
    for (i in 1:(n-length(data_gen_EM))){
      S_0 = S_0 + data_obs[i,1]^2
      S_1 = S_1 + data_obs[i,c(2:5)]*data_obs[i,1]
      S_2 = S_2 + data_obs[i,c(2:5)]%*%t(data_obs[i,c(2:5)])
    }
    
    ## 2.2 Estimation of beta and sigma^2
    beta_est_full = solve(S_2)%*%S_1
    sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/n
    
    # 3. Replace beta and sigma
    beta = beta_est_full
    sigma = sigma_square_est
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.01)&&(norm(beta.previous - beta,"2")<0.01)){
      break
    }else{}
    
    iter_EM_F = iter_EM_F + 1
  }
  
  # loss_beta = norm(beta - c(1,0,1,0),"2")
  loss_mu = norm(data_obs[,1] - data_obs[,c(2:5)]%*%beta,"2")
  
  # Record computing time
  end_time_EMF = Sys.time()
  running_time_EMF = end_time_EMF - start_time_EMF
  
  return(list(beta,sigma,loss_mu,running_time_EMF))
}


# ----------------------------- BIC method with fully complete data -----------------------------
BIC_CC = function(response,              # Fully observed responses
                  covariate,             # Fully observed covariates
                  missing.data,          # Missing data
                  CandidateModel.list    # Candidate model set
){
  
  data_obs = missing.data[which(complete.cases(missing.data)),]
  number.candidatemodels = 15
  
  beta.list = list()
  loss.summary = matrix(0,nrow = number.candidatemodels,ncol = 1)
  
  for (k in 1:number.candidatemodels){
    if (k <= 4){
      num_covariate = 1
      model_cur_loc = k
    }else if ((k >= 5)&&(k <= 10)){
      num_covariate = 2
      model_cur_loc = k - 4
    }else if ((k >= 11)&&(k <= 14)){
      num_covariate = 3
      model_cur_loc = k - 10
    }else{
      num_covariate = 4
      model_cur_loc = 1
    }
    
    X_tem = covariate[,CandidateModel.list[[num_covariate]][,model_cur_loc]]
    model = glm(response~-1 + X_tem,family = gaussian,maxit = 200)
    loss.summary[k,1] = BIC(model)
    beta_tem = c(0,0,0,0)
    beta_tem[CandidateModel.list[[num_covariate]][,model_cur_loc]] = coef(model)
    beta.list[[k]] = beta_tem
  }
  
  # Selection with BIC
  model_optimal = which.min(loss.summary)
  beta_optimal = unlist(beta.list[[model_optimal]])
  # loss_beta = norm(beta_optimal - c(1,0,1,0),"2")
  loss_mu = norm(data_obs[,1] - data_obs[,c(2:5)]%*%beta_optimal,"2")
  
  return(list(beta_optimal,model_optimal,loss_mu))
}


# ----------------------------- BIC method with observed complete data -----------------------------
BIC_OC = function(missing.data,          # Missing data
                  CandidateModel.list    # Candidate model set
){
  
  data_obs = missing.data[which(complete.cases(missing.data)),]
  number.candidatemodels = 15
  
  beta.list = list()
  loss.summary = matrix(0,nrow = number.candidatemodels,ncol = 1)
  
  for (k in 1:number.candidatemodels){
    if (k <= 4){
      num_covariate = 1
      model_cur_loc = k
    }else if ((k >= 5)&&(k <= 10)){
      num_covariate = 2
      model_cur_loc = k - 4
    }else if ((k >= 11)&&(k <= 14)){
      num_covariate = 3
      model_cur_loc = k - 10
    }else{
      num_covariate = 4
      model_cur_loc = 1
    }
    
    X_tem = data_obs[,CandidateModel.list[[num_covariate]][,model_cur_loc] + 1]
    Y = data_obs[,1]
    model = glm(Y~-1 + X_tem,family = gaussian,maxit = 200)
    loss.summary[k,1] = BIC(model)
    beta_tem = c(0,0,0,0)
    beta_tem[CandidateModel.list[[num_covariate]][,model_cur_loc]] = coef(model)
    beta.list[[k]] = beta_tem
  }
  
  # Selection with BIC
  model_optimal = which.min(loss.summary)
  beta_optimal = unlist(beta.list[[model_optimal]])
  # loss_beta = norm(beta_optimal - c(1,0,1,0),"2")
  loss_mu = norm(data_obs[,1] - data_obs[,c(2:5)]%*%beta_optimal,"2")
  
  return(list(beta_optimal,model_optimal,loss_mu))
}


# ------------------------ Model selection and averaging with imputed data ------------------------
MS_Im = function(missing.data,          # Missing data
                 CandidateModel.list    # Candidate model set
){
  
  data_obs = missing.data[which(complete.cases(missing.data)),]
  number.candidatemodels = 15
  
  beta.selection = matrix(0, nrow = 4, ncol = 1)
  beta.averaging = matrix(0, nrow = 4, ncol = 1)
  
  # Implement imputation to missing data
  initial = mice(missing.data, maxit = 0)
  imputed.data.tem = as.matrix(complete(initial, action = 1))
  
  selection.cri = matrix(0,nrow = 1, ncol = number.candidatemodels)
  averaging.cri = matrix(0,nrow = 1, ncol = number.candidatemodels)
  beta.selection.tem = matrix(0,nrow = 4, ncol = number.candidatemodels)
  
  for (k in 1:number.candidatemodels){
    if (k <= 4){
      num_covariate = 1
      model_cur_loc = k
    }else if ((k >= 5)&&(k <= 10)){
      num_covariate = 2
      model_cur_loc = k - 4
    }else if ((k >= 11)&&(k <= 14)){
      num_covariate = 3
      model_cur_loc = k - 10
    }else{
      num_covariate = 4
      model_cur_loc = 1
    }
    
    X_tem = imputed.data.tem[,c(CandidateModel.list[[num_covariate]][,model_cur_loc] + 1)]
    Y = imputed.data.tem[,1]
    model = glm(Y~-1 + X_tem,family = gaussian, maxit = 200)
    selection.cri[1,k] = BIC(model)
    beta_tem = rep(0,4)
    beta_tem[CandidateModel.list[[num_covariate]][,model_cur_loc]] = coef(model)
    beta.selection.tem[,k] = beta_tem
    
    # Averaging with SBIC
    SBIC = t(matrix(rep(-selection.cri,number.candidatemodels), ncol = number.candidatemodels, byrow = TRUE)
             + c(selection.cri))
    w.SBIC = 1/colSums(apply(SBIC/2, 2,exp)) 
    
    # Selection with BIC
    model_optimal = which.min(selection.cri)
    beta.selection[,1] = beta.selection.tem[,model_optimal]
    beta.averaging[,1] = beta.selection.tem%*%(w.SBIC)
  }
  
  beta_optimal_selection = rowSums(beta.selection)
  beta_optimal_averaging = rowSums(beta.averaging)
  
  return(list(beta_optimal_selection, beta_optimal_averaging))
}


# ----------------------- Main function -----------------------------------
Solve_fun = function(n,rep_time,r){
  
  # ------------------ Complete data generating process -------------------------
  set.seed(1 + r)
  p_c = 2
  p_d = 2
  Sigma = matrix(c(1,0,0,1),nrow = 2, ncol = 2)
  X = matrix(0,nrow = n,ncol = p_c + p_d)
  number.random = runif(n, min = 0, max = 1)
  eps = rnorm(n, mean = 0, sd = sqrt(3.25))
  mu.sum = matrix(c(0,0,1,1,0,1,0,1),nrow = 2, ncol = 4,byrow = TRUE)
  for (i in 1:n){
    if (number.random[i]>=0&&number.random[i]<0.25){
      X[i,c(3,4)] = c(0,0)
      X[i,c(1,2)] = mvrnorm(1, mu.sum[,1], Sigma)
    }else if (number.random[i]>=0.25&&number.random[i]<0.5){
      X[i,c(3,4)] = c(0,1)
      X[i,c(1,2)] = mvrnorm(1, mu.sum[,2], Sigma)
    }else if (number.random[i]>=0.5&&number.random[i]<0.75){
      X[i,c(3,4)] = c(1,0)
      X[i,c(1,2)] = mvrnorm(1, mu.sum[,3], Sigma)
    }else if (number.random[i]>=0.75&&number.random[i]<=1){
      X[i,c(3,4)] = c(1,1)
      X[i,c(1,2)] = mvrnorm(1, mu.sum[,4], Sigma)
    }
  }
  beta_true = c(1,0,1,0)
  Y = X%*%beta_true + 2*X[,2]^2 + eps
  
  
  # ------------------------- Generate missing data ------------------------------
  prob_m = 0.1
  miss = matrix(rbinom(n*5,size = 1,prob = prob_m),nrow = n,ncol = 5)
  miss.index = which(miss!=0) # Find the indexes of non-zero values of vector miss
  data.tem = as.vector(cbind(Y,X))
  data.tem[miss.index] = NA
  data = matrix(data.tem,nrow = n,ncol = 5,byrow = FALSE)
  # Delete those samples with more than three missing variables
  miss.index.four = c(0)
  for (i in 1:n){
    if (sum(is.na(data[i,]))>3){
      miss.index.four = append(miss.index.four,i)
    }
  }
  if (length(miss.index.four) == 1){
    data = data
  }else{
    miss.index.four = miss.index.four[-1]
    data[miss.index.four,1] = Y[miss.index.four,1]
    data[miss.index.four,c(2:5)] = X[miss.index.four,]
  }
  
  
  # -------------------- Initial setting for the EMS algorithm -------------------
  data.obs = data[which(complete.cases(data)),]
  num.complete = nrow(data.obs)
  # Individuals that have both complete X and Y
  X.obs = data.obs[,c(2:5)]
  Y.obs = data.obs[,c(1)]
  
  
  # Initial estimates for parameters \beta, sigma_square, \mu, \Pi and Sigma matrix
  beta.ini = solve(t(X.obs)%*%X.obs)%*%t(X.obs)%*%Y.obs
  sigma.square.ini = norm(Y.obs - X.obs%*%beta.ini,type = "F")/(num.complete - 4)
  num.X.dis.complete = nrow(data[complete.cases(data[,c(4,5)]),])
  mu.ini = matrix(0,nrow = 2, ncol = 4)
  Pi.ini = matrix(0,nrow = 4, ncol = 1)
  n.ro = matrix(0,nrow = 4,ncol = 1)
  X.dis.obs.index = list()
  I.o = list()
  Sigma.matrix.ini = matrix(c(0,0,0,0), nrow = 2, ncol = 2)
  for (i in 1:4){
    X.dis.obs.index[[i]] = intersect(which(data[,c(4)] == mu.sum[1,i]),which(data[,c(5)] == mu.sum[2,i]))
    Pi.ini[i,1] = length(X.dis.obs.index[[i]])/num.X.dis.complete
    I.o[[i]] = intersect(which(complete.cases(data[,c(2:5)])),X.dis.obs.index[[i]])
    n.ro[i,1] = length(I.o[[i]])
    mu.ini[,i] = colSums(data[I.o[[i]],c(2,3)])/n.ro[i,1]
    for (j in 1:length(I.o[[i]])){
      Sigma.matrix.ini = Sigma.matrix.ini + (data[I.o[[i]][j],c(2,3)] - mu.ini[,i])%*%t((data[I.o[[i]][j],c(2,3)] - mu.ini[,i]))
    }
  }
  Sigma.matrix.ini = Sigma.matrix.ini/sum(n.ro)
  
  
  EM.result = EM(data, mu = mu.ini, Pi = Pi.ini, Sigma.matrix = Sigma.matrix.ini)
  
  
  # Build non-nested candidate models
  CandidateModel.index.list = list()
  number.candidatemodels = 0
  p = p_c + p_d
  for(i in 1:p){
    CandidateModel.index.list[[i]] = combn(1:p, i)
    number.candidatemodels = number.candidatemodels + ncol(CandidateModel.index.list[[i]])
  }
  
  
  MS.Im.result = MS_Im(missing.data = data, CandidateModel.list = CandidateModel.index.list)
  
  
  EMS.result = EMS(data = data, mu = EM.result[[1]],Pi = EM.result[[2]],Sigma.matrix = EM.result[[3]],
                   beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)
  
  
  EMA_A.result = EMA_A(data = data, mu = EM.result[[1]],Pi = EM.result[[2]],Sigma.matrix = EM.result[[3]],
                   beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)
  
  
  EMA_B.result = EMA_B(data = data, mu = EM.result[[1]],Pi = EM.result[[2]],Sigma.matrix = EM.result[[3]],
                       beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)
  
  
  EMA_EW.result = EMA_EW(data = data, mu = EM.result[[1]],Pi = EM.result[[2]],Sigma.matrix = EM.result[[3]],
                         beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)


  EM_F.result = EM_F(data = data, mu = EM.result[[1]],Pi = EM.result[[2]],Sigma.matrix = EM.result[[3]],
                     beta = beta.ini,sigma = sigma.square.ini)


  BIC.CC.result = BIC_CC(response = Y, covariate = X, missing.data = data, CandidateModel.list = CandidateModel.index.list)


  BIC.OC.result = BIC_OC(missing.data = data, CandidateModel.list = CandidateModel.index.list)
  
  
  address2 = paste0(address,"results/","Misspecified setting/","tempdata/")
  if(!dir.exists(address2)){
    dir.create(paste0(address,"results/"))
    dir.create(paste0(address,"results/","Misspecified setting/"))
    dir.create(address2)
  }
  
  
  
  # Generate fully-observed test data
  n_test = 100
  X_test = matrix(0,nrow = n_test,ncol = p_c + p_d)
  number.random = runif(n_test, min = 0, max = 1)
  eps_test = rnorm(n_test, mean = 0, sd = sqrt(3.25))
  mu.sum = matrix(c(0,0,1,1,0,1,0,1),nrow = 2, ncol = 4,byrow = TRUE)
  for (i in 1:n_test){
    if (number.random[i]>=0&&number.random[i]<0.25){
      X_test[i,c(3,4)] = c(0,0)
      X_test[i,c(1,2)] = mvrnorm(1, mu.sum[,1], Sigma)
    }else if (number.random[i]>=0.25&&number.random[i]<0.5){
      X_test[i,c(3,4)] = c(0,1)
      X_test[i,c(1,2)] = mvrnorm(1, mu.sum[,2], Sigma)
    }else if (number.random[i]>=0.5&&number.random[i]<0.75){
      X_test[i,c(3,4)] = c(1,0)
      X_test[i,c(1,2)] = mvrnorm(1, mu.sum[,3], Sigma)
    }else if (number.random[i]>=0.75&&number.random[i]<=1){
      X_test[i,c(3,4)] = c(1,1)
      X_test[i,c(1,2)] = mvrnorm(1, mu.sum[,4], Sigma)
    }
  }
  Y_test = X_test%*%beta_true + 2*X_test[,2]^2 + eps_test
  
  # Comparison of KL divergences between EMA and EMS
  KL.EMS = norm(Y_test - X_test%*%c(unlist(EMS.result[[1]])),"2")/(2*unlist(EMS.result[[2]])) + log(sqrt(2*pi*unlist(EMS.result[[2]])))
  KL.EMA_A = norm(Y_test - X_test%*%c(unlist(EMA_A.result[[1]])),"2")/(2*unlist(EMA_A.result[[2]])) + log(sqrt(2*pi*unlist(EMA_A.result[[2]])))
  KL.EMA_B = norm(Y_test - X_test%*%c(unlist(EMA_B.result[[1]])),"2")/(2*unlist(EMA_B.result[[2]])) + log(sqrt(2*pi*unlist(EMA_B.result[[2]])))
  KL.EMA.EW = norm(Y_test - X_test%*%c(unlist(EMA_EW.result[[1]])),"2")/(2*unlist(EMA_EW.result[[2]])) + log(sqrt(2*pi*unlist(EMA_EW.result[[2]])))
  KL.EM_F = norm(Y_test - X_test%*%c(unlist(EM_F.result[[1]])),"2")/(2*unlist(EM_F.result[[2]])) + log(sqrt(2*pi*unlist(EM_F.result[[2]])))
  
  
  # Verification of asymptotic optimality
  Dmat_inf_A = 2*t(X_test%*%unlist(EMA_A.result[[6]]))%*%(X_test%*%unlist(EMA_A.result[[6]])) + (1e-5)*diag(number.candidatemodels)
  dvec_inf_A = 2*t(X_test%*%unlist(EMA_A.result[[6]]))%*%(Y_test)
  Amat_inf = t(rbind(matrix(1,1,number.candidatemodels),diag(number.candidatemodels)))
  bvec_inf = matrix(c(1,rep(1e-10,number.candidatemodels)))
  W.inf_A = solve.QP(Dmat = Dmat_inf_A,dvec = dvec_inf_A,Amat = Amat_inf,bvec = bvec_inf,meq=1)$solution
  KL.EMA.inf_A = norm(Y_test - X_test%*%unlist(EMA_A.result[[6]])%*%W.inf_A,"2")/(2*unlist(EMA_A.result[[2]])) + log(sqrt(2*pi*unlist(EMA_A.result[[2]])))
  ratio.loss_A = KL.EMA_A/KL.EMA.inf_A

  Dmat_inf_B = 2*t(X_test%*%unlist(EMA_B.result[[6]]))%*%(X_test%*%unlist(EMA_B.result[[6]])) + (1e-5)*diag(number.candidatemodels)
  dvec_inf_B = 2*t(X_test%*%unlist(EMA_B.result[[6]]))%*%(Y_test)
  W.inf_B = solve.QP(Dmat = Dmat_inf_B,dvec = dvec_inf_B,Amat = Amat_inf,bvec = bvec_inf,meq=1)$solution
  KL.EMA.inf_B = norm(Y_test - X_test%*%unlist(EMA_B.result[[6]])%*%W.inf_B,"2")/(2*unlist(EMA_B.result[[2]])) + log(sqrt(2*pi*unlist(EMA_B.result[[2]])))
  ratio.loss_B = KL.EMA_B/KL.EMA.inf_B
  

  # Calculate the loss function of mu on the test dataset
  loss.mu.EMS = norm(Y_test - X_test%*%c(unlist(EMS.result[[1]])),"2")/n_test
  loss.mu.EMA_A = norm(Y_test - X_test%*%c(unlist(EMA_A.result[[1]])),"2")/n_test
  loss.mu.EMA_B = norm(Y_test - X_test%*%c(unlist(EMA_B.result[[1]])),"2")/n_test
  loss.mu.EMA.EW = norm(Y_test - X_test%*%c(unlist(EMA_EW.result[[1]])),"2")/n_test
  loss.mu.EM.F = norm(Y_test - X_test%*%c(unlist(EM_F.result[[1]])),"2")/n_test
  loss.mu.BIC.CC = norm(Y_test - X_test%*%c(unlist(BIC.CC.result[[1]])),"2")/n_test
  loss.mu.BIC.OC = norm(Y_test - X_test%*%c(unlist(BIC.OC.result[[1]])),"2")/n_test
  loss.mu.MS.Im.SE = norm(Y_test - X_test%*%c(unlist(MS.Im.result[[1]])),"2")/n_test 
  loss.mu.MS.Im.MA = norm(Y_test - X_test%*%c(unlist(MS.Im.result[[2]])),"2")/n_test 
  loss.mu.summary = c(loss.mu.EMS,loss.mu.EMA_A,loss.mu.EMA_B,loss.mu.EMA.EW,loss.mu.EM.F,loss.mu.BIC.CC,loss.mu.BIC.OC,
                      loss.mu.MS.Im.SE,loss.mu.MS.Im.MA)
  # loss.mu.summary = c(KL.EMS,KL.EMA_A,KL.EMA_B,KL.EMA.EW,KL.EM_F,ratio.loss_A,ratio.loss_B)
  resultframe = data.frame(t(c(n,loss.mu.summary)))
  
  write.csv(resultframe, paste0(address2,"n=",n,"_rept=",r,".csv"),row.names = FALSE)
  return(resultframe)
}

