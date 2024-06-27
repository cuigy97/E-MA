#---------------------------------------------------------------------
# Title: The EMA Algorithm - Real data analysis on air quality dataset
#---------------------------------------------------------------------

# -------- Remove all variables ---------
if(sys.nframe() == 0L){rm(list=ls()); gc()}

# address = ''
# setwd('')

# ---- Necessary packages ----
pacman::p_load(expm,mvtnorm,MASS,quadprog,mice,snowfall,parallel,ggplot2)

# ---- Functions required ----
source("MIS_covariate_RD.R")
source("MIS_RD.R")

#------ Function that adjusts the order of the means and covariance matrix based on the index of missing values--------

cov_split = function(mu,              # mean vector of covariates
                     covariance,      # covariance matrix of covariates
                     missing.index    # the index set of missing covariates
){
  
  obs_index = rep(1:6)[-missing.index]
  index_re = c(missing.index,obs_index)
  mu_re = as.matrix(mu[index_re])
  covariance_re = covariance[index_re,index_re]
  
  return(list(mu_re,covariance_re))
  
}


# -------- Perform the EM algorithm to obtain estimates of \mu, \Sigma and \Pi --------

EM = function(data,          # dataset, including covariates and responses
              mu,            # initial value of mu (conditional mean vector of continuous variables given the value of discrete variables)
              Sigma.matrix   # initial value of Sigma.matrix (conditional covariance matrix of continuous variables given the value of discrete ones)
){
  
  # Initial parameter setting
  data_EM = data[-which(complete.cases(data[,-1]) == TRUE),]
  data_ori = data[which(complete.cases(data[,-1]) == TRUE),-1]
  mu.current = mu
  Sigma.matrix.current = Sigma.matrix
  iter = 0
  
  
  while (iter < 100){
    data_gen = list()
    for (z in 1:dim(data_EM)[1]){
      data_gen[[z]] = MIS_RD(data_EM[z,],mu = mu.current, Sigma.matrix = Sigma.matrix.current,
                             beta = rep(1000,6), sigma = 1000, model_index = 31)
    }
    
    mu.previous = mu.current
    Sigma.matrix.previous = Sigma.matrix.current
    
    # Calculate the estimates of mu by the EM algorithm
    mu.current = rep(0,6)
    for (i in 1:length(data_gen)){
      data_tem = matrix(unlist(data_gen[[i]]),nrow = 50,ncol = 7, byrow = FALSE)
      mu.current = mu.current + colSums(data_tem[,-1])/50
    }
    mu.current = (mu.current + colSums(data_ori))/dim(data)[1]
    
    # Calculate the estimates of Omega by the EM algorithm
    Sigma.matrix.current = matrix(0,ncol = 6, nrow = 6)
    # for (i in 1:length(data_gen)){
    for (i in 1:3){
      data_tem = matrix(unlist(data_gen[[i]]),nrow = 50,ncol = 7, byrow = FALSE)
      Sigma.matrix.current_tem = matrix(0,ncol = 6, nrow = 6)
      for (j in 1:50){
        Sigma.matrix.current_tem = Sigma.matrix.current_tem + (data_tem[j,-1] - mu.current)%*%t(data_tem[j,-1] - mu.current)
      }
      Sigma.matrix.current = Sigma.matrix.current + Sigma.matrix.current_tem/50
    }
    for (k in 1:(dim(data)[1] - length(data_gen))){
      Sigma.matrix.current = Sigma.matrix.current + (data_ori[k,] - mu.current)%*%t(data_ori[k,] - mu.current)
    }
    Sigma.matrix.current = Sigma.matrix.current/dim(data)[1]
    
    
    # Stop criterion
    if ((norm(mu.previous - mu.current,"2")<0.001)&&(norm(Sigma.matrix.previous - Sigma.matrix.current,"2")<0.001)){
      break
    }else{}
    
    iter = iter + 1
  }
  return(list(mu.current,Sigma.matrix.current))
}


# ----------------------------- E-MS function -----------------------------

EMS = function(data,          # data set, responses of all subjects
               mu,            # estimate of mu produced by the EM algorithm
               Sigma.matrix,  # estimate of covariance matrix produced by the EM algorithm
               beta,          # initial estimate of the coefficient
               sigma,         # initial estimate of the variance
               CandidateModel.list
){
  
  # Record computing time
  start_time_EMS = Sys.time()
  
  iter_EMS = 0
  model_index = 31 # Assume the initial model for the EMS algorithm is the full model
  number.candidatemodels = 31
  
  data_miss = data[-which(complete.cases(data) == TRUE),]
  data_obs = data[which(complete.cases(data) == TRUE),]
  
  while (iter_EMS < 50){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current model M_c
    data_gen_EMS = apply(data_miss, MARGIN = 1, MIS_RD, mu = mu,
                         Sigma.matrix = Sigma.matrix,beta = beta,sigma = sigma, model_index = model_index)
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma under different candidate models using the obtained full sample
    parameter.est = list()
    loss.summary = matrix(0,nrow = number.candidatemodels,ncol = 1)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 5){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 6)&&(k <= 15)){
        num_covariate = 2
        model_cur_loc = k - 5
      }else if ((k >= 16)&&(k <= 25)){
        num_covariate = 3
        model_cur_loc = k - 15
      }else if ((k >= 26)&&(k <= 30)){
        num_covariate = 4
        model_cur_loc = k - 25
      }else{
        num_covariate = 5
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate + 1,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
      for (i in 1:length(data_gen_EMS)){
        data_gen_tem = matrix(unlist(data_gen_EMS[[i]]),nrow = 50, ncol = 7, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate + 1,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
        for (j in 1:50){
          S_1_tem = S_1_tem + data_gen_tem[j,c(2, c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_gen_tem[j,1]
          S_2_tem = S_2_tem + data_gen_tem[j,c(2, c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
            t(data_gen_tem[j,c(2, c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(dim(data)[1] - length(data_gen_EMS))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(2, c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(2, c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
          t(data_obs[i,c(2, c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
      }
      
      ## 2.3 Estimation of beta, sigma^2 and loss function
      beta_est = solve(S_2)%*%S_1
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/dim(data)[1]
      
      RSS_est = 0
      for (i in 1:length(data_gen_EMS)){
        data_gen_tem = matrix(unlist(data_gen_EMS[[i]]),nrow = 50, ncol = 7, byrow = FALSE)
        RSS_est_tem = 0
        for (j in 1:50){
          RSS_est_tem = RSS_est_tem + (data_gen_tem[j,1] - t(data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc])+1)])%*%beta_est)^2
        }
        RSS_est = RSS_est + RSS_est_tem/50
      }
      for (i in 1:(dim(data)[1] - length(data_gen_EMS))){
        RSS_est = RSS_est + (data_obs[i,1] - t(data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])%*%beta_est)^2
      }
      loss_est = -dim(data)[1]*log(sigma_square_est)/2 - 1/(2*sigma_square_est)*RSS_est
      parameter.est[[k]] = list(beta_est,sigma_square_est)
      loss.summary[k,1] = loss_est
    }
    
    # 3. Find the optimal model based on (E_c{Q(M)|y_0}+\lambda_n|M|)
    loss.summary = loss.summary - log(dim(data)[1])/2*c(rep(2,5),rep(3,10),rep(4,10),rep(5,5),6)
    model_index = which.max(loss.summary)
    if (model_index <= 5){
      num_covariate = 1
      model_cur_loc = model_index
    }else if ((model_index >= 6)&&(model_index <= 15)){
      num_covariate = 2
      model_cur_loc = model_index - 5
    }else if ((model_index >= 16)&&(model_index <= 25)){
      num_covariate = 3
      model_cur_loc = model_index - 15
    }else if ((model_index >= 26)&&(model_index <= 30)){
      num_covariate = 4
      model_cur_loc = model_index - 25
    }else{
      num_covariate = 5
      model_cur_loc = 1
    }
    beta = rep(0,6)
    beta[c(1,CandidateModel.list[[num_covariate]][,model_cur_loc])] = c(unlist(parameter.est[[model_index]][1]))
    sigma = unlist(parameter.est[[model_index]][2])
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.001)&&(norm(beta.previous - beta,"2")<0.001)){
      break
    }else{}
    
    iter_EMS = iter_EMS + 1
  }
  
  
  # Record computing time
  end_time_EMS = Sys.time()
  running_time_EMS = end_time_EMS - start_time_EMS
  
  return(list(beta,sigma,min(loss.summary),model_index,running_time_EMS))
}


# ----------------------------- E-MA_A function -----------------------------
EMA_A = function(data,          # data set, responses of all subjects
                 mu,            # estimate of mu produced by the EM algorithm
                 Sigma.matrix,  # estimate of covariance matrix produced by the EM algorithm
                 beta,          # initial estimate of the coefficient
                 sigma,         # initial estimate of the variance
                 CandidateModel.list
){
  
  # Record computing time
  start_time_EMA = Sys.time()
  
  iter_EMA = 0
  model_index = 31
  number.candidatemodels = 31
  
  data_miss = data[-which(complete.cases(data) == TRUE),]
  data_obs = data[which(complete.cases(data) == TRUE),]
  
  while (iter_EMA < 100){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EMA = list()
    data_gen_EMA_matrix = matrix(0,nrow = dim(data_miss)[1]*50,ncol = 7)
    for (z in 1:dim(data_miss)[1]){
      data_gen_EMA[[z]] = MIS_RD(data_miss[z,],mu = mu,Sigma.matrix = Sigma.matrix,
                                 beta = beta, sigma = sigma, model_index = model_index)
      data_gen_EMA_matrix[c((50*z-49):(50*z)),] = matrix(unlist(data_gen_EMA[[z]]),nrow = 50,ncol = 7,byrow = FALSE)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma of different candidate models using the obtained full sample
    parameter_est_EMA = list()
    parameter_est_EMA_matrix = matrix(0,nrow = 6,ncol = 31)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 5){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 6)&&(k <= 15)){
        num_covariate = 2
        model_cur_loc = k - 5
      }else if ((k >= 16)&&(k <= 25)){
        num_covariate = 3
        model_cur_loc = k - 15
      }else if ((k >= 26)&&(k <= 30)){
        num_covariate = 4
        model_cur_loc = k - 25
      }else{
        num_covariate = 5
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate + 1,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
      for (i in 1:length(data_gen_EMA)){
        data_gen_tem = matrix(unlist(data_gen_EMA[[i]]),nrow = 50, ncol = 7, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate + 1,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
        for (j in 1:50){
          S_1_tem = S_1_tem + data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_gen_tem[j,1]
          S_2_tem = S_2_tem + data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
            t(data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(dim(data)[1] - length(data_gen_EMA))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
          t(data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
      }
      
      ## 2.3 Estimation of beta and sigma^2
      beta_est_tem = solve(S_2)%*%S_1
      beta_est_full = rep(0,6)
      beta_est_full[c(1,CandidateModel.list[[num_covariate]][,model_cur_loc])] = beta_est_tem
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/dim(data)[1]
      parameter_est_EMA_matrix[,k] = beta_est_full
      parameter_est_EMA[[k]] = list(beta_est_full,sigma_square_est)
    }
    
    # 3. Find the optimal model weight based on the weight choice criterion by using solve.QP function
    ## Initialize parameter setting
    sigma_square_EMA = unlist(parameter_est_EMA[[31]][2])
    A_1 = data_gen_EMA_matrix[,1]
    B_1 = data_gen_EMA_matrix[,c(2:7)]%*%parameter_est_EMA_matrix
    A_2 = data_obs[,1]
    B_2 = data_obs[,c(2:7)]%*%parameter_est_EMA_matrix
    nu = c(rep(2,5),rep(3,10),rep(4,10),rep(5,5),6)
    
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
    if ((abs(sigma - sigma.previous)<0.001)&&(norm(beta.previous - beta,"2")<0.001)){
      break
    }else{}
    
    iter_EMA = iter_EMA + 1
  }
  
  # Record computing time
  end_time_EMA = Sys.time()
  running_time_EMA = end_time_EMA - start_time_EMA
  
  return(list(beta,sigma,W_MA,running_time_EMA))
}


# ----------------------------- E-MA_B function -----------------------------
EMA_B = function(data,          # data set, responses of all subjects
                 mu,            # estimate of mu produced by the EM algorithm
                 Sigma.matrix,  # estimate of covariance matrix produced by the EM algorithm
                 beta,          # initial estimate of the coefficient
                 sigma,         # initial estimate of the variance
                 CandidateModel.list
){
  
  # Record computing time
  start_time_EMA = Sys.time()
  
  iter_EMA = 0
  model_index = 31
  number.candidatemodels = 31
  
  data_miss = data[-which(complete.cases(data) == TRUE),]
  data_obs = data[which(complete.cases(data) == TRUE),]
  
  while (iter_EMA < 100){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EMA = list()
    data_gen_EMA_matrix = matrix(0,nrow = dim(data_miss)[1]*50,ncol = 7)
    for (z in 1:dim(data_miss)[1]){
      data_gen_EMA[[z]] = MIS_RD(data_miss[z,],mu = mu,Sigma.matrix = Sigma.matrix,
                                 beta = beta, sigma = sigma, model_index = model_index)
      data_gen_EMA_matrix[c((50*z-49):(50*z)),] = matrix(unlist(data_gen_EMA[[z]]),nrow = 50,ncol = 7,byrow = FALSE)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma of different candidate models using the obtained full sample
    parameter_est_EMA = list()
    parameter_est_EMA_matrix = matrix(0,nrow = 6,ncol = 31)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 5){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 6)&&(k <= 15)){
        num_covariate = 2
        model_cur_loc = k - 5
      }else if ((k >= 16)&&(k <= 25)){
        num_covariate = 3
        model_cur_loc = k - 15
      }else if ((k >= 26)&&(k <= 30)){
        num_covariate = 4
        model_cur_loc = k - 25
      }else{
        num_covariate = 5
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate + 1,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
      for (i in 1:length(data_gen_EMA)){
        data_gen_tem = matrix(unlist(data_gen_EMA[[i]]),nrow = 50, ncol = 7, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate + 1,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
        for (j in 1:50){
          S_1_tem = S_1_tem + data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_gen_tem[j,1]
          S_2_tem = S_2_tem + data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
            t(data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(dim(data)[1] - length(data_gen_EMA))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
          t(data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
      }
      
      ## 2.3 Estimation of beta and sigma^2
      beta_est_tem = solve(S_2)%*%S_1
      beta_est_full = rep(0,6)
      beta_est_full[c(1,CandidateModel.list[[num_covariate]][,model_cur_loc])] = beta_est_tem
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/dim(data)[1]
      parameter_est_EMA_matrix[,k] = beta_est_full
      parameter_est_EMA[[k]] = list(beta_est_full,sigma_square_est)
    }
    
    # 3. Find the optimal model weight based on the weight choice criterion by using solve.QP function
    ## Initialize parameter setting
    sigma_square_EMA = unlist(parameter_est_EMA[[31]][2])
    A_1 = data_gen_EMA_matrix[,1]
    B_1 = data_gen_EMA_matrix[,c(2:7)]%*%parameter_est_EMA_matrix
    A_2 = data_obs[,1]
    B_2 = data_obs[,c(2:7)]%*%parameter_est_EMA_matrix
    nu = c(rep(2,5),rep(3,10),rep(4,10),rep(5,5),6)
    
    ## Calculate the optimal model weight
    Dmat = (1/50*t(B_1)%*%B_1 + t(B_2)%*%B_2)/sigma_square_EMA + (1e-5)*diag(number.candidatemodels)
    dvec = (1/50*t(B_1)%*%A_1 + t(B_2)%*%A_2)/sigma_square_EMA - log(dim(data)[1])/2*nu
    Amat = t(rbind(matrix(1,1,number.candidatemodels),diag(number.candidatemodels)))
    bvec = matrix(c(1,rep(1e-10,number.candidatemodels)))
    W_MA = solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)$solution
    
    ## Calculate the optimal model averaging estimator for beta
    beta_MA = parameter_est_EMA_matrix%*%W_MA
    
    # 4. Replace beta and sigma
    beta = beta_MA
    sigma = sigma_square_EMA
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.001)&&(norm(beta.previous - beta,"2")<0.001)){
      break
    }else{}
    
    iter_EMA = iter_EMA + 1
  }
  
  # Record computing time
  end_time_EMA = Sys.time()
  running_time_EMA = end_time_EMA - start_time_EMA
  
  return(list(beta,sigma,W_MA,running_time_EMA))
}

# ----------------------------- E-MA function with equal weighting -----------------------------
EMA_EW = function(data,          # data set, responses of all subjects
                  mu,            # estimate of mu produced by the EM algorithm
                  Sigma.matrix,  # estimate of covariance matrix produced by the EM algorithm
                  beta,          # initial estimate of the coefficient
                  sigma,         # initial estimate of the variance
                  CandidateModel.list
){
  
  # Record computing time
  start_time_EMAEW = Sys.time()
  
  iter_EMA = 0
  model_index = 31
  number.candidatemodels = 31
  
  data_miss = data[-which(complete.cases(data) == TRUE),]
  data_obs = data[which(complete.cases(data) == TRUE),]
  
  while (iter_EMA < 50){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EMA = list()
    data_gen_EMA_matrix = matrix(0,nrow = dim(data_miss)[1]*50,ncol = 7)
    for (z in 1:dim(data_miss)[1]){
      data_gen_EMA[[z]] = MIS_RD(data_miss[z,],mu = mu, Sigma.matrix = Sigma.matrix,
                                 beta = beta,sigma = sigma,model_index = model_index)
      data_gen_EMA_matrix[c((50*z-49):(50*z)),] = matrix(unlist(data_gen_EMA[[z]]),nrow = 50,ncol = 7,byrow = FALSE)
    }
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_c{Q(M)|y_0} and estimate paramters beta and sigma of different candidate models using the obtained full sample
    parameter_est_EMA = list()
    parameter_est_EMA_matrix = matrix(0,nrow = 6,ncol = 31)
    for (k in 1:number.candidatemodels){
      ## 2.1 Find the index of each candidate model in CandidateModel.list 
      if (k <= 5){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 6)&&(k <= 15)){
        num_covariate = 2
        model_cur_loc = k - 5
      }else if ((k >= 16)&&(k <= 25)){
        num_covariate = 3
        model_cur_loc = k - 15
      }else if ((k >= 26)&&(k <= 30)){
        num_covariate = 4
        model_cur_loc = k - 25
      }else{
        num_covariate = 5
        model_cur_loc = 1
      }
      
      ## 2.2 Compute (1) S_0 = \sum_{i=1}^n E_c(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_c(X_{i,c}Y_i|y_0,x_0)
      ##             (3) S_2 = \sum_{i=1}^n E_c(X_{i,c}X_{i,c}^\prime|y_0,x_0)
      S_0 = 0
      S_1 = matrix(0,nrow = num_covariate + 1,ncol = 1)
      S_2 = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
      for (i in 1:length(data_gen_EMA)){
        data_gen_tem = matrix(unlist(data_gen_EMA[[i]]),nrow = 50, ncol = 7, byrow = FALSE)
        S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
        S_1_tem = matrix(0,nrow = num_covariate + 1,ncol = 1)
        S_2_tem = matrix(0,nrow = num_covariate + 1,ncol = num_covariate + 1)
        for (j in 1:50){
          S_1_tem = S_1_tem + data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_gen_tem[j,1]
          S_2_tem = S_2_tem + data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
            t(data_gen_tem[j,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
        }
        S_1 = S_1 + S_1_tem/50
        S_2 = S_2 + S_2_tem/50
      }
      for (i in 1:(dim(data)[1] - length(data_gen_EMA))){
        S_0 = S_0 + data_obs[i,1]^2
        S_1 = S_1 + data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]*data_obs[i,1]
        S_2 = S_2 + data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)]%*%
          t(data_obs[i,c(2,c(CandidateModel.list[[num_covariate]][,model_cur_loc]) + 1)])
      }
      
      ## 2.3 Estimation of beta and sigma^2
      beta_est_tem = solve(S_2)%*%S_1
      beta_est_full = rep(0,6)
      beta_est_full[c(1,CandidateModel.list[[num_covariate]][,model_cur_loc])] = beta_est_tem
      sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/dim(data)[1]
      parameter_est_EMA_matrix[,k] = beta_est_full
      parameter_est_EMA[[k]] = list(beta_est_full,sigma_square_est)
    }
    
    # 3. Obtain the model averaging estimator with equal weights
    ## Calculate the EW estimator for beta
    beta_MA = parameter_est_EMA_matrix%*%rep(1/number.candidatemodels,number.candidatemodels)
    
    # 4. Replace beta and sigma
    beta = beta_MA
    sigma = unlist(parameter_est_EMA[[31]][2])
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.001)&&(norm(beta.previous - beta,"2")<0.001)){
      break
    }else{}
    
    iter_EMA = iter_EMA + 1
  }
  
  # Record computing time
  end_time_EMAEW = Sys.time()
  running_time_EMAEW = end_time_EMAEW - start_time_EMAEW
  
  return(list(beta,sigma,running_time_EMAEW))
}


# ----------------------------- EM algorithm using the full model -----------------------------
EM_F = function(data,          # data set, responses of all subjects
                mu,            # estimate of mu produced by the EM algorithm
                Sigma.matrix,  # estimate of covariance matrix produced by the EM algorithm
                beta,          # initial estimate of the coefficient
                sigma          # initial estimate of the variance
){
  
  # Record computing time
  start_time_EMF = Sys.time()
  
  iter_EM_F = 0
  model_index = 31
  
  data_miss = data[-which(complete.cases(data) == TRUE),]
  data_obs = data[which(complete.cases(data) == TRUE),]
  
  while (iter_EM_F < 50){
    
    # 1. Resample using MIS.R to obtain the resampling result under the current combined model 
    data_gen_EM = apply(data_miss, MARGIN = 1, MIS_RD, mu = mu,
                        Sigma.matrix = Sigma.matrix,beta = beta,sigma = sigma, model_index = model_index)
    
    beta.previous = beta
    sigma.previous = sigma
    
    # 2. Calculate E_f{Q(M)|y_0} and estimate paramters beta and sigma under the full model
    ## 2.1 Compute (1) S_0 = \sum_{i=1}^n E_f(Y_i^2|y_0,x_0); (2) S_1 = \sum_{i=1}^n E_f(X_{i,c}Y_i|y_0,x_0)
    ##             (3) S_2 = \sum_{i=1}^n E_f(X_{i,c}X_{i,c}^\prime|y_0,x_0)
    S_0 = 0
    num_covariate = 6
    S_1 = matrix(0,nrow = num_covariate,ncol = 1)
    S_2 = matrix(0,nrow = num_covariate,ncol = num_covariate)
    for (i in 1:length(data_gen_EM)){
      data_gen_tem = matrix(unlist(data_gen_EM[[i]]),nrow = 50, ncol = 7, byrow = FALSE)
      S_0 = S_0 + sum(data_gen_tem[,1]^2)/50
      S_1_tem = matrix(0,nrow = num_covariate,ncol = 1)
      S_2_tem = matrix(0,nrow = num_covariate,ncol = num_covariate)
      for (j in 1:50){
        S_1_tem = S_1_tem + data_gen_tem[j,c(2:7)]*data_gen_tem[j,1]
        S_2_tem = S_2_tem + data_gen_tem[j,c(2:7)]%*%t(data_gen_tem[j,c(2:7)])
      }
      S_1 = S_1 + S_1_tem/50
      S_2 = S_2 + S_2_tem/50
    }
    for (i in 1:(dim(data)[1] - length(data_gen_EM))){
      S_0 = S_0 + data_obs[i,1]^2
      S_1 = S_1 + data_obs[i,c(2:7)]*data_obs[i,1]
      S_2 = S_2 + data_obs[i,c(2:7)]%*%t(data_obs[i,c(2:7)])
    }
    
    ## 2.2 Estimation of beta and sigma^2
    beta_est_full = solve(S_2)%*%S_1
    sigma_square_est = (S_0 - t(S_1)%*%solve(S_2)%*%S_1)/dim(data)[1]
    
    # 3. Replace beta and sigma
    beta = beta_est_full
    sigma = sigma_square_est
    
    # Stop criterion
    if ((abs(sigma - sigma.previous)<0.001)&&(norm(beta.previous - beta,"2")<0.001)){
      break
    }else{}
    
    iter_EM_F = iter_EM_F + 1
  }
  
  
  # Record computing time
  end_time_EMF = Sys.time()
  running_time_EMF = end_time_EMF - start_time_EMF
  
  return(list(beta,sigma,running_time_EMF))
}


# ----------------------------- BIC method with observed complete data -----------------------------
BIC_OC = function(missing.data,          # Missing data
                  CandidateModel.list    # Candidate model set
){
  # Initial parameter setting
  # missing.data = train_data
  # CandidateModel.list = CandidateModel.index.list
  
  data_obs = missing.data[which(complete.cases(missing.data)),]
  number.candidatemodels = 31
  
  beta.list = list()
  sigma.summary = matrix(0,nrow = number.candidatemodels,ncol = 1)
  loss.summary = matrix(0,nrow = number.candidatemodels,ncol = 1)
  
  for (k in 1:number.candidatemodels){
    if (k <= 5){
      num_covariate = 1
      model_cur_loc = k
    }else if ((k >= 6)&&(k <= 15)){
      num_covariate = 2
      model_cur_loc = k - 5
    }else if ((k >= 16)&&(k <= 25)){
      num_covariate = 3
      model_cur_loc = k - 15
    }else if ((k >= 26)&&(k <= 30)){
      num_covariate = 4
      model_cur_loc = k - 25
    }else{
      num_covariate = 5
      model_cur_loc = 1
    }
    
    X_tem = data_obs[,c(2,CandidateModel.list[[num_covariate]][,model_cur_loc] + 1)]
    Y = data_obs[,1]
    model = glm(Y~-1 + X_tem,family = gaussian,maxit = 200)
    loss.summary[k,1] = BIC(model)
    sigma.summary[k,1] = norm(residuals(model),"2")/(dim(data_obs)[1] - num_covariate)
    beta_tem = rep(0,6)
    beta_tem[c(1,CandidateModel.list[[num_covariate]][,model_cur_loc])] = coef(model)
    beta.list[[k]] = beta_tem
  }
  
  # Selection with BIC
  model_optimal = which.min(loss.summary)
  beta_optimal = unlist(beta.list[[model_optimal]])
  sigma_optimal = sigma.summary[model_optimal,1]
  
  return(list(beta_optimal,model_optimal,sigma_optimal))
}


# ------------------------ Model selection and averaging based on imputed data ------------------------
MS_Im = function(missing.data,          # Missing data
                 CandidateModel.list    # Candidate model set
){
  
  data_obs = missing.data[which(complete.cases(missing.data)),]
  number.candidatemodels = 31
  
  # Implement imputation to missing data
  initial = mice(missing.data, maxit = 0)
  predM = initial$predictorMatrix
  imputed_data = mice(missing.data, predictorMatrix = predM)
  
  beta.selection = matrix(0, nrow = 6, ncol = 5)
  beta.averaging = matrix(0, nrow = 6, ncol = 5)
  
  for (i in 1:5){
    
    imputed.data.tem = as.matrix(complete(imputed_data, action = i))
    
    selection.cri = matrix(0,nrow = 1, ncol = number.candidatemodels)
    averaging.cri = matrix(0,nrow = 1, ncol = number.candidatemodels)
    beta.selection.tem = matrix(0,nrow = 6, ncol = number.candidatemodels)
    
    for (k in 1:number.candidatemodels){
      if (k <= 5){
        num_covariate = 1
        model_cur_loc = k
      }else if ((k >= 6)&&(k <= 15)){
        num_covariate = 2
        model_cur_loc = k - 5
      }else if ((k >= 16)&&(k <= 25)){
        num_covariate = 3
        model_cur_loc = k - 15
      }else if ((k >= 26)&&(k <= 30)){
        num_covariate = 4
        model_cur_loc = k - 25
      }else{
        num_covariate = 5
        model_cur_loc = 1
      }
      
      X_tem = imputed.data.tem[,c(2,CandidateModel.list[[num_covariate]][,model_cur_loc] + 1)]
      Y = imputed.data.tem[,1]
      model = glm(Y~-1 + X_tem,family = gaussian, maxit = 200)
      selection.cri[1,k] = BIC(model)
      beta_tem = rep(0,6)
      beta_tem[c(1,CandidateModel.list[[num_covariate]][,model_cur_loc])] = coef(model)
      beta.selection.tem[,k] = beta_tem
    }
    
    # Averaging with SBIC
    SBIC = t(matrix(rep(-selection.cri,number.candidatemodels), ncol = number.candidatemodels, byrow = TRUE)
             + c(selection.cri))
    w.SBIC = 1/colSums(apply(SBIC/2, 2,exp))
    # w.SBIC = exp(-selection.cri/2)/sum(exp(-selection.cri/2))
    
    # Selection with BIC
    model_optimal = which.min(selection.cri)
    beta.selection[,i] = beta.selection.tem[,model_optimal]
    beta.averaging[,i] = beta.selection.tem%*%(w.SBIC)
  }
  
  beta_optimal_selection = rowSums(beta.selection)/5
  beta_optimal_averaging = rowSums(beta.averaging)/5
  
  return(list(beta_optimal_selection, beta_optimal_averaging))
}


# -------- Data input and processing --------
data_full = read.csv("AirQualityUCI.csv",header = TRUE)
data = data_full[which(data_full$Time == "08:00:00"),]
missing_count = apply(data, 1, function(row) sum(is.na(row)))
data_to_delete = which(missing_count > 8)
data1 = data[-data_to_delete,c(8,10,3,5,6,13,15)]
data_scale = scale(data1)

# -------- Data splitting --------
complete_data = data_scale[which(complete.cases(data_scale) == TRUE),]
missing_data = data_scale[-which(complete.cases(data_scale) == TRUE),]

# -------- Randomly divide the sample into a training set and a test set --------
Solve_fun = function(n,rep_time,r){
  
  set.seed(1 + r)
  bootstrap_data = complete_data[c(1:10),]
  train_indices_1 = sample(dim(complete_data)[1], 0.8 * dim(complete_data)[1])
  train_indices_2 = sample(dim(missing_data)[1], 0.5 * dim(missing_data)[1])
  train_data = as.matrix(rbind(complete_data[train_indices_1,],missing_data[train_indices_2,]), ncol = 7, byrow = TRUE)
  test_data = as.matrix(complete_data[-train_indices_1,], ncol = 7, byrow = TRUE)
  response_train = train_data[,1]
  covariate_train = train_data[,-1]
  response_train_com = response_train[c(1:length(train_indices_1))]
  covariate_train_com = covariate_train[c(1:length(train_indices_1)),]
  response_test = test_data[,1]
  covariate_test = test_data[,-1]
  num.complete = length(train_indices_1)
  
  # -------- Initial setting for the forthcoming algorithms --------
  # Initial estimates for parameters \beta, sigma_square, \mu, \Pi and Sigma matrix
  beta.ini = solve(t(covariate_train_com)%*%covariate_train_com)%*%t(covariate_train_com)%*%response_train_com
  sigma.square.ini = norm(response_train_com - covariate_train_com%*%beta.ini,type = "F")/(num.complete - dim(covariate_train_com)[2])
  mu.ini = colSums(covariate_train_com)/num.complete
  Sigma.matrix.ini = matrix(rep(0,36), nrow = 6, ncol = 6)
  for (i in 1:num.complete){
    Sigma.matrix.ini = Sigma.matrix.ini + (covariate_train_com[i,] - mu.ini)%*%t(covariate_train_com[i,] - mu.ini)
  }
  Sigma.matrix.ini = Sigma.matrix.ini/num.complete
  
  
  EM.result = EM(train_data, mu = mu.ini, Sigma.matrix = Sigma.matrix.ini)
  
  
  # Build non-nested candidate models
  CandidateModel.index.list = list()
  number.candidatemodels = 0
  for(i in 1:5){
    CandidateModel.index.list[[i]] = combn(2:6,i)
    number.candidatemodels = number.candidatemodels + ncol(CandidateModel.index.list[[i]])
  }
  
  
  
  EMS.result = EMS(data = train_data, mu = EM.result[[1]],Sigma.matrix = EM.result[[2]],
                   beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)
  
  
  EMA_A.result = EMA_A(data = train_data, mu = EM.result[[1]],Sigma.matrix = EM.result[[2]],
                       beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)
  
  
  EMA_B.result = EMA_B(data = train_data, mu = EM.result[[1]],Sigma.matrix = EM.result[[2]],
                       beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)
  
  
  EMA_EW.result = EMA_EW(data = train_data, mu = EM.result[[1]],Sigma.matrix = EM.result[[2]],
                         beta = beta.ini,sigma = sigma.square.ini,CandidateModel.list = CandidateModel.index.list)
  
  
  EM_F.result = EM_F(data = train_data, mu = EM.result[[1]],Sigma.matrix = EM.result[[2]],
                     beta = beta.ini,sigma = sigma.square.ini)
  
  
  BIC.OC.result = BIC_OC(missing.data = train_data, CandidateModel.list = CandidateModel.index.list)
  
  
  MS.Im.result = MS_Im(missing.data = train_data, CandidateModel.list = CandidateModel.index.list)
  
  
  address2 = paste0(address,"results/","tempdata/")
  if(!dir.exists(address2)){
    dir.create(paste0(address,"results/"))
    dir.create(address2)
  }
  
  n_test = length(response_test)
  loss.mu.EMS = norm(response_test - covariate_test%*%c(unlist(EMS.result[[1]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  loss.mu.EMA_A = norm(response_test - covariate_test%*%c(unlist(EMA_A.result[[1]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  loss.mu.EMA_B = norm(response_test - covariate_test%*%c(unlist(EMA_B.result[[1]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  loss.mu.EMA.EW = norm(response_test - covariate_test%*%c(unlist(EMA_EW.result[[1]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  loss.mu.EM.F = norm(response_test - covariate_test%*%c(unlist(EM_F.result[[1]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  loss.mu.BIC.OC = norm(response_test - covariate_test%*%c(unlist(BIC.OC.result[[1]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  loss.mu.MS.Im.SE = norm(response_test - covariate_test%*%c(unlist(MS.Im.result[[1]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  loss.mu.MS.Im.MA = norm(response_test - covariate_test%*%c(unlist(MS.Im.result[[2]])),"2")/n_test - unlist(BIC.OC.result[[3]])
  resultframe = data.frame(t(c(loss.mu.EMS,loss.mu.EMA_A,loss.mu.EMA_B,loss.mu.EMA.EW,loss.mu.EM.F,loss.mu.BIC.OC,loss.mu.MS.Im.SE,loss.mu.MS.Im.MA)))
  
  write.csv(resultframe, paste0(address2,"n=",n,"_rept=",r,".csv"),row.names = FALSE)
  return(resultframe)
}

