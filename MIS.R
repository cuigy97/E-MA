# Output a set of fully-observed samples that can be directly used in approximating the conditional expectation 
# Two different cases: 1. The resampling only involves covariates;
#                      2. The resampling involves both covariates and responses.

MIS = function(data.full,          # Observations of each individual
               mu,                 # Initial value of mu in the current iteration
               Pi,                 # Initial value of Pi in the current iteration
               Sigma.matrix,       # Initial value of covariance in the current iteration
               beta = c(1000,1000,1000,1000),    # The default value indicates that the MIS only targets the covariates
               sigma = 1000,                     # Variance
               model_index         # Indicate on which candidate model the resampling is based
){
  
  data.individual = as.matrix(data.full)
  nu = matrix(c(0,0,1,1,0,1,0,1),nrow = 2, ncol = 4,byrow = TRUE)
  r.nu = c(1,2,3,4)
  miss.index.individual = c(which(is.na(data.individual)))
  num_sample = 300
  samples = matrix(0,nrow = num_sample, ncol = 5)
  
  
  # Find the location of the current candidate model in CandidateModel.index.list
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
  CandidateModel.index.list = list()
  for(i in 1:4){
    CandidateModel.index.list[[i]] = combn(1:4, i)
  }
  index = CandidateModel.index.list[[num_covariate]][,model_cur_loc]
  index_com = setdiff(c(1:5),c(1,c(index)+1))
  
  
  if (sigma == 1000){
    Output_MIS = c(-20,0,0,1,1)
  }else{
    Output_MIS = c(1,0,0,1,1)
  }
  
  for (i in 1:num_sample){
    if ((beta[1] == 1000)&&(sigma == 1000)){ # Resampling only covariates
      data.tem = MIS_covariate(data.covariate = data.individual,
                             mu = mu,
                             Pi = Pi,
                             Sigma.matrix = Sigma.matrix)
    }else{ # Resampling all variables
      if (1 %in% miss.index.individual){
        if (((c(2)%in% miss.index.individual)|(c(3)%in% miss.index.individual))&&(!(c(4)%in% miss.index.individual))&&(!(c(5)%in% miss.index.individual))){
          
          # 1. Response and continuous variables
          if (c(2)%in% miss.index.individual){
            if (c(3)%in% miss.index.individual){ # Both continuous variables are missing
              if ((data.individual[4,1] == 0)&&(data.individual[5,1] == 0)){
                data.con.tem = mvrnorm(n = 1,mu = mu[,1],Sigma = Sigma.matrix)
                data.tem = c(data.con.tem,c(0,0))
              }else if ((data.individual[4,1] == 0)&&(data.individual[5,1] == 1)){
                data.con.tem = mvrnorm(n = 1,mu = mu[,2],Sigma = Sigma.matrix)
                data.tem = c(data.con.tem,c(0,0))
              }else if ((data.individual[4,1] == 1)&&(data.individual[5,1] == 0)){
                data.con.tem = mvrnorm(n = 1,mu = mu[,3],Sigma = Sigma.matrix)
                data.tem = c(data.con.tem,c(0,0))
              }else{
                data.con.tem = mvrnorm(n = 1,mu = mu[,4],Sigma = Sigma.matrix)
                data.tem = c(data.con.tem,c(0,0))
              }
            }else{ # Only the first continuous variable is missing
              if ((data.individual[4,1] == 0)&&(data.individual[5,1] == 0)){
                data.con.tem = rnorm(n = 1,mean = mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,1]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2]))
                data.tem = c(data.con.tem,data.individual[3,1],c(0,0))
              }else if ((data.individual[4,1] == 0)&&(data.individual[5,1] == 1)){
                data.con.tem = rnorm(n = 1,mean = mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,2]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2]))
                data.tem = c(data.con.tem,data.individual[3,1],c(0,1))
              }else if ((data.individual[4,1] == 1)&&(data.individual[5,1] == 0)){
                data.con.tem = rnorm(n = 1,mean = mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,3]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2]))
                data.tem = c(data.con.tem,data.individual[3,1],c(1,0))
              }else{
                data.con.tem = rnorm(n = 1,mean = mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,4]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2]))
                data.tem = c(data.con.tem,data.individual[3,1],c(1,1))
              }
            }
          }else{ # Only the second continuous variable is missing
            if ((data.individual[4,1] == 0)&&(data.individual[5,1] == 0)){
              data.con.tem = rnorm(n = 1,mean = mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,1]),
                                   sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1]))
              data.tem = c(data.individual[2,1],data.con.tem,c(0,0))
            }else if ((data.individual[4,1] == 0)&&(data.individual[5,1] == 1)){
              data.con.tem = rnorm(n = 1,mean = mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,2]),
                                   sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1]))
              data.tem = c(data.individual[2,1],data.con.tem,c(0,1))
            }else if ((data.individual[4,1] == 1)&&(data.individual[5,1] == 0)){
              data.con.tem = rnorm(n = 1,mean = mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,3]),
                                   sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1]))
              data.tem = c(data.individual[2,1],data.con.tem,c(1,0))
            }else{
              data.con.tem = rnorm(n = 1,mean = mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,4]),
                                   sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1]))
              data.tem = c(data.individual[2,1],data.con.tem,c(1,1))
            }
          }
        }else if (((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual))&&(!(c(2)%in% miss.index.individual))&&(!(c(3)%in% miss.index.individual))){ 
          
          # 2. Response and discrete variables
          if (c(4)%in% miss.index.individual){
            if (c(5)%in% miss.index.individual){ # Both discrete variables are missing
              data.discrete.tem = nu[,sample(r.nu,size = 1,replace = TRUE,prob = Pi)]
              data.tem = c(data.individual[c(2,3),1],data.discrete.tem)
            }else{ # The first discrete variable is missing
              if (data.individual[5,] == 1){
                data.discrete.tem = sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[2,1]/(Pi[2,1]+Pi[4,1]),Pi[4,1]/(Pi[2,1]+Pi[4,1])))
                data.tem = c(data.individual[c(2,3),1],data.discrete.tem,data.individual[5,1])
              }else{
                data.discrete.tem = sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),Pi[3,1]/(Pi[1,1]+Pi[3,1])))
                data.tem = c(data.individual[c(2,3),1],data.discrete.tem,data.individual[5,1])
              }
            }
          }else{ # The second discrete variable is missing
            if (data.individual[4,] == 1){
              data.discrete.tem = sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[3,1]/(Pi[3,1]+Pi[4,1]),Pi[4,1]/(Pi[3,1]+Pi[4,1])))
              data.tem = c(data.individual[c(2:4),1],data.discrete.tem)
            }else{
              data.discrete.tem = sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[2,1]),Pi[2,1]/(Pi[1,1]+Pi[2,1])))
              data.tem = c(data.individual[c(2:4),1],data.discrete.tem)
            }
          }
        }else if (((c(2)%in% miss.index.individual)|(c(3)%in% miss.index.individual))&&((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual))){
          
          # 3. Response, continuous and discrete variables
          if (c(4)%in% miss.index.individual){
            if (c(2)%in% miss.index.individual){ # 2,4
              if (data.individual[5,1] == 1){ # The value of the observed discrete variable is 1
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[2,1]/(Pi[2,1]+Pi[4,1]),Pi[4,1]/(Pi[2,1]+Pi[4,1])))
                if (data.discrete.tem == 1){ # Discrete variables c(1,1)
                  data.tem = c(rnorm(n = 1,mean = mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,4]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,1))
                }else{ # Discrete variables c(0,1)
                  data.tem = c(rnorm(n = 1,mean = mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,2]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,1))
                }
              }else{ # The value of the observed discrete variable is 0
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),(Pi[3,1]/(Pi[1,1]+Pi[3,1]))))
                if (data.discrete.tem == 1){ # Discrete variables c(1,0)
                  data.tem = c(rnorm(n = 1,mean = mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,3]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,0))
                }else{ # Discrete variables c(0,0)
                  data.tem = c(rnorm(n = 1,mean = mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,1]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,0))
                }
              }
            }else{ # 3,4
              if (data.individual[5,1] == 1){ # The value of the observed discrete variable (5) is 1
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[2,1]/(Pi[2,1]+Pi[4,1]),Pi[4,1]/(Pi[2,1]+Pi[4,1])))
                if (data.discrete.tem == 1){ # Discrete variables c(1,1)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,4]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,1))
                }else{ # Discrete variables c(0,1)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,2]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,1))
                }
              }else{ # The value of the observed discrete variable is 0
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),(Pi[3,1]/(Pi[1,1]+Pi[3,1]))))
                if (data.discrete.tem == 1){ # Discrete variables c(1,0)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,3]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,0))
                }else{ # Discrete variables c(0,0)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,1]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,0))
                }
              }
            }
          }else{
            if (c(2)%in% miss.index.individual){ # 2,5
              if (data.individual[4,1] == 1){ # The value of the observed discrete variable is 1
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[3,1]/(Pi[3,1]+Pi[4,1]),Pi[4,1]/(Pi[3,1]+Pi[4,1])))
                if (data.discrete.tem == 1){ # Discrete variables c(1,1)
                  data.tem = c(rnorm(n = 1,mean = mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,4]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,1))
                }else{ # Discrete variables c(1,0)
                  data.tem = c(rnorm(n = 1,mean = mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,3]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,0))
                }
              }else{ # The value of the observed discrete variable is 0
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[2,1]),Pi[2,1]/(Pi[1,1]+Pi[2,1])))
                if (data.discrete.tem == 1){ # Discrete variables c(0,1)
                  data.tem = c(rnorm(n = 1,mean = mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,2]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,1))
                }else{ # Discrete variables c(0,0)
                  data.tem = c(rnorm(n = 1,mean = mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,1]),
                                     sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,0))
                }
              }
            }else{ # 3,5
              if (data.individual[4,1] == 1){ # The value of the observed discrete variable (4) is 1
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[3,1]/(Pi[3,1]+Pi[4,1]),Pi[4,1]/(Pi[3,1]+Pi[4,1])))
                if (data.discrete.tem == 1){ # Discrete variables c(1,1)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,4]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,1))
                }else{ # Discrete variables c(1,0)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,3]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,0))
                }
              }else{ # The value of the observed discrete variable (4) is 0
                data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[2,1]),Pi[2,1]/(Pi[1,1]+Pi[2,1])))
                if (data.discrete.tem == 1){ # Discrete variables c(0,1)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,2]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,1))
                }else{ # Discrete variables c(0,0)
                  data.tem = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,1]),
                                                          sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,0))
                }
              }
            }
          }
        }else{ # 4. Response ¡Ì
          data.tem = data.individual[-1,]
        }
      }else{ # Response not included
        data.tem = MIS_covariate(data.covariate = data.individual,
                                 mu = mu,
                                 Pi = Pi,
                                 Sigma.matrix = Sigma.matrix)
      }
    }
    
    Output_tem = c(0,data.tem)
    Output_tem[index_com] = 0
    
    
    if (sigma == 1000){ 
      # Resample covariates only
      Output = data.tem
    }else{ 
      # Resample all variables
      if (c(1)%in% miss.index.individual){ 
        # Response missing
        y.tem = rnorm(n = 1, mean = t(Output_tem[c(2:5)])%*%beta, sd = sqrt(sigma))
        Output = c(y.tem,data.tem)
      }else{ 
        # Response observed
        Output = c(data.individual[1,1],data.tem)
      }
    }
    
    # MIS algorithm
    ## Extend the dimension of Output to 5 in order to implement density function_full.R
    if (length(Output) == 4){
      Output = c(-20,Output)
    }
    ## Implement MIS algorithm
    a1 = density_fun(Output,Pi,mu,Sigma.matrix,beta,sigma,miss.index.individual,model_index,target = 1)
    a2 = density_fun(Output,Pi,mu,Sigma.matrix,beta,sigma,miss.index.individual,model_index,target = 0)
    a3 = density_fun(Output_MIS,Pi,mu,Sigma.matrix,beta,sigma,miss.index.individual,model_index,target = 1)
    a4 = density_fun(Output_MIS,Pi,mu,Sigma.matrix,beta,sigma,miss.index.individual,model_index,target = 0)
    acceptance_prob = min(1,(a1/a2)/(a3/a4))
    if (runif(1) < acceptance_prob) {
      Output_MIS = Output
    }
    samples[i,] = Output_MIS
  }
  # }
  
  # Return Output_MIS
  if (Output_MIS[1] == -20){
    samples = samples[c(251:300),-1]
  }else{
    samples = samples[c(251:300),]
  }
  return(list(samples))
}


