# Generate density function for the joint distribution of the response and covariates under 
# the true model f(x) or under an assumed model g(x).

density_fun = function(Output,
                       Pi,
                       mu,
                       Sigma.matrix,
                       beta,
                       sigma,
                       miss.index.individual,   # Index indicating the missing values of the original individual data
                       model_index,             # Indicate on which model the calculation is based
                       target                   # target = 1: f(x); target = 0: g(x)
){
  
  if (target == 1){ # f(x)
    if (sigma == 1000){ # Resample covariates only
        if ((Output[4] == 0)&&(Output[5] == 0)){
          result_covariate = exp(log(Pi[1])-0.5*(t(Output[c(2,3)] - mu[,1])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,1])))
        }else if ((Output[4] == 0)&&(Output[5] == 1)){
          result_covariate = exp(log(Pi[2])-0.5*(t(Output[c(2,3)] - mu[,2])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,2])))
        }else if ((Output[4] == 1)&&(Output[5] == 0)){
          result_covariate = exp(log(Pi[3])-0.5*(t(Output[c(2,3)] - mu[,3])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,3])))
        }else{
          result_covariate = exp(log(Pi[4])-0.5*(t(Output[c(2,3)] - mu[,4])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,4])))
        }
    }else{ # Resample all variables
        if ((Output[4] == 0)&&(Output[5] == 0)){
          result_covariate = exp(log(Pi[1])-0.5*(t(Output[c(2,3)] - mu[,1])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,1])))
        }else if ((Output[4] == 0)&&(Output[5] == 1)){
          result_covariate = exp(log(Pi[2])-0.5*(t(Output[c(2,3)] - mu[,2])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,2])))
        }else if ((Output[4] == 1)&&(Output[5] == 0)){
          result_covariate = exp(log(Pi[3])-0.5*(t(Output[c(2,3)] - mu[,3])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,3])))
        }else{
          result_covariate = exp(log(Pi[4])-0.5*(t(Output[c(2,3)] - mu[,4])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,4])))
        }
    }
  }else{ # g(x)
    if (sigma == 1000){ # Resample covariates only
      result_covariate = density_covariate_fun(Output,
                                     Pi,
                                     mu,
                                     Sigma.matrix,
                                     beta,
                                     sigma,
                                     miss.index.individual
                                     )
    }else{ # Resample all variables
      
      if (1 %in% miss.index.individual){
        if (((c(2)%in% miss.index.individual)|(c(3)%in% miss.index.individual))&&(!(c(4)%in% miss.index.individual))&&(!(c(5)%in% miss.index.individual))){
          
          # 1. Response and continuous variables
          if (c(2)%in% miss.index.individual){
            if (c(3)%in% miss.index.individual){ # Both continuous variables are missing
              if ((Output[4] == 0)&&(Output[5] == 0)){
                result_dis = 1
                result_con = exp(-0.5*t(Output[c(2,3)] - mu[,1])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,1]))
                result_covariate = result_dis*result_con
              }else if ((Output[4] == 0)&&(Output[5] == 1)){
                result_dis = 1
                result_con = exp(-0.5*t(Output[c(2,3)] - mu[,2])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,2]))
                result_covariate = result_dis*result_con
              }else if ((Output[4] == 1)&&(Output[5] == 0)){
                result_dis = 1
                result_con = exp(-0.5*t(Output[c(2,3)] - mu[,3])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,3]))
                result_covariate = result_dis*result_con
              }else{
                result_dis = 1
                result_con = exp(-0.5*t(Output[c(2,3)] - mu[,4])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,4]))
                result_covariate = result_dis*result_con
              }
            }else{ # Only the first continuous variable is missing
              if ((Output[4] == 0)&&(Output[5] == 0)){
                result_dis = 1
                result_con = exp(-0.5*(Output[2] - (mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,1])))^2
                                 *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                result_covariate = result_dis*result_con
              }else if ((Output[4] == 0)&&(Output[5] == 1)){
                result_dis = 1
                result_con = exp(-0.5*(Output[2] - (mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,2])))^2
                                 *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                result_covariate = result_dis*result_con
              }else if ((Output[4] == 1)&&(Output[5] == 0)){
                result_dis = 1
                result_con = exp(-0.5*(Output[2] - (mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,3])))^2
                                 *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                result_covariate = result_dis*result_con
              }else{
                result_dis = 1
                result_con = exp(-0.5*(Output[2] - (mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,4])))^2
                                 *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                result_covariate = result_dis*result_con
              }
            }
          }else{ # Only the second continuous variable is missing
            if ((Output[4] == 0)&&(Output[5] == 0)){
              result_dis = 1
              result_con = exp(-0.5*(Output[3] - (mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,1])))^2
                               *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
              result_covariate = result_dis*result_con
            }else if ((Output[4] == 0)&&(Output[5] == 1)){
              result_dis = 1
              result_con = exp(-0.5*(Output[3] - (mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,2])))^2
                               *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
              result_covariate = result_dis*result_con
            }else if ((Output[4] == 1)&&(Output[5] == 0)){
              result_dis = 1
              result_con = exp(-0.5*(Output[3] - (mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,3])))^2
                               *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
              result_covariate = result_dis*result_con
            }else{
              result_dis = 1
              result_con = exp(-0.5*(Output[3] - (mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,4])))^2
                               *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
              result_covariate = result_dis*result_con
            }
          }
        }else if (((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual))&&(!(c(2)%in% miss.index.individual))&&(!(c(3)%in% miss.index.individual))){ 
          
          # 2. Response and discrete variables
          if (c(4)%in% miss.index.individual){
            if (c(5)%in% miss.index.individual){ # Both discrete variables are missing
              if ((Output[4] == c(0))&&(Output[5] == c(0))){
                result_dis = Pi[1,1]
                result_con = 1
                result_covariate = result_dis*result_con
              }else if ((Output[4] == c(0))&&(Output[5] == c(1))){
                result_dis = Pi[1,1]
                result_con = 1
                result_covariate = result_dis*result_con
              }else if ((Output[4] == c(1))&&(Output[5] == c(0))){
                result_dis = Pi[1,1]
                result_con = 1
                result_covariate = result_dis*result_con
              }else{
                result_dis = Pi[1,1]
                result_con = 1
                result_covariate = result_dis*result_con
              }
            }else{ # The first discrete variable is missing
              if (Output[5] == 1){
                if (Output[4] == 1){
                  result_dis = Pi[4,1]/(Pi[2,1]+Pi[4,1])
                  result_con = 1
                  result_covariate = result_dis*result_con
                }else{
                  result_dis = Pi[2,1]/(Pi[2,1]+Pi[4,1])
                  result_con = 1
                  result_covariate = result_dis*result_con
                }
              }else{
                if (Output[4] == 1){
                  result_dis = Pi[3,1]/(Pi[1,1]+Pi[3,1])
                  result_con = 1
                  result_covariate = result_dis*result_con
                }else{
                  result_dis = Pi[1,1]/(Pi[1,1]+Pi[3,1])
                  result_con = 1
                  result_covariate = result_dis*result_con
                }
              }
            }
          }else{ # The second discrete variable is missing
            if (Output[4] == 1){
              if (Output[5] == 1){
                result_dis = Pi[4,1]/(Pi[3,1]+Pi[4,1])
                result_con = 1
                result_covariate = result_dis*result_con
              }else{
                result_dis = Pi[3,1]/(Pi[3,1]+Pi[4,1])
                result_con = 1
                result_covariate = result_dis*result_con
              }
            }else{
              if (Output[5] == 1){
                result_dis = Pi[2,1]/(Pi[1,1]+Pi[2,1])
                result_con = 1
                result_covariate = result_dis*result_con
              }else{
                result_dis = Pi[1,1]/(Pi[1,1]+Pi[2,1])
                result_con = 1
                result_covariate = result_dis*result_con
              }
            }
          }
        }else if (((c(2)%in% miss.index.individual)|(c(3)%in% miss.index.individual))&&((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual))){
          
          # 3. Response, continuous and discrete variables
          if (c(4)%in% miss.index.individual){
            if (c(2)%in% miss.index.individual){ # 2,4 missing
              if (Output[5] == 1){ # The value of variable 5 is 1
                if (Output[4] == 1){ # Discrete variables c(1,1)
                  result_dis = Pi[4,1]/(Pi[2,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,4])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(0,1)
                  result_dis = Pi[2,1]/(Pi[2,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,2])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }
              }else{ # The value of the observed discrete variable is 0
                if (Output[4] == 1){ # Discrete variables c(1,0)
                  result_dis = Pi[3,1]/(Pi[1,1]+Pi[3,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,3])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(0,0)
                  result_dis = Pi[1,1]/(Pi[1,1]+Pi[3,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,1])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }
              }
            }else{ # 3,4 missing
              if (Output[5] == 1){ # The value of the observed discrete variable (5) is 1
                if (Output[4] == 1){ # Discrete variables c(1,1)
                  result_dis = Pi[4,1]/(Pi[2,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,4])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(0,1)
                  result_dis = Pi[2,1]/(Pi[2,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,2])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }
              }else{ # The value of the observed discrete variable is 0
                if (Output[4] == 1){ # Discrete variables c(1,0)
                  result_dis = Pi[3,1]/(Pi[1,1]+Pi[3,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,3])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(0,0)
                  result_dis = Pi[1,1]/(Pi[1,1]+Pi[3,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,1])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }
              }
            }
          }else{
            if (c(2)%in% miss.index.individual){ # 2,5
              if (Output[4] == 1){ # The value of the observed discrete variable is 1
                if (Output[5] == 1){ # Discrete variables c(1,1)
                  result_dis = Pi[4,1]/(Pi[3,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,4])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(1,0)
                  result_dis = Pi[3,1]/(Pi[3,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,3])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }
              }else{ # The value of the observed discrete variable is 0
                if (Output[5] == 1){ # Discrete variables c(0,1)
                  result_dis = Pi[2,1]/(Pi[1,1]+Pi[2,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,2])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(0,0)
                  result_dis = Pi[1,1]/(Pi[1,1]+Pi[2,1])
                  result_con = exp(-0.5*(Output[2]-(mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,1])))^2
                                   *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
                  result_covariate = result_dis*result_con
                }
              }
            }else{ # 3,5
              if (Output[4] == 1){ # The value of the observed discrete variable (4) is 1
                if (Output[5] == 1){ # Discrete variables c(1,1)
                  result_dis = Pi[4,1]/(Pi[3,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,4])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(1,0)
                  result_dis = Pi[3,1]/(Pi[3,1]+Pi[4,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,3])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }
              }else{ # The value of the observed discrete variable (4) is 0
                if (Output[5] == 1){ # Discrete variables c(0,1)
                  result_dis = Pi[2,1]/(Pi[1,1]+Pi[2,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,2])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }else{ # Discrete variables c(0,0)
                  result_dis = Pi[1,1]/(Pi[1,1]+Pi[2,1])
                  result_con = exp(-0.5*(Output[3]-(mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,1])))^2
                                   *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
                  result_covariate = result_dis*result_con
                }
              }
            }
          }
        }else{ # 4. Response
          result_dis = 1
          result_con = 1
          result_covariate = result_dis*result_con
        }
      }else{ # Response not included
        result_covariate = density_covariate_fun(Output,
                                       Pi,
                                       mu,
                                       Sigma.matrix,
                                       beta,
                                       sigma,
                                       miss.index.individual
        )
      }
    }
  }
  
  if (model_index <= 4){
    num_covariate = 1
    model_cur_loc = model_index
  }else if ((model_index >= 4)&&(model_index <= 10)){
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
  number.candidatemodels = 0
  for(i in 1:4){
    CandidateModel.index.list[[i]] = combn(1:4, i)
    number.candidatemodels = number.candidatemodels + ncol(CandidateModel.index.list[[i]])
  }
  index = CandidateModel.index.list[[num_covariate]][,model_cur_loc]
  index_com = setdiff(c(1:5),c(1,c(index)+1))
  Output[index_com] = 0
  
  if (sigma == 1000){ # Resample covariates only
    result_res = 1
    result = result_covariate*result_res
  }else{ # Resample all variables
    if (c(1)%in% miss.index.individual){ # Response missing
      result_res = exp(-(Output[1] - t(Output[c(2:5)])%*%beta)^2/(2*sigma))
      result = result_covariate*result_res
    }else{ # Response observed
      result_res = 1
      result = result_covariate*result_res
    }
  }
  return(result)
}
