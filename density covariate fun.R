# Generate density function for the distribution of covariates under the true model f(x) 
# or under an assumed model g(x).

density_covariate_fun = function(Output,
                                 Pi,
                                 mu,
                                 Sigma.matrix,
                                 beta,
                                 sigma,
                                 miss.index.individual
){
  if ((c(2)%in% miss.index.individual)|(c(3)%in% miss.index.individual)){
    if ((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual)){ 
      # 1. Continuous variables and discrete ones
      if ((c(2)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))&&(!(c(5)%in% miss.index.individual))){ 
        
        # 1.1 The first continuous and discrete variables are missing
        if (Output[5] == 1){ # The value of the observed discrete variable is 1
          if (Output[4] == 1){ # Discrete variables c(1,1)
            result_dis = Pi[4,1]/(Pi[2,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,4])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,1)
            result_dis = Pi[2,1]/(Pi[2,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,2])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }
        }else{ # The value of the observed discrete variable is 0
          if (Output[4] == 1){ # Discrete variables c(1,0)
            result_dis = Pi[3,1]/(Pi[1,1]+Pi[3,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,3])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,0)
            result_dis = Pi[1,1]/(Pi[1,1]+Pi[3,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,1])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))&&(!(c(4)%in% miss.index.individual))){ 
        
        # 1.2 The first continuous variable and the second discrete variable are missing
        if (Output[4] == 1){ # The value of the observed discrete variable is 1
          if (Output[5] == 1){ # Discrete variables c(1,1)
            result_dis = Pi[4,1]/(Pi[3,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,4])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(1,0)
            result_dis = Pi[3,1]/(Pi[3,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,3])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }
        }else{ # The value of the observed discrete variable is 0
          if (Output[5] == 1){ # Discrete variables c(0,1)
            result_dis = Pi[2,1]/(Pi[1,1]+Pi[2,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,2])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,0)
            result_dis = Pi[1,1]/(Pi[1,1]+Pi[2,1])
            result_con = exp(-0.5*(Output[2]-(mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,1])))^2
                             *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
            result = result_dis*result_con
          }
        }
      }else if ((c(3)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(!(c(2)%in% miss.index.individual))&&(!(c(5)%in% miss.index.individual))){ 
        
        # 1.3 The first discrete and the second continuous variables are missing
        if (Output[5] == 1){ # The value of the observed discrete variable (5) is 1
          if (Output[4] == 1){ # Discrete variables c(1,1)
            result_dis = Pi[4,1]/(Pi[2,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,4])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,1)
            result_dis = Pi[2,1]/(Pi[2,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,2])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }
        }else{ # The value of the observed discrete variable is 0
          if (Output[4] == 1){ # Discrete variables c(1,0)
            result_dis = Pi[3,1]/(Pi[1,1]+Pi[3,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,3])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,0)
            result_dis = Pi[1,1]/(Pi[1,1]+Pi[3,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,1])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }
        }
      }else if ((c(3)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(2)%in% miss.index.individual))&&(!(c(4)%in% miss.index.individual))){ 
        
        # 1.4 The second continuous and discrete variables are missing
        if (Output[4] == 1){ # The value of the observed discrete variable (4) is 1
          if (Output[5] == 1){ # Discrete variables c(1,1)
            result_dis = Pi[4,1]/(Pi[3,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,4])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(1,0)
            result_dis = Pi[3,1]/(Pi[3,1]+Pi[4,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,3])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }
        }else{ # The value of the observed discrete variable (4) is 0
          if (Output[5] == 1){ # Discrete variables c(0,1)
            result_dis = Pi[2,1]/(Pi[1,1]+Pi[2,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,2])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,0)
            result_dis = Pi[1,1]/(Pi[1,1]+Pi[2,1])
            result_con = exp(-0.5*(Output[3]-(mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,1])))^2
                             *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
            result = result_dis*result_con
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(!(c(5)%in% miss.index.individual))){ 
        
        # 1.5 The second discrete variable is observed (Other covariates are missing)
        if (Output[5] == 1){
          if (Output[4] == 1){ # Discrete variables c(1,1)
            result_dis = Pi[4,1]/(Pi[2,1]+Pi[4,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,4])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,4]))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,1)
            result_dis = Pi[2,1]/(Pi[2,1]+Pi[4,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,2])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,2]))
            result = result_dis*result_con
          }
        }else{
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),Pi[3,1]/(Pi[1,1]+Pi[3,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,0)
            result_dis = Pi[3,1]/(Pi[1,1]+Pi[3,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,3])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,3]))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,0)
            result_dis = Pi[1,1]/(Pi[1,1]+Pi[3,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,1])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,1]))
            result = result_dis*result_con
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(4)%in% miss.index.individual))){ 
        
        # 1.6 The first discrete variable is observed
        if (Output[4] == 1){
          if (Output[5] == 1){ # Discrete variables c(1,1)
            result_dis = Pi[4,1]/(Pi[3,1]+Pi[4,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,4])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,4]))
            result = result_dis*result_con
          }else{ # Discrete variables c(1,0)
            result_dis = Pi[3,1]/(Pi[3,1]+Pi[4,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,3])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,3]))
            result = result_dis*result_con
          }
        }else{
          if (Output[5] == 1){ # Discrete variables c(0,1)
            result_dis = Pi[2,1]/(Pi[1,1]+Pi[2,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,2])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,2]))
            result = result_dis*result_con
          }else{ # Discrete variables c(0,0)
            result_dis = Pi[1,1]/(Pi[1,1]+Pi[2,1])
            result_con = exp(-0.5*t(Output[c(2,3)] - mu[,1])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,1]))
            result = result_dis*result_con
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))){ 
        
        # 1.7 The second continuous variable is observed
        if ((Output[4] == c(0))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,1])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(0))&&(Output[5] == c(1))){
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,2])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(1))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,3])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }else {
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,4])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }
      }else{ # 1.8 The first continuous variable is observed
        if ((Output[4] == c(0))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,1])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(0))&&(Output[5] == c(1))){
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,2])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(1))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,3])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }else {
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,4])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }
      }
    }else{ # 2. Only continuous variables are missing
      if ((c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)){ # 2.1 Both continuous variables are missing ¡Ì
        if ((Output[4] == c(0))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*t(Output[c(2,3)] - mu[,1])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,1]))
          result = result_dis*result_con
        }else if ((Output[4] == c(0))&&(Output[5] == c(1))){
          result_dis = 1
          result_con = exp(-0.5*t(Output[c(2,3)] - mu[,2])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,2]))
          result = result_dis*result_con
        }else if ((Output[4] == c(1))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*t(Output[c(2,3)] - mu[,3])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,3]))
          result = result_dis*result_con
        }else{
          result_dis = 1
          result_con = exp(-0.5*t(Output[c(2,3)] - mu[,4])%*%solve(Sigma.matrix)%*%(Output[c(2,3)] - mu[,4]))
          result = result_dis*result_con
        }
      }else if((c(2)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))){ 
        # 2.2 The first continuous variable is missing ¡Ì
        if ((Output[4] == c(0))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,1])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(0))&&(Output[5] == c(1))){
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,2])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(1))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,3])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }else{
          result_dis = 1
          result_con = exp(-0.5*(Output[2] - (mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(Output[3]-mu[2,4])))^2
                           *(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])^(-1))
          result = result_dis*result_con
        }
      }else{ # 2.3 The second continuous variable is missing ¡Ì
        if ((Output[4] == c(0))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,1])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(0))&&(Output[5] == c(1))){
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,2])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }else if ((Output[4] == c(1))&&(Output[5] == c(0))){
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,3])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }else{
          result_dis = 1
          result_con = exp(-0.5*(Output[3] - (mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(Output[2]-mu[1,4])))^2
                           *(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])^(-1))
          result = result_dis*result_con
        }
      }
    }
  }else{ # 3. Only discrete variables
    if ((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual)){
      if ((c(4)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)){ 
        # 3.1 Both discrete variables are missing 
        if ((Output[4] == c(0))&&(Output[5] == c(0))){
          result = Pi[1,1]
        }else if ((Output[4] == c(0))&&(Output[5] == c(1))){
          result = Pi[2,1]
        }else if ((Output[4] == c(1))&&(Output[5] == c(0))){
          result = Pi[3,1]
        }else{
          result = Pi[4,1]
        }
      }else if ((c(4)%in% miss.index.individual)&&(!(c(5)%in% miss.index.individual))){ 
        # 3.2 The first discrete variable is missing 
        if (Output[5] == 1){
          if (Output[4] == 1){
            result = Pi[4,1]/(Pi[2,1]+Pi[4,1])
          }else{
            result = Pi[2,1]/(Pi[2,1]+Pi[4,1])
          }
        }else{
          if (Output[4] == 1){
            result = Pi[3,1]/(Pi[1,1]+Pi[3,1])
          }else{
            result = Pi[1,1]/(Pi[1,1]+Pi[3,1])
          }
        }
      }else if((c(5)%in% miss.index.individual)&&(!(c(4)%in% miss.index.individual))){ 
        # 3.3 The second discrete variable is missing ¡Ì
        if (Output[4] == 1){
          if (Output[5] == 1){
            result = Pi[4,1]/(Pi[3,1]+Pi[4,1])
          }else{
            result = Pi[3,1]/(Pi[3,1]+Pi[4,1])
          }
        }else{
          if (Output[5] == 1){
            result = Pi[2,1]/(Pi[1,1]+Pi[2,1])
          }else{
            result = Pi[1,1]/(Pi[1,1]+Pi[2,1])
          }
        }
      }
    }
  }
  return(result)
}