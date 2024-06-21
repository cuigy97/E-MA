# Output a set of generated covariates that can be directly used to approximate conditional expectations.

MIS_covariate <- function(data.covariate,   # Observations of each individual
                         mu,                 # Initial value of mu in the current iteration
                         Pi,                 # Initial value of Pi in the current iteration
                         Sigma.matrix
){
  
  data.individual = as.matrix(data.covariate)
  nu = matrix(c(0,0,1,1,0,1,0,1),nrow = 2, ncol = 4,byrow = TRUE)
  r.nu = c(1,2,3,4)
  miss.index.individual = c(which(is.na(data.individual)))
  Output = matrix(0,nrow = 1,ncol = 4)
  
  if ((c(2)%in% miss.index.individual)|(c(3)%in% miss.index.individual)){
    if ((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual)){ 
      # 1. Continuous variables and discrete ones
      if ((c(2)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))&&(!(c(5)%in% miss.index.individual))){ 
        
        # 1.1 The first continuous and discrete variables are missing
        if (data.individual[5,1] == 1){ # The value of the observed discrete variable is 1
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[2,1]/(Pi[2,1]+Pi[4,1]),Pi[4,1]/(Pi[2,1]+Pi[4,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,1)
            Output = c(rnorm(n = 1,mean = mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,4]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,1))
          }else{ # Discrete variables c(0,1)
            Output = c(rnorm(n = 1,mean = mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,2]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,1))
          }
        }else{ # The value of the observed discrete variable is 0
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),(Pi[3,1]/(Pi[1,1]+Pi[3,1]))))
          if (data.discrete.tem == 1){ # Discrete variables c(1,0)
            Output = c(rnorm(n = 1,mean = mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,3]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,0))
          }else{ # Discrete variables c(0,0)
            Output = c(rnorm(n = 1,mean = mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,1]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,0))
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))&&(!(c(4)%in% miss.index.individual))){ 
        
        # 1.2 The first continuous variable and the second discrete variable are missing
        if (data.individual[4,1] == 1){ # The value of the observed discrete variable is 1
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[3,1]/(Pi[3,1]+Pi[4,1]),Pi[4,1]/(Pi[3,1]+Pi[4,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,1)
            Output = c(rnorm(n = 1,mean = mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,4]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,1))
          }else{ # Discrete variables c(1,0)
            Output = c(rnorm(n = 1,mean = mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,3]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,0))
          }
        }else{ # The value of the observed discrete variable is 0
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[2,1]),Pi[2,1]/(Pi[1,1]+Pi[2,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(0,1)
            Output = c(rnorm(n = 1,mean = mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,2]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,1))
          }else{ # Discrete variables c(0,0)
            Output = c(rnorm(n = 1,mean = mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,1]),
                             sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,0))
          }
        }
      }else if ((c(3)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(!(c(2)%in% miss.index.individual))&&(!(c(5)%in% miss.index.individual))){ 
        
        # 1.3 The first discrete and the second continuous variables are missing
        if (data.individual[5,1] == 1){ # The value of the observed discrete variable (5) is 1
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[2,1]/(Pi[2,1]+Pi[4,1]),Pi[4,1]/(Pi[2,1]+Pi[4,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,1)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,4]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,1))
          }else{ # Discrete variables c(0,1)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,2]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,1))
          }
        }else{ # The value of the observed discrete variable is 0
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),(Pi[3,1]/(Pi[1,1]+Pi[3,1]))))
          if (data.discrete.tem == 1){ # Discrete variables c(1,0)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,3]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,0))
          }else{ # Discrete variables c(0,0)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,1]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,0))
          }
        }
      }else if ((c(3)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(2)%in% miss.index.individual))&&(!(c(4)%in% miss.index.individual))){ 
        
        # 1.4 The second continuous and discrete variables are missing
        if (data.individual[4,1] == 1){ # The value of the observed discrete variable (4) is 1
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[3,1]/(Pi[3,1]+Pi[4,1]),Pi[4,1]/(Pi[3,1]+Pi[4,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,1)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,4]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,1))
          }else{ # Discrete variables c(1,0)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,3]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,0))
          }
        }else{ # The value of the observed discrete variable (4) is 0
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[2,1]),Pi[2,1]/(Pi[1,1]+Pi[2,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(0,1)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,2]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,1))
          }else{ # Discrete variables c(0,0)
            Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,1]),
                                                  sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,0))
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(!(c(5)%in% miss.index.individual))){ 
        
        # 1.5 The second discrete variable is observed
        if (data.individual[5,1] == 1){
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[2,1]/(Pi[2,1]+Pi[4,1]),Pi[4,1]/(Pi[2,1]+Pi[4,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,1)
            Output = c(mvrnorm(n = 1,mu = mu[,4],Sigma = Sigma.matrix),c(1,1))
          }else{ # Discrete variables c(0,1)
            Output = c(mvrnorm(n = 1,mu = mu[,2],Sigma = Sigma.matrix),c(0,1))
          }
        }else{
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),Pi[3,1]/(Pi[1,1]+Pi[3,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,0)
            Output = c(mvrnorm(n = 1,mu = mu[,3],Sigma = Sigma.matrix),c(1,0))
          }else{ # Discrete variables c(0,0)
            Output = c(mvrnorm(n = 1,mu = mu[,1],Sigma = Sigma.matrix),c(0,0))
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(4)%in% miss.index.individual))){ 
        
        # 1.6 The first discrete variable is observed
        if (data.individual[4,1] == 1){
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[3,1]/(Pi[3,1]+Pi[4,1]),Pi[4,1]/(Pi[3,1]+Pi[4,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(1,1)
            Output = c(mvrnorm(n = 1,mu = mu[,4],Sigma = Sigma.matrix),c(1,1))
          }else{ # Discrete variables c(1,0)
            Output = c(mvrnorm(n = 1,mu = mu[,3],Sigma = Sigma.matrix),c(0,1))
          }
        }else{
          data.discrete.tem = sample(c(0,1),size = 1,replace = FALSE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[2,1]),Pi[2,1]/(Pi[1,1]+Pi[2,1])))
          if (data.discrete.tem == 1){ # Discrete variables c(0,1)
            Output = c(mvrnorm(n = 1,mu = mu[,2],Sigma = Sigma.matrix),c(0,1))
          }else{ # Discrete variables c(0,0)
            Output = c(mvrnorm(n = 1,mu = mu[,1],Sigma = Sigma.matrix),c(0,0))
          }
        }
      }else if ((c(2)%in% miss.index.individual)&&(c(4)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))){ 
        
        # 1.7 The second continuous variable is observed
        data.discrete.tem = nu[,sample(r.nu,size = 1,replace = TRUE,prob = Pi)]
        if ((data.discrete.tem[1] == c(0))&&(data.discrete.tem[2] == c(0))){
          Output = c(rnorm(n = 1,mean = mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,1]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,0))
        }else if ((data.discrete.tem[1] == c(0))&&(data.discrete.tem[2] == c(1))){
          Output = c(rnorm(n = 1,mean = mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,2]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(0,1))
        }else if ((data.discrete.tem[1] == c(1))&&(data.discrete.tem[2] == c(0))){
          Output = c(rnorm(n = 1,mean = mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,3]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,0))
        }else {
          Output = c(rnorm(n = 1,mean = mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,4]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[3,1],c(1,1))
        }
      }else{ # 1.8 The first continuous variable is observed
        data.discrete.tem = nu[,sample(r.nu,size = 1,replace = TRUE,prob = Pi)]
        if ((data.discrete.tem[1] == c(0))&&(data.discrete.tem[2] == c(0))){
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,1]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,0))
        }else if ((data.discrete.tem[1] == c(0))&&(data.discrete.tem[2] == c(1))){
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,2]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(0,1))
        }else if ((data.discrete.tem[1] == c(1))&&(data.discrete.tem[2] == c(0))){
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,3]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,0))
        }else {
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,4]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),c(1,1))
        }
      }
    }else{ # 2. Only continuous variables 
      if ((c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)){ 
        # 2.1 Both continuous variables are missing
        if ((data.individual[4,1] == c(0))&&(data.individual[5,1] == c(0))){
          Output = c(mvrnorm(n = 1,mu = mu[,1],Sigma = Sigma.matrix),data.individual[c(4,5),1])
        }else if ((data.individual[4,1] == c(0))&&(data.individual[5,1] == c(1))){
          Output = c(mvrnorm(n = 1,mu = mu[,2],Sigma = Sigma.matrix),data.individual[c(4,5),1])
        }else if ((data.individual[4,1] == c(1))&&(data.individual[5,1] == c(0))){
          Output = c(mvrnorm(n = 1,mu = mu[,3],Sigma = Sigma.matrix),data.individual[c(4,5),1])
        }else{
          Output = c(mvrnorm(n = 1,mu = mu[,4],Sigma = Sigma.matrix),data.individual[c(4,5),1])
        }
      }else if((c(2)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))){ # 2.2 The first continuous variable is missing ¡Ì
        if ((data.individual[4,1] == c(0))&&(data.individual[5,1] == c(0))){
          Output = c(rnorm(n = 1,mean = mu[1,1]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,1]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[c(3:5),1])
        }else if ((data.individual[4,1] == c(0))&&(data.individual[5,1] == c(1))){
          Output = c(rnorm(n = 1,mean = mu[1,2]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,2]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[c(3:5),1])
        }else if ((data.individual[4,1] == c(1))&&(data.individual[5,1] == c(0))){
          Output = c(rnorm(n = 1,mean = mu[1,3]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,3]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[c(3:5),1])
        }else{
          Output = c(rnorm(n = 1,mean = mu[1,4]+Sigma.matrix[1,2]/Sigma.matrix[2,2]*(data.individual[3,1]-mu[2,4]),
                           sd = sqrt(Sigma.matrix[1,1]-Sigma.matrix[1,2]^2/Sigma.matrix[2,2])),data.individual[c(3:5),1])
        }
      }else{ # 2.3 The second continuous variable is missing
        if ((data.individual[4,1] == c(0))&&(data.individual[5,1] == c(0))){
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,1]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,1]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),data.individual[c(4:5),1])
        }else if ((data.individual[4,1] == c(0))&&(data.individual[5,1] == c(1))){
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,2]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,2]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),data.individual[c(4:5),1])
        }else if ((data.individual[4,1] == c(1))&&(data.individual[5,1] == c(0))){
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,3]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,3]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),data.individual[c(4:5),1])
        }else{
          Output = c(data.individual[2,1],rnorm(n = 1,mean = mu[2,4]+Sigma.matrix[1,2]/Sigma.matrix[1,1]*(data.individual[2,1]-mu[1,4]),
                                                sd = sqrt(Sigma.matrix[2,2]-Sigma.matrix[1,2]^2/Sigma.matrix[1,1])),data.individual[c(4:5),1])
        }
      }
    }
  }else{ # 3. Only discrete variables
    if ((c(4)%in% miss.index.individual)|(c(5)%in% miss.index.individual)){
      if ((c(4)%in% miss.index.individual)&&(c(5)%in% miss.index.individual)){ 
        # 3.1 Both discrete variables are missing ¡Ì
        Output = c(data.individual[c(2,3),1],nu[,sample(r.nu,size = 1,replace = TRUE,prob = Pi)])
      }else if ((c(4)%in% miss.index.individual)&&(!(c(5)%in% miss.index.individual))){ 
        # 3.2 The first discrete variable is missing ¡Ì
        if (data.individual[5,] == 1){
          Output = c(data.individual[c(2,3),1],sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[2,1]/(Pi[2,1]+Pi[4,1]),Pi[4,1]/(Pi[2,1]+Pi[4,1]))),data.individual[5,1])
        }else{
          Output = c(data.individual[c(2,3),1],sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[3,1]),Pi[3,1]/(Pi[1,1]+Pi[3,1]))),data.individual[5,1])
        }
      }else if((c(5)%in% miss.index.individual)&&(!(c(4)%in% miss.index.individual))){ 
        # 3.3 The second discrete variable is missing ¡Ì
        if (data.individual[4,] == 1){
          Output = c(data.individual[-c(1,5),1],sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[3,1]/(Pi[3,1]+Pi[4,1]),Pi[4,1]/(Pi[3,1]+Pi[4,1]))))
        }else{
          Output = c(data.individual[-c(1,5),1],sample(c(0,1),size = 1,replace = TRUE,prob = c(Pi[1,1]/(Pi[1,1]+Pi[2,1]),Pi[2,1]/(Pi[1,1]+Pi[2,1]))))
        }
      }
    }else{ # 4. Covariates fully observed ¡Ì
      Output = data.individual[-1,1]
    }
  }
  return(Output)
}