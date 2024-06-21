MIS_covariate_RD <- function(data.covariate,   # Observations of each individual, dimension: 6 times 1
                             mu,               # Initial value of mu in the current iteration
                             Sigma.matrix      # Initial value of covariance matrix in the current iteration
){
  
  data.individual = as.matrix(data.covariate)
  miss.index.individual = c(which(is.na(data.individual)))
  covariate_full = matrix(0,nrow = 50,ncol = 6)
  
  if ((c(1)%in% miss.index.individual)&&(c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)){
    
    # Covariate 1, 2, 3 are missing
    mu_tem = (cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual))[[1]]
    Sigma.matrix_tem = (cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual))[[2]]
    mu_1_tem = mu_tem[c(1:3)]
    mu_2_tem = mu_tem[c(4:6)]
    Sigma.matrix11_tem = Sigma.matrix_tem[c(1:3),c(1:3)]
    Sigma.matrix12_tem = Sigma.matrix_tem[c(1:3),c(4:6)]
    Sigma.matrix22_tem = Sigma.matrix_tem[c(4:6),c(4:6)]
    covariate_obs = data.individual[-miss.index.individual,]
    mu_con = mu_1_tem + Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(covariate_obs - mu_2_tem)
    Sigma.matrix_con = Sigma.matrix11_tem - Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%t(Sigma.matrix12_tem)
    covariate_miss = mvrnorm(n = 50,mu = mu_con,Sigma = Sigma.matrix_con)
    covariate_full_adjust = cbind(covariate_miss,matrix(rep(covariate_obs,50),nrow = 50, byrow = TRUE))
    covariate_full = covariate_full_adjust
    
  }else if ((c(1)%in% miss.index.individual)&&(c(2)%in% miss.index.individual)&&(!(c(3)%in% miss.index.individual))){
    
    # Covariate 1, 2 are missing
    mu_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[1]]
    Sigma.matrix_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[2]]
    mu_1_tem = mu_tem[c(1:2)]
    mu_2_tem = mu_tem[c(3:6)]
    Sigma.matrix11_tem = Sigma.matrix_tem[c(1:2),c(1:2)]
    Sigma.matrix12_tem = Sigma.matrix_tem[c(1:2),c(3:6)]
    Sigma.matrix22_tem = Sigma.matrix_tem[c(3:6),c(3:6)]
    covariate_obs = data.individual[-miss.index.individual,]
    mu_con = mu_1_tem + Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(covariate_obs - mu_2_tem)
    Sigma.matrix_con = Sigma.matrix11_tem - Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%t(Sigma.matrix12_tem)
    covariate_miss = mvrnorm(n = 50,mu = mu_con,Sigma = Sigma.matrix_con)
    covariate_full_adjust = cbind(covariate_miss,matrix(rep(covariate_obs,50),nrow = 50, byrow = TRUE))
    covariate_full = covariate_full_adjust
    
  }else if ((c(1)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)&&(!(c(2)%in% miss.index.individual))){
    
    # Covariate 1, 3 are missing
    mu_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[1]]
    Sigma.matrix_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[2]]
    mu_1_tem = mu_tem[c(1:2)]
    mu_2_tem = mu_tem[c(3:6)]
    Sigma.matrix11_tem = Sigma.matrix_tem[c(1:2),c(1:2)]
    Sigma.matrix12_tem = Sigma.matrix_tem[c(1:2),c(3:6)]
    Sigma.matrix22_tem = Sigma.matrix_tem[c(3:6),c(3:6)]
    covariate_obs = data.individual[-miss.index.individual,]
    mu_con = mu_1_tem + Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(covariate_obs - mu_2_tem)
    Sigma.matrix_con = Sigma.matrix11_tem - Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%t(Sigma.matrix12_tem)
    covariate_miss = mvrnorm(n = 50,mu = mu_con,Sigma = Sigma.matrix_con)
    covariate_full_adjust = cbind(covariate_miss,matrix(rep(covariate_obs,50),nrow = 50, byrow = TRUE))
    covariate_full = covariate_full_adjust[,c(1,3,2,4,5,6)]
    
  }else if ((c(2)%in% miss.index.individual)&&(c(3)%in% miss.index.individual)&&(!(c(1)%in% miss.index.individual))){
    
    # Covariate 2, 3 are missing
    mu_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[1]]
    Sigma.matrix_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[2]]
    mu_1_tem = mu_tem[c(1:2)]
    mu_2_tem = mu_tem[c(3:6)]
    Sigma.matrix11_tem = Sigma.matrix_tem[c(1:2),c(1:2)]
    Sigma.matrix12_tem = Sigma.matrix_tem[c(1:2),c(3:6)]
    Sigma.matrix22_tem = Sigma.matrix_tem[c(3:6),c(3:6)]
    covariate_obs = data.individual[-miss.index.individual,]
    mu_con = mu_1_tem + Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(covariate_obs - mu_2_tem)
    Sigma.matrix_con = Sigma.matrix11_tem - Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%t(Sigma.matrix12_tem)
    covariate_miss = mvrnorm(n = 50,mu = mu_con,Sigma = Sigma.matrix_con)
    covariate_full_adjust = cbind(covariate_miss,matrix(rep(covariate_obs,50),nrow = 50, byrow = TRUE))
    covariate_full = covariate_full_adjust[,c(3,1,2,4,5,6)]
    
  }else if ((c(1)%in% miss.index.individual)&&(!(c(2)%in% miss.index.individual))&&(!(c(3)%in% miss.index.individual))){
    
    # Covariate 1 is missing
    mu_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[1]]
    Sigma.matrix_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[2]]
    mu_1_tem = mu_tem[c(1)]
    mu_2_tem = mu_tem[c(2:6)]
    Sigma.matrix11_tem = Sigma.matrix_tem[c(1),c(1)]
    Sigma.matrix12_tem = Sigma.matrix_tem[c(1),c(2:6)]
    Sigma.matrix22_tem = Sigma.matrix_tem[c(2:6),c(2:6)]
    covariate_obs = data.individual[-miss.index.individual,]
    mu_con = mu_1_tem + Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(covariate_obs - mu_2_tem)
    Sigma.matrix_con = Sigma.matrix11_tem - Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(Sigma.matrix12_tem)
    covariate_miss = mvrnorm(n = 50,mu = mu_con,Sigma = Sigma.matrix_con)
    covariate_full_adjust = cbind(covariate_miss,matrix(rep(covariate_obs,50),nrow = 50, byrow = TRUE))
    covariate_full = covariate_full_adjust
    
  }else if ((c(2)%in% miss.index.individual)&&(!(c(1)%in% miss.index.individual))&&(!(c(3)%in% miss.index.individual))){
    
    # Covariate 2 is missing
    mu_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[1]]
    Sigma.matrix_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[2]]
    mu_1_tem = mu_tem[c(1)]
    mu_2_tem = mu_tem[c(2:6)]
    Sigma.matrix11_tem = Sigma.matrix_tem[c(1),c(1)]
    Sigma.matrix12_tem = Sigma.matrix_tem[c(1),c(2:6)]
    Sigma.matrix22_tem = Sigma.matrix_tem[c(2:6),c(2:6)]
    covariate_obs = data.individual[-miss.index.individual,]
    mu_con = mu_1_tem + Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(covariate_obs - mu_2_tem)
    Sigma.matrix_con = Sigma.matrix11_tem - Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(Sigma.matrix12_tem)
    covariate_miss = mvrnorm(n = 50,mu = mu_con,Sigma = Sigma.matrix_con)
    covariate_full_adjust = cbind(covariate_miss,matrix(rep(covariate_obs,50),nrow = 50, byrow = TRUE))
    covariate_full = covariate_full_adjust[,c(2,1,3,4,5,6)]
    
  }else if ((c(3)%in% miss.index.individual)&&(!(c(2)%in% miss.index.individual))&&(!(c(1)%in% miss.index.individual))){
    
    # Covariate 3 is missing
    mu_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[1]]
    Sigma.matrix_tem = cov_split(mu = mu, covariance = Sigma.matrix, missing.index = miss.index.individual)[[2]]
    mu_1_tem = mu_tem[c(1)]
    mu_2_tem = mu_tem[c(2:6)]
    Sigma.matrix11_tem = Sigma.matrix_tem[c(1),c(1)]
    Sigma.matrix12_tem = Sigma.matrix_tem[c(1),c(2:6)]
    Sigma.matrix22_tem = Sigma.matrix_tem[c(2:6),c(2:6)]
    covariate_obs = data.individual[-miss.index.individual,]
    mu_con = mu_1_tem + Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(covariate_obs - mu_2_tem)
    Sigma.matrix_con = Sigma.matrix11_tem - Sigma.matrix12_tem%*%solve(Sigma.matrix22_tem)%*%(Sigma.matrix12_tem)
    covariate_miss = mvrnorm(n = 50,mu = mu_con,Sigma = Sigma.matrix_con)
    covariate_full_adjust = cbind(covariate_miss,matrix(rep(covariate_obs,50),nrow = 50, byrow = TRUE))
    covariate_full = covariate_full_adjust[,c(2,3,1,4,5,6)]
    
  }
  return(covariate_full)
}
