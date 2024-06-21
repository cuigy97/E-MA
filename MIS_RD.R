MIS_RD = function(data.full,            # Observations of each individual
                 mu,                    # Initial value of mu in the current iteration
                 Sigma.matrix,          # Initial value of covariance in the current iteration
                 beta = rep(1000,6),    # The default value indicates that the MIS only targets the covariates
                 sigma = 1000,          # Variance
                 model_index            # Indicate on which candidate model the resampling is based
){
  
  data.individual = as.matrix(data.full)
  miss.index.individual = c(which(is.na(data.individual)))
  samples_full = matrix(0,nrow = 50, ncol = 7)
  
  
  # Find the location of the current candidate model in CandidateModel.index.list
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
  CandidateModel.index.list = list()
  for(i in 1:5){
    CandidateModel.index.list[[i]] = combn(2:6, i)
  }
  index = c(1,CandidateModel.index.list[[num_covariate]][,model_cur_loc])
  index_com = setdiff(c(1:7),c(1,c(index)+1))
  
  # Resample the covariates
  data.tem = MIS_covariate_RD(data.covariate = data.individual[-1,],
                              mu = mu,
                              Sigma.matrix = Sigma.matrix)
  
  # Adjust the value of beta based on the selected candidate model
  Output_tem = cbind(matrix(0,ncol = 1,nrow = 50),data.tem)
  Output_tem[,index_com] = matrix(0,ncol = length(index_com), nrow = 50)
  
  
  if ((beta[1] == 1000)&&(sigma == 1000)){
    # Resample only the covariates
    samples_full = as.matrix(cbind(rep(data.individual[1,1],50),data.tem))
  }else{ 
    # Resampling all variables
    if (c(1)%in% miss.index.individual){ 
      # Response missing
      y.tem = rnorm(n = 50, mean = t(Output_tem[c(2:7)])%*%beta, sd = sqrt(sigma))
      samples_full = cbind(y.tem,data.tem)
    }else{ 
      # Response observed
      samples_full = cbind(matrix(rep(data.individual[1,1],50),ncol = 1),data.tem)
    }
  }

  return(list(samples_full))
}


