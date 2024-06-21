# The E-MA algorithm
This file contains code and the dataset that are necessary for reproducing experimental results of simulation studies and the real data analysis of the E-MA algorithm.

## Code for simulation
### EMA-Linear-Misspecified models.R 
This R file contains code to reproduce the simulation results of different comparison methods, including the E-MS, E-MA, EMA-EW, EM-F, MS-CC, MS-CO, MS-Im and MA-Im.

### MIS.R
An R package that performs Metropolis independence sampler to generate responses and covariates from the true model.

### MIS_covariate.R
An R package that performs Metropolis independence sampler to generate covariates from the true model.

### density function_full.R
An R package that calculates the value of the joint density function of the response and covariates under the true model f(x) or an assumed model g(x).

### density covariate fun.R
An R package that calculates the value of the density function of the covariates under the true model f(x) or an assumed model g(x).


## Code for real data analysis
### Real data analysis.R
This R file contains code to reproduce the results of the real data analysis. Alternative methods include the E-MS, E-MA, EMA-EW, EM-F, MS-CC, MS-Im and MA-Im. The reason why we exclude the MS-CO method is that the missing data mechanism of the dataset clearly violates the MCAR pattern, under which the complete-case study has a poor performance.

### MIS_RD.R
A similar version of MIS.R that adapts to the real data analysis.

### MIS_covariate_RD.R
A similar version of MIS_covariate.R that adapts to the real data analysis.

### AirQualityUCI.csv
The dataset that contains the air quality information, which can be downloaded from the Machine Learning Repository at the University of California Irvine (https://archive.ics.uci.edu/dataset/360/air+quality).
