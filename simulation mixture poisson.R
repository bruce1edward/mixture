library(robustbase)
setwd("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/simulations")
source("LL_C_Armijo_Poi.R")
source("P_GLM.R")
source("fit_mixture_GLM.R")
source("data_generation.R")
###
set.seed(456)
###Poisson regression
n <- 1000
d <- 10
b <- 2
b0 <- 0
k_n <- seq(0.1, 0.6, by = 0.1) 
iteration = 1
mse = array(0, c(iteration, 7, 6))

#tic("sleeping")
for (iters in 1:iteration) {
   nmethod = 1
for (alph in k_n){
  #generate data
  #alph = 0.6
  data = dat_gen(n,d,b,b0,alph)
  X = data$X
  yperm = data$yperm
  beta = data$beta
  beta_oracle = data$beta_oracle
  
  #naive
  glm_naive <- glm(yperm ~ X, family = "poisson")
  
  #Mixture Method
  glm_mixture <- fit_mixture(X, yperm, 0.5, control = list(init = "IRLS", family = "poisson"))
  
  ###Robust-GLM method
  glm_robust <- glmrob(yperm ~ X, family = poisson, method= "Mqle", control = glmrobMqle.control(tcc= 1.2))
  
  ###Penalized-GLM method
  glm_penalized <- P_GLM_Poi(X,yperm,0.1,coef(glm_naive))
  
  ###LL and Chamber 
  #Q
  Q = (1 - alph - alph/(n-1)) * diag(n) + (alph/(n-1)) * matrix(1, n, n)
  glm_LL <- LL_Armijo_Poi(d, X, Q, yperm, coef(glm_naive))
  glm_Chamber <- Chamber_Armijo_Poi(d, X, Q, yperm, coef(glm_naive))
  
  ####Naive, mixture, Oralce
  mse[iters,1,nmethod] = sqrt(sum((coef(glm_naive) - beta)^2))
  mse[iters,2,nmethod] = sqrt(sum((glm_mixture$betahat - beta)^2))
  mse[iters,3,nmethod] = sqrt(sum(glm_penalized - beta)^2)
  mse[iters,4,nmethod] = sqrt(sum((coef(glm_robust) - beta)^2))
  mse[iters,5,nmethod] = sqrt(sum((as.numeric(unlist(glm_LL[1])) - beta)^2))
  mse[iters,6,nmethod] = sqrt(sum((as.numeric(unlist(glm_Chamber[1])) - beta)^2))
  mse[iters,7,nmethod] = sqrt(sum((beta_oracle - beta)^2))
  nmethod = nmethod + 1
}
  iters = iters + 1
}

#toc()
mse_m = colMeans(mse, dims = 1)

library(R.matlab)
filename <- paste("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/simulations/poisson", ".mat", sep = "")
writeMat(filename, Y = mse_m)