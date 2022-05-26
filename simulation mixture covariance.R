library(mvtnorm)
library(MASS)
set.seed(132)
setwd("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/simulations")
###Covariance Estimation
n <- 1000
d <- 3
k_n <- seq(0.1, 0.6, by = 0.1) 
iteration = 50
mse = array(0, c(iteration, 5, 6))
rho = 0.1

for (iters in 1:iteration) {
  nmethod = 1
  for (alph in k_n){
    #alph = 0.3
    k <- round(n * alph)
    var_covar <- matrix(nrow = d, ncol = d, rho)
    diag(var_covar) <- 1
    Z <- mvrnorm(n, mu = rep(0,d), Sigma = var_covar)
    X <- Z[,1]
    Y <- Z[,2:3]
    Xperm <- X
    Xperm[1:k] <- X[sample(1:k)]
    Zperm <-cbind(Xperm,Y)
    #naive
    cov_naive <- cov(Zperm, use="all.obs")
    #oracle
    cov_oracle <- cov(Z, use="all.obs")
    maxiter <- 10000
    alphacur <- 0.5
    tol <- 1E-4
    objs <- numeric(maxiter)
    iter <- 1
    ######
    f_xy <- function(omega) dmvnorm(Zperm, mean = rep(0, d), sigma = solve(omega), log = FALSE, checkSymmetry = TRUE)
    nloglik <- function(omega_sigma, omega_gamma, alpha) sum(-log((1-alpha) * f_xy(omega_sigma) + alpha * f_xy(omega_gamma)))
    S <- rbind(cbind(t(Xperm)%*%Xperm, t(Xperm)%*%Y), cbind(t(Y)%*%Xperm, t(Y)%*%Y))
    S_prime <- rbind(cbind(t(Xperm)%*%Xperm, matrix(nrow = 1, ncol = 2, 0)), cbind(matrix(nrow = 2, ncol = 1, 0), t(Y)%*%Y))
    #####initialization###############
    omega <- solve(cov_naive)
    omegacur <- solve(cov_naive)
    omega_p_cur <- omegacur
    omega_p_cur[1,2:3] <- matrix(nrow = 1, ncol = 2, 0)
    omega_p_cur[2:3,1] <- matrix(nrow = 2, ncol = 1, 0)
    objs[iter] <- nloglik(omegacur, omega_p_cur, alphacur)
    while(iter < maxiter){
      num <- alphacur * f_xy(omega_p_cur)
      denom <- num + (1 - alphacur) * f_xy(omegacur)
      pcur <- num/denom
      alphacur <- mean(pcur)
      w_sigma_covarfit <- cov.wt(Zperm, wt = 1 - pcur, cor = FALSE, center = FALSE,
                                 method = c("ML"))
      omegacur <- solve(w_sigma_covarfit$cov)
      w_gamma_covarfit <- cov.wt(Zperm, wt = pcur, cor = FALSE, center = FALSE,
                                 method = c("ML"))
      w_gamma_covarfit$cov[1,2:3]<- matrix(nrow = 1, ncol = 2, 0)
      w_gamma_covarfit$cov[2:3,1]<- matrix(nrow = 2, ncol = 1, 0)
      omega_p_cur <- solve(w_gamma_covarfit$cov)
      iter <- iter + 1
      objs[iter] <- nloglik(omegacur, omega_p_cur, alphacur)
      if(objs[iter] + tol > objs[iter-1])
        break
    }
    ###Robust-GLM method
    covariance_robust <- cov.rob(Zperm)
    
    ###LL and Chamber 
    #Q
    Q = (1 - alph - alph/(n-1)) * diag(n) + (alph/(n-1)) * matrix(1, n, n)
    Zperm_Q <-cbind(Q%*%Xperm,Y)
    cov_LL <- cov(Zperm_Q, use="all.obs")
    cov_LL[1,1] <- cov_naive[1,1]
    ####Naive, Mixture, Oracle
    mse[iters,1,nmethod] = norm(cov_naive  - var_covar, type = c("F"))
    mse[iters,2,nmethod] = norm(solve(omegacur) - var_covar, type = c("F"))
    mse[iters,3,nmethod] = norm(covariance_robust$cov - var_covar, type = c("F"))
    mse[iters,4,nmethod] = norm(cov_LL - var_covar, type = c("F"))
    mse[iters,5,nmethod] = norm(cov_oracle  - var_covar, type = c("F"))
    nmethod = nmethod + 1
  }
  iters = iters + 1
}

mse_m = colMeans(mse, dims = 1)

library(R.matlab)
filename <- paste("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/covariance", ".mat", sep = "")
writeMat(filename, Y = mse_m, rho = rho)