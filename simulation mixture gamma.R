library(robustbase)
set.seed(132)
setwd("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/simulations")
source("LL_C_IRLS_Gam.R")
source("P_GLM.R")
#source("fit_mixture_GLM.R")
#source("data_generation.R")
###Gamma Regression
n <- 1000
d <- 10
k_n <- seq(0.3, 0.8, by = 0.1) 
b <- 1
b0 <- 4
shape <- 10000
iteration = 100
mse = array(0, c(iteration, 6, 6))

for (iters in 1:iteration) {
  nmethod = 1
  for (alph in k_n){
    #alph = 0.3
    k <- round(n * alph)
    X <- abs(matrix(rnorm(n*d), nrow = n, ncol = d))
    beta <- abs(rnorm(d))
    beta <- b*c(beta)/ sqrt(sum(c(beta)^2))
    beta <- c(b0,beta)
    mu <- 1/(cbind(rep(1,n),X)%*%beta)  # use cannocial link for convenience
    y <- rgamma(n, shape = shape, scale = mu/shape)
    dens <- density(y)
    yperm <- y
    yperm[1:k] <- y[sample(1:k)]
    ############################################
    fymu <- function(mu) dgamma(yperm, shape = shape, scale = mu/shape)
    fy <- approx(dens$x,dens$y,xout=yperm)$y
    nloglik <- function(mu, alpha) sum(-log((1-alpha) * fymu(mu) + alpha * fy))
    ############################################
    maxiter <- 1000
    tol <- 1E-4
    objs <- numeric(maxiter)
    ##############################################
    glm_oracle <- glm(y ~ X, family = Gamma)
    glm0 <- glm(yperm ~ X, family = Gamma)
    betacur <- coef(glm0)
    mu <- 1/(cbind(rep(1,n),X)%*%betacur) 
    alpha <- 0.5
    iter <- 1
    objs[iter] <- nloglik(mu, alpha)
    
    while(iter < maxiter){
      num <- (1-alpha) * fymu(mu)
      denom <- num + alpha * fy
      pcur <- num/denom
      alpha <- 1 - mean(pcur)
      wglmfit <- glm(yperm ~ X, family = Gamma, weights = pcur)
      betacur <- coef(wglmfit)
      mu <- 1/(cbind(rep(1,n),X) %*% betacur)
      iter <- iter + 1
      objs[iter] <- nloglik(mu, alpha)
      if(objs[iter] + tol > objs[iter-1])
        break
    }
    
    
    ###Robust-GLM method
    glm_robust <- glmrob(yperm ~ X, family = Gamma)
    
    ###Penalized-GLM method
    glm_penalized <- P_GLM_Gam(cbind(rep(1,n),X), yperm, 0.01, coef(glm0), 1)
    
    ###LL and Chamber 
    #Q
    Q = (1 - alph - alph/(n-1)) * diag(n) + (alph/(n-1)) * matrix(1, n, n)
    glm_LL <- LL_IRLS_Gamma(d, X, Q, yperm, coef(glm0))
    glm_Chamber <- Chamber_IRLS_Gamma(d, X, Q, yperm, coef(glm0))
    
    ####Naive, Mixture, Oracle
    mse[iters,1,nmethod] = sqrt(sum((coef(glm0) - beta)^2))/sqrt(sum((coef(glm_oracle) - beta)^2))
    mse[iters,2,nmethod] = sqrt(sum((betacur - beta)^2))/sqrt(sum((coef(glm_oracle) - beta)^2))
    mse[iters,3,nmethod] = sqrt(sum(glm_penalized - beta)^2)/sqrt(sum((coef(glm_oracle) - beta)^2))
    mse[iters,4,nmethod] = sqrt(sum((coef(glm_robust) - beta)^2))/sqrt(sum((coef(glm_oracle) - beta)^2))
    mse[iters,5,nmethod] = sqrt(sum((as.numeric(unlist(glm_LL[1])) - beta)^2))/sqrt(sum((coef(glm_oracle) - beta)^2))
    mse[iters,6,nmethod] = sqrt(sum((as.numeric(unlist(glm_Chamber[1])) - beta)^2))/sqrt(sum((coef(glm_oracle) - beta)^2))
    nmethod = nmethod + 1
  }
  iters = iters + 1
}

mse_m = colMeans(mse, dims = 1)

library(R.matlab)
filename <- paste("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/simulations/gamma2", ".mat", sep = "")
writeMat(filename, Y = mse_m, b = b, b0 = b0, nu = shape)
