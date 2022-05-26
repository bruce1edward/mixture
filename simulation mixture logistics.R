library(robustbase)
set.seed(132)
setwd("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/simulations")
source("LL_C_Armijo_Bin.R")
source("P_GLM.R")
###Logistic Regression
n <- 1000
d <- 1
k_n <- seq(0.3, 0.8, by = 0.1) 
b <- 3
iteration = 100
mse = array(0, c(iteration, 7, 6))

for (iters in 1:iteration) {
  nmethod = 1
  for (alph in k_n){
    #alph = 0.3
    k <- round(n * alph)
    X <- matrix(rnorm(n*d), nrow = n, ncol = d)
    beta <- rnorm(d)
    beta <- beta/ sqrt(sum(beta^2))
    beta <- beta*b
    mu <- plogis(X %*% beta)
    y <- rbinom(n, size = 1, prob = mu)
    yperm <- y
    yperm[1:k] <- y[sample(1:k)]
    #
    p <- mean(yperm)
    fymu <- function(mu) (mu^yperm) * ((1 - mu)^(1 - yperm))   
    fy <- (p^yperm) * (1-p)^(1-yperm)
    nloglik <- function(mu, alpha) sum(-log((1-alpha) * fymu(mu) + alpha * fy))
    
    #
    maxiter <- 1000
    tol <- 1E-4
    objs <- numeric(maxiter)
    #
    
    glm_oracle <- glm(y ~ X - 1, family = "binomial")
    glm0 <- glm(yperm ~ X - 1, family = "binomial")
    betacur <- coef(glm0)
    mu <- plogis(X %*% betacur)
    alpha <- 0.5
    iter <- 1
    objs[iter] <- nloglik(mu, alpha)
    
    while(iter < maxiter){
      
      num <- (1-alpha) * fymu(mu)
      denom <- num + alpha * fy
      pcur <- num/denom
      alpha <- 1 - mean(pcur)
      wglmfit <- glm(yperm ~ X - 1, family = binomial, weights = pcur)
      betacur <- coef(wglmfit)
      mu <- plogis(X %*% betacur)
      iter <- iter + 1
      objs[iter] <- nloglik(mu, alpha)
      if(objs[iter] + tol > objs[iter-1])
        break
    }
    
    ###Robust-GLM method
    glm_robust <- glmrob(yperm ~ X - 1, family=binomial, method="Mqle", control = glmrobMqle.control(tcc=3.5))
    
    ###Penalized-GLM method
    glm_penalized <- P_GLM_bin(X, yperm, 0.008, coef(glm0), 1)
    
    ###LL and Chamber 
    #Q
    Q = (1 - alph - alph/(n-1)) * diag(n) + (alph/(n-1)) * matrix(1, n, n)
    glm_LL <- LL_Armijo_bin(d, X, Q, yperm, coef(glm0), 1)
    glm_Chamber <- Chamber_Armijo_bin(d, X, Q, yperm, coef(glm0), 1)
    
    ####Naive, Mixture, Oracle
    mse[iters,1,nmethod] = sqrt(sum((coef(glm0) - beta)^2))
    mse[iters,2,nmethod] = sqrt(sum((betacur - beta)^2))
    mse[iters,3,nmethod] = sqrt(sum(glm_penalized - beta)^2)
    mse[iters,4,nmethod] = sqrt(sum((coef(glm_robust) - beta)^2))
    mse[iters,5,nmethod] = sqrt(sum((as.numeric(unlist(glm_LL[1])) - beta)^2))
    mse[iters,6,nmethod] = sqrt(sum((as.numeric(unlist(glm_Chamber[1])) - beta)^2))
    mse[iters,7,nmethod] = sqrt(sum((coef(glm_oracle) - beta)^2))
    nmethod = nmethod + 1
  }
  iters = iters + 1
}

#
invalid <- mse[,2,5] > 4 | mse[,2,6] > 4

mse_m = colMeans(mse[!invalid, ,], dims = 1)

library(R.matlab)
filename <- paste("C:/Users/zwang39/OneDrive - George Mason University - O365 Production/Paper(Mixture)/simulations/logistics", ".mat", sep = "")
writeMat(filename, Y = mse_m)



