##################################################################################
### the ACT model with sigma unknown
##################################################################################
library(survival)
set.seed(132)
n <- 1000
d <- 10 #1
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
beta <- rnorm(d)
beta <-  beta/ sqrt(sum(beta^2))
mu <- X %*% beta
sigma <- 0.25
sigmastar <- sigma
xi_mult <- rweibull(n, shape = 1/sigma)
xi <- log(xi_mult)
y <- exp(drop(mu)) * xi_mult#
iteration <- 1
mse = array(0, c(iteration, 7, 6))
quat <- seq(0.1 , 0.5, by = 0.1)
cens <- (y >= 3)
###################Oracle#############
act_oracle <- survreg(Surv(y, event = 1 - cens) ~ X - 1, scale = sigma, dist = "weibull")
#################################################################################
for (iters in 1:iteration) {
  nmethod = 1
  for (alp in 1:length(quat)) {
  #alp = 1
  #####################################
  alpha <- quat[alp]
  alphastar <- alpha 
  k <- round(n * alpha)
  yperm <- y
  yperm[1:k] <- y[sample(1:k)]
  cens <- (yperm >= 3) # 11% censored #sum(cens)/n
  ###########################initial estimator (log normal) ################################## better
  lm_naive <- lm(log(yperm[!cens]) ~ X[!cens,] - 1, )
  betacur <- coef(lm_naive)
  mu <- X %*% betacur
  alpha <- 0.5
  sigma <- mean(residuals(lm_naive)^2)
  ###################Naive#############
  act_naive <- survreg(Surv(yperm, event = 1 - cens) ~ X - 1, scale = sigma, dist = "weibull")
  #####################################
  fymu <- function(mu, cens, sigma) dweibull(yperm, shape = 1/sigma, scale = exp(mu))^(1 - cens)  * (1 - pweibull(yperm, shape = 1/sigma, scale = exp(mu)))^cens
  dens_y  <- function(y, cens, sigma) rowMeans(sweep(dweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma))^(1-cens) * (1 - rowMeans(pweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma)))^cens
  fy <- function(sigma) dens_y(yperm, cens, sigma)
  nloglik_marginal <- function(mu, cens, sigma, alpha) sum(-log((1-alpha) * fymu(mu, cens, sigma) + alpha * fy))
  
  iter <- 1
  objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
  while(iter < maxiter){
    num <- (1-alpha) * fymu(mu, cens, sigma)
    denom <- num + alpha * fy(sigma)
    pcur <- num/denom
    alpha <- 1 - mean(pcur)
    lm_pcur <- survreg(Surv(yperm, event = 1 - cens) ~ X - 1, scale = sigma, dist = "weibull", weights = pmax(pcur, 1E-6))
    betacur <- coef(lm_pcur)
    mu <- X %*% betacur
    sigma <- lm_pcur$scale
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
    if(objs[iter] + tol > objs[iter-1])
      break
  }
  
  mse[iters,1,nmethod] = sqrt(sum((betacur - beta)^2))
  mse[iters,2,nmethod] = sqrt(sum((coef(act_naive) - beta)^2))
  mse[iters,3,nmethod] = sqrt(sum((coef(act_naive) - beta)^2))
  mse[iters,4,nmethod] = abs(sigmastar - sigma)
  mse[iters,5,nmethod] = abs(act_naive$scale - sigmastar)
  mse[iters,6,nmethod] = abs(act_oracle$scale - sigmastar)
  mse[iters,7,nmethod] = abs(alphastar - alpha)
  nmethod = nmethod + 1
  }
  iters = iters + 1
}
