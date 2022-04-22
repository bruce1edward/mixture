##################################################################################
### cox proportional hazard model with baseline unknown
##################################################################################
library(survival)
set.seed(132)
n <- 1000
d <- 10 #1
alpha <- 0.2
k <- round(n * alpha)
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
beta <- rnorm(d)
beta <-  beta/ sqrt(sum(beta^2))
mu <- X %*% beta
gamma <- exp(drop(mu))
sigma <- 0.25
xi_mult <- rweibull(n, shape = 1/sigma)
xi <- log(xi_mult)
y <- exp(drop(mu)) * xi_mult#
yperm <- y
yperm[1:k] <- y[sample(1:k)]
cens <- (yperm >= 3) # 11% censored
alpha <- 0.5
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)

########################Naive###################################
creg_naive <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1)
mu <- X %*% creg_naive$coef

# now estimate null component

survf0 <- survfit(Surv(yperm, event = 1 - cens) ~ 1)
Delta <- c(survf0$cumhaz[1], diff(survf0$cumhaz))
times <- survf0$time

# kernel estimation of the baseline hazard function
lambdahat <- function(t, h){
  
  rowSums(sweep(dnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

# corresponding cumulative hazard function
Lambdahat <- function(t, h){
  
  rowSums(sweep(pnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

fy <- lambdahat(yperm, 0.1)^(1-cens) * exp(-Lambdahat(yperm, 0.1))^cens

################## initialize baseline function #############
#################baseline for corrected matched data#################
lambdahat_0 <- 4*yperm^3
Lambdahat_0 <- yperm^4
#################baseline for  mismatched data################
g_lambdahat_0 <- lambdastar(yperm)# lambdahat(yperm, 0.1) #lambdastar(yperm)
g_Lambdahat_0 <- Lambdastar(yperm)# Lambdahat(yperm, 0.1) #Lambdahat(yperm)

fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *  exp(-exp(mu) * L0)
fy <- function(g_l0, g_L0) g_l0^(1-cens) * exp(-g_L0)^cens
nloglik_marginal <- function(mu, cens, alpha, l0, L0, g_l0, g_L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy(g_l0, g_L0)), na.rm = TRUE)
# update weights and alpha; then repeat
iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0, g_lambdahat_0, g_Lambdahat_0)
beta_cur <- creg_naive$coef
# marginal approach
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, cens, lambdahat_0, Lambdahat_0)
  denom <- num + alpha * fy(g_lambdahat_0, g_Lambdahat_0)
  pcur <- num/denom
  #####sum(is.na(pcur)); pcur[pcur != 0]; sum(pcur)
  ###############Since coxph can't handle 0 weight, so assign 10^-5 to the observation with 0 weight
  pcur[pcur == 0] = 10^-5 #sum(pcur == 0)
  alpha <- 1 - mean(pcur)
  #################partially likelihood for estimating beta#############################
  #creg <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1, weights = pcur, subset = (pcur != 0))
  #creg <- coxph(Surv(yperm[pcur != 0], event = 1 - cens[pcur != 0]) ~ X[pcur != 0,] - 1, weights = pcur[pcur != 0])
  creg <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1, weights = as.vector(pcur))
  #################Breslow Estimator of the Cumulative Baseline Hazard Rate#############
  #Breslow_Estimator <- survfit(creg, type="aalen")
  Breslow_Estimator1 <- basehaz(creg, centered=FALSE)
  #Lambdahat_0 <- Breslow_Estimator$cumhaz[match(yperm, Breslow_Estimator$time)]
  Lambdahat_0_ <- Breslow_Estimator1$hazard
  #plot(Breslow_Estimator$time, Breslow_Estimator1$time)
  #plot(Lambdahat_0, Lambdahat_0_)
  lambdahat_0 <- diff(c(0,Lambdahat_0_))
  
  ##################Compared with true baseline harazrd######################
  plot(lambdahat_0[match(yperm, Breslow_Estimator1$time)], 4*yperm^3, xlab = "Breslow Estimator of baseline harzard function", ylab = "True baseline harzard function for weibull")
  plot(Lambdahat_0_[match(yperm, Breslow_Estimator1$time)], yperm^4, xlab = "Breslow Estimator of baseline cumulative harzard function", ylab = "True baseline cumulative harzard function for weibull")
  
  #plot(Lambdahat_0[Lambdahat_0 < 1], (yperm^4)[yperm<1], xlab = "Breslow Estimator of baseline cumulative harzard function", ylab = "True baseline cumulative harzard function for weibull")
  #cbind(Lambdahat_0, yperm^4)
  
  #plot(Breslow_Estimator$time,  Lambdahat_0)
  mu <- X %*% creg$coefficients
  beta_cur <- creg$coefficients
  iter <- iter + 1
  objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0_, g_lambdahat_0, g_Lambdahat_0)
  #####check#################
  #(exp(mu) * lambdahat_0)^(1 - cens) *  exp(-exp(mu) * Lambdahat_0)
  #length(mu); length(lambdahat_0)
  if(objs[iter] + tol > objs[iter-1])
    break
}
#objs
abs(alpha - 0.2)
sqrt(sum(-beta_cur*sigma - beta)^2)