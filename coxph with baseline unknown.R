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

### baseline hazard
cens <- (yperm >= 3) # 11% censored
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

alpha <- 0.5
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)

# compute naive estimator
?coxph
?Surv

creg <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1)
bhaz <- basehaz(creg)
mu <- X %*% creg$coef

######################kernel density estimation#################

lambdahat_0 <- lambdahat(yperm, 0.1)
Lambdahat_0 <- Lambdahat(yperm, 0.1)

###################### smoothed #################

lambdahat_0 <- 4*yperm^3#lambdahat(yperm, 0.1)
Lambdahat_0 <- yperm^4#Lambdahat(yperm, 0.1)

fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *  exp(-exp(mu) * L0)
fy <- lambdahat(yperm, 0.1)^(1-cens) * exp(-Lambdahat(yperm, 0.1))^cens
nloglik_marginal <- function(mu, cens, alpha, l0, L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy))
# update weights and alpha; then repeat
iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
beta_cur <- creg$coef
# marginal approach
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, cens, lambdahat_0, Lambdahat_0)
  denom <- num + alpha * fy
  pcur <- num/denom
  alpha <- 1 - mean(pcur)
  
  t1_grad <- -colSums(sweep(X, MARGIN = 1, FUN = "*", STATS = (1-cens)*pcur))
  obj_cox_star <- function(beta, X, cens) sum(pcur* ((1 - cens) * (-X%*%beta -log(lambdahat_0)) + exp(X %*% beta) * Lambdahat_0))
  grad_cox_star <- function(beta, X, cens)  t1_grad + t(X) %*% (exp(X %*% beta) * Lambdahat_0 * pcur) #colSums(sweep(X, MARGIN = 1, FUN = "*", STATS = exp(X %*% beta) * Lambdastar_t))
  
  #####################issue optim initial value is not finite###########################
  #obj_cox_star(beta_cur, X, cens)
  #(-X%*%beta -log(lambdahat_0))
  
  # estimate regression parameters
  res_star <- optim(par = beta_cur, fn = obj_cox_star, gr = grad_cox_star, method = "BFGS", X = X, cens = cens)
  mu <- X %*% res_star$par
  beta_cur <- res_star$par
  iter <- iter + 1
  objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
  if(objs[iter] + tol > objs[iter-1])
    break
}
abs(alpha - 0.2)
sqrt(sum(-beta_cur*sigma - beta)^2)
#naive
sqrt(sum(-creg$coef*sigma - beta)^2)

######################Breslow Estimator of the Cumulative Baseline Hazard Rate################################

# compare with "true" underlying baseline hazard
tgrid2 <- seq(from = 0, to = 20, by = 0.01)

lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))))
Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

# Breslow estimator with survfit function (does not work)

#survfit(creg,type="aalen")
#plot(survfit(creg,type="aalen"),fun="cumhaz")
#plot(survfit(creg,type="aalen"))

#library(gbm3)
#?baseline_hazard
#edit(baseline_hazard)

#####################################Baseline Hazard####################################
#
Breslow_lambda <- function(t){
  #0=right censored, 1=event at time
  surv_times = y; status = 1 - cens; X = X; beta = beta 
  unique_death_times <- sort(unique(surv_times[status == 1]))
  lambda <- 0
  di <- sum(surv_times[status == 1] == unique_death_times[which.min(abs(t - unique_death_times))])
  lambda <- lambda + di/sum(exp(X[surv_times >= t,]%*%beta))                                                                       
  return(lambda)
}
Breslow_eval_lambdahat <- lapply(tgrid2, Breslow_lambda)

#####################################Cumulative Baseline Hazard####################################

#Breslow_Lambda <- function(surv_times, status, X, beta, t){
Breslow_Lambda <- function(t){
  #0=right censored, 1=event at time
  surv_times = y; status = 1 - cens; X = X; beta = beta 
  unique_death_times <- sort(unique(surv_times[status == 1]))
  alpha <- length(unique_death_times)
  ti <- unique_death_times[unique_death_times <= t]
  Lambda <- 0
  for (i in ti) {
    #print(i);i = ti[1];i
    di <- sum(surv_times[status == 1] == i)
    Lambda <- Lambda + di/sum(exp(X[surv_times >= i,]%*%beta))                                                                       
  }
  return(Lambda)
}
Breslow_eval_Lambdahat <- lapply(tgrid2, Breslow_Lambda)

#################################Compared with true baseline function##########################(not good)
eval_lambdastar <- lambdastar(tgrid2)
eval_Lambdastar <- Lambdastar(tgrid2)

plot(tgrid2, eval_lambdastar, type = "l", col = "blue", lwd = 2)
lines(tgrid2, Breslow_eval_lambdahat, col = "red", lwd = 1.5)

plot(lambdastar(tgrid2), Breslow_eval_lambdahat, cex = 0.2)
abline(0,1, lwd = 2, col = "blue")

plot(tgrid2, eval_Lambdastar, type = "l", col = "blue", lwd = 2)
lines(tgrid2, Breslow_eval_Lambdahat, col = "red", lwd = 1.5)

plot(Lambdastar(tgrid2), Breslow_eval_Lambdahat, cex = 0.2)
abline(0,1, lwd = 2, col = "blue")


###############################true baseline ####################################
lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))))
Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )
