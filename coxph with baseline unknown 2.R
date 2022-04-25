##################################################################################
### cox proportional hazard model with baseline unknown (marginal g(y) is known)
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
h <- 0.1


# kernel estimation of the baseline hazard function
lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(dnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

# corresponding cumulative hazard function
Lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(pnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

########################Naive###################################
survf0 <- survfit(Surv(yperm, event = 1 - cens) ~ 1)
Delta <- c(survf0$cumhaz[1], diff(survf0$cumhaz))
times <- survf0$time
##sum(sort(yperm) == times)
################## initialize baseline function #############
#################baseline for corrected matched data#################
lambdahat_0 <- lambdahat(yperm, 0.1, Delta)
Lambdahat_0 <- Lambdahat(yperm, 0.1, Delta)
#################baseline for  mismatched data################
g_lambdahat_0 <- lambdastar(yperm)
g_Lambdahat_0 <- Lambdastar(yperm)

fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *  exp(-exp(mu) * L0)
fy <- g_lambdahat_0^(1-cens) * exp(-g_Lambdahat_0)^cens
nloglik_marginal <- function(mu, cens, alpha, l0, L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy), na.rm = TRUE)
# update weights and alpha; then repeat
iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
########################Naive###################################
creg_naive <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1)
mu <- X %*% creg_naive$coef
beta_cur <- creg_naive$coef
##########################
track_alpha <- rep(0, maxiter)
track_beta <- matrix(0, nrow = maxiter, ncol = d)
track_alpha[1] <- alpha
track_beta[1,] <- beta_cur
# marginal approach
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, cens, lambdahat_0, Lambdahat_0)
  denom <- num + alpha * fy
  pcur <- num/denom
  alpha <- 1 - mean(pcur)
  track_alpha[iter+1] <- alpha
  #################partially likelihood for estimating beta#############################
  creg <- coxph(Surv(yperm, event = 1 - cens) ~ X, weights = pmax(drop(pcur),1E-6))
  mu <- X %*% creg$coefficients
  beta_cur <- creg$coefficients
  track_beta[iter+1,] <- beta_cur
  #################Breslow Estimator of the Cumulative Baseline Hazard Rate#############
  Breslow_Estimator <- basehaz(creg, centered=FALSE)
  #edit(basehaz)
  cumhazard <- Breslow_Estimator$hazard
  Delta1 <- diff(c(0,cumhazard))
  times <- Breslow_Estimator$time
  #sum(Breslow_Estimator$time == times)
  lambdahat_0_ <- lambdahat(yperm, h, Delta1)
  Lambdahat_0_ <- Lambdahat(yperm, h, Delta1)
  iter <- iter + 1
  objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0_, Lambdahat_0_)
  #while(objs[iter] + tol > objs[iter-1]){
    ###line search
  #  h = 0.8*h
  #  lambdahat_0_ <- lambdahat(yperm, h, Delta1)
  #  Lambdahat_0_ <- Lambdahat(yperm, h, Delta1)
  #  objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0_, Lambdahat_0_)
  #  if(h < tol)
  #   break
  #}
  step_size = 0.9
  if(objs[iter] + tol > objs[iter-1]){
    break
    
    ####################################step size selection############################
    comment = function() {
      
      Breslow_Lambda <- function(surv_times, status, X, beta, cumulative, weight){
        #0=right censored, 1=event at time #sum(status)
        #surv_times = y; status = 1 - cens; X = X; beta = creg_oracle$coefficients
        unique_death_times <- sort(unique(surv_times))
        alpha <- length(unique_death_times)
        coxph_preds = (X%*%beta)*weight
        for (i in seq_along(unique_death_times)) {
          alpha[i] <- sum(surv_times[status == 1] == unique_death_times[i])/sum(exp(coxph_preds[surv_times >= 
                                                                                                  unique_death_times[i]]))
        }
        if (cumulative) 
          alpha <- cumsum(alpha)
        return(alpha)
      }
      
      beta_cur_cov = beta_cur*(1 - step_size) + track_beta[iter-1,]*step_size
      mu <- X %*% beta_cur_cov
      cumhazard1 <- Breslow_Lambda(yperm, 1-cens, X, beta_cur, TRUE, pmax(drop(pcur),1E-6))
      #plot(cumhazard, cumhazard1)
      Delta1 <- diff(c(0,cumhazard))
      times <- sort(unique(yperm))
      #sum(Breslow_Estimator$time == times)
      lambdahat_0_ <- lambdahat(yperm, h, Delta1)
      Lambdahat_0_ <- Lambdahat(yperm, h, Delta1)
      #objs[iter] <- 
      nloglik_marginal(mu, cens, alpha, lambdahat_0_, Lambdahat_0_)
      
    }
    
  }
}

abs(alpha - 0.2)
sqrt(sum(-beta_cur*sigma - beta)^2)

abs(track_alpha[iter-1] - 0.2)
sqrt(sum(-track_beta[iter-1,]*sigma - beta)^2)

##################################################################################
### cox proportional hazard model with baseline unknown (marginal g(y) is unknown)
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


# kernel estimation of the baseline hazard function
lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(dnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

# corresponding cumulative hazard function
Lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(pnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

########################Naive###################################
survf0 <- survfit(Surv(yperm, event = 1 - cens) ~ 1)
Delta <- c(survf0$cumhaz[1], diff(survf0$cumhaz))
times <- survf0$time
##sum(sort(yperm) == times)
################## initialize baseline function #############
#################baseline for corrected matched data#################
lambdahat_0 <- lambdahat(yperm, 0.1, Delta)
Lambdahat_0 <- Lambdahat(yperm, 0.1, Delta)
#################baseline for  mismatched data################
g_lambdahat_0 <- lambdahat(yperm, 0.1, Delta) #lambdastar(yperm)
g_Lambdahat_0 <- Lambdahat(yperm, 0.1, Delta) #Lambdahat(yperm)

fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *  exp(-exp(mu) * L0)
fy <- g_lambdahat_0^(1-cens) * exp(-g_Lambdahat_0)^cens
nloglik_marginal <- function(mu, cens, alpha, l0, L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy), na.rm = TRUE)
# update weights and alpha; then repeat
iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
########################Naive###################################
creg_naive <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1)
mu <- X %*% creg_naive$coef
beta_cur <- creg_naive$coef
##########################
track_alpha <- rep(0, maxiter)
track_beta <- matrix(0, nrow = maxiter, ncol = d)
track_alpha[1] <- alpha
track_beta[1,] <- beta_cur
h = 0.1
# marginal approach
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, cens, lambdahat_0, Lambdahat_0)
  denom <- num + alpha * fy
  pcur <- num/denom
  pcur[pcur == 0] = 10^-5 #sum(pcur == 0)
  alpha <- 1 - mean(pcur)
  track_alpha[iter+1] <- alpha
  #################partially likelihood for estimating beta#############################
  creg <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1, weights = as.vector(pcur))
  mu <- X %*% creg$coefficients
  beta_cur <- creg$coefficients
  track_beta[iter+1,] <- beta_cur
  #################Breslow Estimator of the Cumulative Baseline Hazard Rate#############
  Breslow_Estimator <- basehaz(creg, centered=FALSE)
  cumhazard <- Breslow_Estimator$hazard
  Delta1 <- diff(c(0,cumhazard))
  times <- cumhazard$time
  #sum(Breslow_Estimator$time == times)
  lambdahat_0_ <- lambdahat(yperm, h, Delta1)
  Lambdahat_0_ <- Lambdahat(yperm, h, Delta1)
  iter <- iter + 1
  objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0_, Lambdahat_0_)
  while(objs[iter] + tol > objs[iter-1]){
    ###line search
    h = 0.8*h
    lambdahat_0_ <- lambdahat(yperm, h, Delta1)
    Lambdahat_0_ <- Lambdahat(yperm, h, Delta1)
    objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0_, Lambdahat_0_)
    if(h < tol)
      break
  }
}

abs(alpha - 0.2)
sqrt(sum(-beta_cur*sigma - beta)^2)

abs(track_alpha[iter-1] - 0.2)
sqrt(sum(-track_beta[iter-1,]*sigma - beta)^2)