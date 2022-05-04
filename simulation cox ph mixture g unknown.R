##################################################################################
### cox proportional hazard model with baseline unknown (marginal g(y) is unknown)
##################################################################################
library(survival)
set.seed(132)
n <- 1000
d <- 10 #1
quat <- seq(0.05 , 0.4, by = 0.05)
k <- round(n * alpha)
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
beta <- rnorm(d)
beta <-  beta/ sqrt(sum(beta^2))
mu <- X %*% beta
gamma <- exp(drop(mu))
sigma <- 0.25
sigmastar <- sigma
xi_mult <- rweibull(n, shape = 1/sigma)
xi <- log(xi_mult)
y <- exp(drop(mu)) * xi_mult#
alpha <- 0.5
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)
h <- 0.1 #0.025 #0.025
iteration = 1
mse = array(0, c(iteration, 4, length(quat)))

# kernel estimation of the baseline hazard function
lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(dnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}
# corresponding cumulative hazard function
Lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(pnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

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
    ###################Oracle#############
    act_oracle <- coxph(Surv(y, event = 1 - cens) ~ X - 1)
    ###################Naive#############
    act_naive <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1)
    ########################Naive###################################
    survf0 <- survfit(Surv(yperm, event = 1 - cens) ~ 1)
    Delta <- c(survf0$cumhaz[1], diff(survf0$cumhaz))
    times <- survf0$time
    ##sum(sort(yperm) == times)
    ################## initialize baseline function #############
    #################baseline for corrected matched data#################
    lambdahat_0 <- lambdahat(yperm, h, Delta)
    Lambdahat_0 <- Lambdahat(yperm, h, Delta)
    #################baseline for  mismatched data################
    g_lambdahat_0 <- lambdahat(yperm, h, Delta)
    g_Lambdahat_0 <- Lambdahat(yperm, h, Delta)
    
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
    track_lam <- matrix(0, nrow = maxiter, ncol = n)
    track_Lam <- matrix(0, nrow = maxiter, ncol = n)
    track_alpha[1] <- alpha
    track_beta[1,] <- beta_cur
    track_lam[1,] <- lambdahat_0
    track_Lam[1,] <- Lambdahat_0
    
    lambdahat_0_ <- lambdahat_0
    Lambdahat_0_ <- Lambdahat_0
    # marginal approach
    while(iter < maxiter){
      num <- (1-alpha) * fymu(mu, cens, lambdahat_0_, Lambdahat_0_)
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
      track_lam[iter+1,] <- lambdahat_0_
      track_Lam[iter+1,] <- Lambdahat_0_
      iter <- iter + 1
      objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0_, Lambdahat_0_)
      if(objs[iter] + tol > objs[iter-1]){
        for(step_size in seq(0.05, 0.5, by = 0.05)){
          ####################################step size selection############################
          beta_cur_cov = track_beta[iter-1,]*(1 - step_size) + (track_beta[iter-1,] - track_beta[iter-2,]) *step_size
          mu_cov <- X %*% beta_cur_cov
          ##sum(mu_cov == mu)
          alpha_cov <- track_alpha[iter-1]*(1 - step_size) + (track_alpha[iter-1] - track_alpha[iter-2])*step_size
          ##alpha_cov == alpha
          lambdahat_0_cov <- track_lam[iter-1,]*(1 - step_size) + (track_lam[iter-1,] - track_lam[iter-2,])*step_size
          ##sum(lambdahat_0_cov == lambdahat_0_)
          Lambdahat_0_cov <- track_Lam[iter-1,]*(1 - step_size) + (track_Lam[iter-1,] - track_Lam[iter-2,])*step_size
          ##########################
          objs[iter] <- nloglik_marginal(mu_cov, cens, alpha_cov, lambdahat_0_cov, Lambdahat_0_cov)
          mu <- mu_cov
          alpha <- alpha_cov
          lambdahat_0_ <- lambdahat_0_cov
          Lambdahat_0_ <- Lambdahat_0_cov
          track_alpha[iter] <- alpha
          track_beta[iter,] <- beta_cur_cov
          track_lam[iter,] <- lambdahat_0_
          track_Lam[iter,] <- Lambdahat_0_
          if (objs[iter] + tol < objs[iter-1]){
            break
          }
          #objs[1:12]
        }
        iter = 5
        if (objs[iter] + tol > objs[iter-1]){
          print(paste0("step size failed at iteration ", iter))
          break
        }
      }
    }
    #objs[1:10]
    #objs[iter-1]
    #objs[iter]
    #abs(track_alpha[iter-1] - 0.2)
    #sqrt(sum(-track_beta[iter-1,]*sigma - beta)^2)
    #abs(alpha_cov - 0.2)
    #sqrt(sum(-beta_cur_cov*sigma - beta)^2)
    #abs(alpha - 0.2)
    #sqrt(sum(-beta_cur*sigma - beta)^2)
    #plot(-log(gamma)/sigma, mu)
    #abline(0,1, col = "blue", lwd = 3)
    
    mse[iters,1,nmethod] = sqrt(sum((-track_beta[iter-1,]*sigma - beta)^2))
    mse[iters,2,nmethod] = sqrt(sum((-coef(act_naive)*sigma- beta)^2))
    mse[iters,3,nmethod] = sqrt(sum((-coef(act_oracle)*sigma - beta)^2))
    mse[iters,4,nmethod] = abs(alphastar - alpha)
    nmethod = nmethod + 1
  }
  #nmethod = nmethod + 1
}

mse_m = colMeans(mse, dims = 1)

