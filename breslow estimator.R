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
cens <- (y >= 3) # 11% censored

########################True baseline hazard function###################
lambdahat_0 <- (1/sigma)*y^(1/sigma - 1)
Lambdahat_0 <- y^(1/sigma)

#################Breslow Estimator of the Cumulative Baseline Hazard Rate#############
creg_oracle <- coxph(Surv(y, event = 1 - cens) ~ X - 1)
Breslow_Estimator1 <- basehaz(creg_oracle, centered=FALSE)
Lambdahat_0_ <- Breslow_Estimator1$hazard
lambdahat_0_ <- diff(c(0,Lambdahat_0_))

#################definition###############################
#Breslow_Lambda <- function(surv_times, status, X, beta, t){
Breslow_Lambda <- function(cumulative){
  #0=right censored, 1=event at time #sum(status)
  surv_times = y; status = 1 - cens; X = X; beta = beta 
  unique_death_times <- sort(unique(surv_times[status == 1]))
  alpha <- length(unique_death_times)
  coxph_preds = X%*%beta
  for (i in seq_along(unique_death_times)) {
    alpha[i] <- sum(surv_times[status == 1] == unique_death_times[i])/sum(exp(coxph_preds[surv_times >= 
                                                                                           unique_death_times[i]]))
  }
  if (cumulative) 
    alpha <- cumsum(alpha)
  return(alpha)
}
Breslow_eval_Lambda <- Breslow_Lambda(TRUE)
Breslow_eval_lambda <- Breslow_Lambda(FALSE)

##################Compared my implementation with basehaz function#########

plot(Breslow_eval_Lambda, Lambdahat_0_[1:893], xlab = "mine", ylab = "R based")

##################Compared with true baseline harazrd######################

surv_times = y; status = 1 - cens; X = X; beta = beta 
unique_death_times <- sort(unique(surv_times[status == 1]))

plot(Breslow_eval_Lambda[match(y, unique_death_times)], y^4, xlab = "mine", ylab = "true")

plot(Lambdahat_0_[match(y, Breslow_Estimator1$time)], y^4, xlab = "R based", ylab = "True")
