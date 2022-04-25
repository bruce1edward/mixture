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
beta_cur <- creg_oracle$coefficients
Breslow_Estimator1 <- basehaz(creg_oracle, centered=FALSE)
#edit(basehaz)

###############
mod <- survfit(creg_oracle, se.fit = FALSE)
Lambdahat_0_ <- Breslow_Estimator1$hazard
lambdahat_0_ <- diff(c(0,Lambdahat_0_))
#/(Breslow_Estimator1[,2] - c(0,Breslow_Estimator1[1:end-1,2]))
#abc <- cbind(c(0,Breslow_Estimator1[1:n-1,2]), Breslow_Estimator1[,2]);
#unique(y)

#plot(Breslow_Estimator1[,2],Breslow_Estimator1[,1], xlab = 'time', ylab = 'Cumulative Baseline Hazard', main = 'basehaz')
#lines(Breslow_Estimator1[,2],Breslow_Estimator1[,2]^4, xlab = 'time', ylab = 'Cumulative Baseline Hazard', main = 'basehaz')

#plot(Breslow_Estimator1[,2],lambdahat_0_, xlab = 'time', ylab = 'Baseline Hazard', main = 'basehaz')
#lines(Breslow_Estimator1[,2],4*Breslow_Estimator1[,2]^3, xlab = 'time', ylab = 'Baseline Hazard', main = 'basehaz')

#################definition###############################
#Breslow_Lambda <- function(surv_times, status, X, beta, t){
Breslow_Lambda <- function(surv_times, status, X, beta, cumulative){
  #0=right censored, 1=event at time #sum(status)
  #surv_times = y; status = 1 - cens; X = X; beta = creg_oracle$coefficients
  unique_death_times <- sort(unique(surv_times))
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
Breslow_eval_Lambda <- Breslow_Lambda(y, 1-cens, X, beta_cur, TRUE)
Breslow_eval_lambda <- Breslow_Lambda(y, 1-cens, X, beta_cur, FALSE)

##################Compared my implementation with basehaz function#########

plot(Breslow_eval_Lambda, Lambdahat_0_, xlab = "mine", ylab = "R based(basehaz)", main="Cumulative Baseline Hazard")



################################################################################################################
























##################Compared with true baseline harazrd######################

surv_times = y; status = 1 - cens; X = X; beta = beta 
unique_death_times <- sort(unique(surv_times[status == 1]))

plot(Breslow_eval_Lambda[match(y, unique_death_times)], y^4, xlab = "mine", ylab = "true")

x1 = Lambdahat_0_[match(y, Breslow_Estimator1$time)]
y1 = y^4

plot(x1[x1<20], y1[x1<20], xlab = "R based(basehaz)", ylab = "True", main="Cumulative Baseline Hazard")

cbind(x1[x1>20], y1[x1>20])


#-------------------------------------------------------------------------------------------

y1 = y[y<=3]
x1 = Lambdahat_0_[match(y1, Breslow_Estimator1$time)]
plot(x1, y1^4, xlab = "R based(basehaz)", ylab = "True", main="Cumulative Baseline Hazard")
abline(0,1)

y1 = y[y>3]
x1 = Lambdahat_0_[match(y1, Breslow_Estimator1$time)]
plot(x1, y1^4, xlab = "R based(basehaz)", ylab = "True", main="Cumulative Baseline Hazard")

y1 = y[y<=3]
x1 = lambdahat_0_[match(y1, Breslow_Estimator1$time)]
plot(x1, 4*y1^3, xlab = "R based(basehaz)", ylab = "True", main="Baseline Hazard")
abline(0,1)
