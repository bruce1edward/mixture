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


creg_naive <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1)
mu <- X %*% creg_naive$coef

###################### smoothed #################

lambdahat_0_g <- 4*yperm^3#lambdahat(yperm, 0.1)
Lambdahat_0_g <- yperm^4#Lambdahat(yperm, 0.1)

################## initialize baseline function #############
lambdahat_0 <- lambdahat_0_g
Lambdahat_0 <- Lambdahat_0_g

fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *  exp(-exp(mu) * L0)
fy <- function(l0, L0) l0^(1-cens) * exp(-L0)^cens
nloglik_marginal <- function(mu, cens, alpha, l0, L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy(lambdahat_0_g, Lambdahat_0_g)), na.rm = TRUE)
# update weights and alpha; then repeat
iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
beta_cur <- creg_naive$coef
# marginal approach
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, cens, lambdahat_0, Lambdahat_0)
  denom <- num + alpha * fy(lambdahat_0_g, Lambdahat_0_g)
  pcur <- num/denom
  alpha <- 1 - mean(pcur)
  #################partially likelihood for estimating beta#############################
  creg <- coxph(Surv(yperm, event = 1 - cens) ~ X - 1, weights = as.vector(pcur))
  #################Breslow Estimator of the Cumulative Baseline Hazard Rate#############
  Breslow_Estimator <- survfit(creg, type="aalen")
  Lambdahat_0 <- Breslow_Estimator$cumhaz[match(yperm, Breslow_Estimator$time)]
  lambdahat_0 <- c(Lambdahat_0[1], diff(Lambdahat_0))
  mu <- X %*% creg$coefficients
  beta_cur <- creg$coefficients
  iter <- iter + 1
  objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
  
  #####check#################
  #(exp(mu) * lambdahat_0)^(1 - cens) *  exp(-exp(mu) * Lambdahat_0)
  #length(mu); length(lambdahat_0)
  if(objs[iter] + tol > objs[iter-1])
    break
}
abs(alpha - 0.2)
sqrt(sum(-beta_cur*sigma - beta)^2)