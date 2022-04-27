##################################################################################
### cox proportional hazard model with baseline known (marginal g(y) is known)
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

lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

########################Naive###################################
survf0 <- survfit(Surv(yperm, event = 1 - cens) ~ 1)
Delta <- c(survf0$cumhaz[1], diff(survf0$cumhaz))
times <- survf0$time
##sum(sort(yperm) == times)
################## initialize baseline function #############
#################baseline for corrected matched data#################
lambdahat_0 <- 4*yperm^3
Lambdahat_0 <- yperm^4
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
lambdahat_0_ <- lambdahat_0
Lambdahat_0_ <- Lambdahat_0
# marginal approach
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, cens, lambdahat_0_, Lambdahat_0_)
  denom <- num + alpha * fy
  pcur <- num/denom
  alpha <- 1 - mean(pcur)
  #################partially likelihood for estimating beta#############################
  creg <- coxph(Surv(yperm, event = 1 - cens) ~ X, weights = pmax(drop(pcur),1E-6))
  mu <- X %*% creg$coefficients
  beta_cur <- creg$coefficients
  iter <- iter + 1
  objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0_, Lambdahat_0_)
  step_size = 0.1
  if(objs[iter] + tol > objs[iter-1]){
    break
    }
}

objs[1:100]
objs[iter]

abs(alpha - 0.2)
sqrt(sum(-beta_cur*sigma - beta)^2)

plot(-log(gamma)/sigma, mu)
abline(0,1, col = "blue", lwd = 3)

######################################################################################
survf0 <- survfit(Surv(yperm, event = 1 - cens) ~ 1)
Delta <- c(survf0$cumhaz[1], diff(survf0$cumhaz))
times <- survf0$time

# kernel estimation of the baseline hazard function
lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(dnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}
# corresponding cumulative hazard function
Lambdahat <- function(t, h, Delta){
  
  rowSums(sweep(pnorm(outer(t, times, "-"), sd = h), MARGIN = 2, FUN = "*", STATS = Delta))
  
}

tgrid2 <- seq(from = 0, to = 3, by = 0.01)
#eval_Lambdahat <- Lambdahat(tgrid2, h, Delta)
#eval_Lambdastar <- Lambdastar(tgrid2)

#lines(tgrid2, eval_Lambdastar, type = "l", col = "blue", lwd = 2)

#plot(Lambdastar(tgrid), Lambdahat(tgrid, 0.1), cex = 0.2)
#abline(0,1, lwd = 2, col = "blue")

plot(tgrid2, Lambdastar(tgrid2),main = "kernel estimation of the baseline cumulative hazard function vs True")
lines(tgrid2, Lambdahat(tgrid2, 0.1, Delta), col = "blue", lwd = 4, pch = 17)
lines(tgrid2, Lambdahat(tgrid2, 0.05, Delta), col = "red", lwd = 4, pch = 18)
lines(tgrid2, Lambdahat(tgrid2, 0.025, Delta), col = "green", lwd = 4, pch = 19)
legend(0, 1.5, c("h = 0.1", "h = 0.05", "h = 0.025"), col = c("blue", "red", "green"), pch = c(17,18,19))
