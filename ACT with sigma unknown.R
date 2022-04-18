##################################################################################
### the ACT model with sigma unknown
##################################################################################
set.seed(132)
n <- 1000
d <- 10 #1
alpha <- 0.2
k <- round(n * alpha)
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
beta <- rnorm(d)
beta <-  beta/ sqrt(sum(beta^2))
mu <- X %*% beta
sigma <- 0.25
xi_mult <- rweibull(n, shape = 1/sigma)
xi <- log(xi_mult)
y <- exp(drop(mu)) * xi_mult#
yperm <- y
yperm[1:k] <- y[sample(1:k)]

############################initial estimators (log nomral) ################ better
lm_naive <- lm(log(yperm) ~ X)
betacur <- coef(lm_naive)[-1]
mu <- X %*% betacur
alpha <- 0.5
sigma <- mean(residuals(lm_naive)^2)

############################initial estimators (act model) ################
#act_naive <- survreg(Surv(yperm, event = rep(1,n)) ~ X - 1, dist = "weibull")
#betacur <- coef(act_naive)
#mu <- X %*% betacur
#alpha <- 0.5
#sigma <- act_naive$scale

# marginal distribution of y's is a scale mixture
gamma <- exp(drop(mu))
dens_y  <- function(y) rowMeans(sweep(dweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma))
fymu <- function(mu, sigma) dweibull(yperm, shape = 1/sigma, scale = exp(mu))
fy <- dens_y(yperm)
nloglik_marginal <- function(mu, sigma, alpha) sum(-log((1-alpha) * fymu(mu, sigma) + alpha * fy))
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)
iter <- 1
objs[iter] <- nloglik_marginal(mu, sigma, alpha)
# marginal approach
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, sigma)
  denom <- num + alpha * fy
  pcur <- num/denom
  alpha <- 1 - mean(pcur)
  lm_pcur <- survreg(Surv(yperm, event = rep(1,n)) ~ X - 1, scale = sigma, dist = "weibull", weights = pmax(pcur, 1E-6))
  betacur <- coef(lm_pcur)
  mu <- X %*% betacur
  sigma <- lm_pcur$scale
  iter <- iter + 1
  objs[iter] <- nloglik_marginal(mu, sigma, alpha)
  if(objs[iter] + tol > objs[iter-1])
    break
}
#sqrt(sum((beta - coef(lm_naive)[-1])^2))
abs(alpha - 0.2)
abs(sigma - 0.25)
sqrt(sum((beta - betacur)^2))

### now, additionally add censoring
cens <- (yperm >= 3) # 11% censored
fymu <- function(mu, cens, sigma) dweibull(yperm, shape = 1/sigma, scale = exp(mu))^(1 - cens)  * (1 - pweibull(yperm, shape = 1/sigma, scale = exp(mu)))^cens
dens_y  <- function(y, cens) rowMeans(sweep(dweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma))^(1-cens) * (1 - rowMeans(pweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma)))^cens
fy <- dens_y(yperm, cens)
nloglik_marginal <- function(mu, cens, sigma, alpha) sum(-log((1-alpha) * fymu(mu, cens, sigma) + alpha * fy))
###########################initial estimator (log normal) ################################## better
lm_naive <- lm(log(yperm[!cens]) ~ X[!cens,] - 1, )
betacur <- coef(lm_naive)
mu <- X %*% betacur
alpha <- 0.5
sigma <- mean(residuals(lm_naive)^2)

###########################initial estimators (act model) ################################
#act_naive <- survreg(Surv(yperm, event = 1 - cens) ~ X - 1, dist = "weibull")
#betacur <- coef(act_naive)
#mu <- X %*% betacur
#alpha <- 0.5
#sigma <- act_naive$scale

iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
while(iter < maxiter){
  num <- (1-alpha) * fymu(mu, cens, sigma)
  denom <- num + alpha * fy
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

abs(alpha - 0.2)
abs(sigma - 0.25)
sqrt(sum((betacur - beta)^2))