library(survival)
set.seed(132)
### Example 1: Laplace noise
dlognormal <- function(x, mean, sd){
   dnorm(log(x), mean = mean, sd = sd) / x
}
n <- 1000
d <- 10 #1
alpha <- 0.2
k <- round(n * alpha)
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
beta <- rnorm(d)
beta <-  beta/ sqrt(sum(beta^2))
mu <- X %*% beta
sigma <- 0.25
xi <- sigma * rnorm(n)  # rlaplace(n)
y <- exp(drop(mu + xi)) # log-normal response
yperm <- y
yperm[1:k] <- y[sample(1:k)]
fymu <- function(mu, sigma) dlognormal(yperm, mu, sigma)
# assume known for simplicity (location mixture)
dens_y  <- function(y) colMeans((sapply(y, function(z) dlognormal(z, mu, sigma))))
kde_obj <- density(yperm)
kde_y <- approxfun(kde_obj$x, kde_obj$y)
#ygrid <- seq(from = 0, to = 12, length = 1000)
#plot(ygrid, kde_y(ygrid), type = "l")
#lines(ygrid, dens_y(ygrid), type = "l", col = "red")
fy <- dens_y(yperm)
# plot likelihood ratios
plot(log(fymu(mu, sigma) / fy))
nloglik_marginal <- function(mu, sigma, alpha) sum(-log((1-alpha) * fymu(mu, sigma) + alpha * fy))
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)
# find some simple initial estimators
lm_naive <- lm(log(yperm) ~ X - 1)
betacur <- coef(lm_naive)
mu <- X %*% betacur
alpha <- 0.5
sigma <- mean(residuals(lm_naive)^2)
iter <- 1
objs[iter] <- nloglik_marginal(mu, sigma, alpha)
# marginal approach
while(iter < maxiter){
    num <- (1-alpha) * fymu(mu, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)
    lm_pcur  <- lm(log(yperm) ~ X - 1, weights = pcur)
    betacur <- coef(lm_pcur)
    mu <- X %*% betacur
    sigma <- sqrt(weighted.mean(residuals(lm_pcur)^2, w = pcur))
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, sigma, alpha)
    if(objs[iter] + tol > objs[iter-1])
        break
}
alpha
sigma
sqrt(sum((betacur - beta)^2))

############################################################################################################

# NEXT: account for possible right-censoring [assuming sigma^2 is known]
cens <- (yperm >= 3) 
fymu <- function(mu, cens, sigma) dlnorm(yperm, mu, sigma)^(1-cens) * (1 - plnorm(yperm, mu, sigma))^cens
# assume known for simplicity (location mixture)
dens_y  <- function(y, cens) colMeans((sapply(y, function(z) dlnorm(z, mu, sigma))))^(1- cens) * (colMeans((sapply(y, function(z) 1-plnorm(z, mu, sigma)))))^cens
fy <- dens_y(yperm, cens)
#plot(log(fymu(mu, cens, sigma) / fy))
#plot(log(fymu(mu, cens, sigma)[cens] / fy[cens]))
nloglik_marginal <- function(mu, cens, sigma, alpha) sum(-log((1-alpha) * fymu(mu, cens, sigma) + alpha * fy))
lm_naive <- lm(log(yperm[!cens]) ~ X[!cens,] - 1, )
betacur <- coef(lm_naive)
mu <- X %*% betacur
alpha <- 0.5
sigma <- .25 #ASSUME KNOWN
iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
# marginal approach
while(iter < maxiter){
    num <- (1-alpha) * fymu(mu, cens, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)
    lm_pcur <- survreg(Surv(yperm, event = 1 - cens) ~ X - 1, scale = 0.25, dist = "lognormal", weights = pmax(pcur, 1E-6))
    betacur <- coef(lm_pcur)
    mu <- X %*% betacur
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
    if(objs[iter] + tol > objs[iter-1])
        break
}
alpha
sqrt(sum((betacur - beta)^2))
##################################################################################
### Moving forward: the PH model
##################################################################################
ls()
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
#hist(xi, nclass = 100, prob = TRUE)

y <- exp(drop(mu)) * xi_mult#
yperm <- y
yperm[1:k] <- y[sample(1:k)]

# marginal distribution of y's is a scale mixture
gamma <- exp(drop(mu))
dens_y  <- function(y) rowMeans(sweep(dweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma))
#ygrid <- seq(from = 0, to = 12, length = 1000)
#hist(y, nclass = 50, prob = TRUE, ylim = c(0,1))
#lines(ygrid, dens_y(ygrid), type = "l")
#kde_obj <- density(yperm)
#kde_y <- approxfun(kde_obj$x, kde_obj$y)
#lines(ygrid, kde_y(ygrid), col = "red")
#integrate(dens_y, lower = 0, upper = 20)
# without censoring, regression parameters should be accessible from linear regression
lm_oracle <- lm(log(y) ~ X)
coef(lm_oracle)
#sqrt(sum((beta - coef(lm_oracle)[-1])^2))
# now, add mismatches [assuming scale is known]
fymu <- function(mu, sigma) dweibull(yperm, shape = 1/sigma, scale = exp(mu))
fy <- dens_y(yperm)
# plot likelihood ratios
#plot(log(fymu(mu, sigma) / fy), ylim = c(0,10))
nloglik_marginal <- function(mu, sigma, alpha) sum(-log((1-alpha) * fymu(mu, sigma) + alpha * fy))
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)

# find some simple initial estimators
lm_naive <- lm(log(yperm) ~ X)
betacur <- coef(lm_naive)[-1]
mu <- X %*% betacur
alpha <- 0.5
#sigma <- mean(residuals(lm_naive)^2)
iter <- 1
objs[iter] <- nloglik_marginal(mu, sigma, alpha)
# marginal approach
while(iter < maxiter){
    num <- (1-alpha) * fymu(mu, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)
    lm_pcur <- survreg(Surv(yperm, event = rep(1,n)) ~ X - 1, scale = sigma, dist = "weibull", weights = pmax(pcur, 1E-6))
    #lm_pcur  <- lm(log(yperm) ~ X - 1, weights = pcur)
    betacur <- coef(lm_pcur)
    mu <- X %*% betacur
    #sigma <- sqrt(weighted.mean(residuals(lm_pcur)^2, w = pcur))
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, sigma, alpha)
    if(objs[iter] + tol > objs[iter-1])
        break
}
sqrt(sum((beta - coef(lm_naive)[-1])^2))
sqrt(sum((beta - betacur)^2))

### now, additionally add censoring
cens <- (yperm >= 3) # 11% censored
fymu <- function(mu, cens, sigma) dweibull(yperm, shape = 1/sigma, scale = exp(mu))^(1 - cens)  * (1 - pweibull(yperm, shape = 1/sigma, scale = exp(mu)))^cens
# assume known for simplicity (scale mixture)
dens_y  <- function(y, cens) rowMeans(sweep(dweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma))^(1-cens) * (1 - rowMeans(pweibull(outer(y, gamma, FUN = "/"), shape = 1/sigma)))^cens
fy <- dens_y(yperm, cens)
#
plot(log(fymu(mu, cens, sigma) / fy), ylim = c(-10, 10), col = cens  +1)
#plot(log(fymu(mu, cens, sigma)[cens] / fy[cens]))


nloglik_marginal <- function(mu, cens, sigma, alpha) sum(-log((1-alpha) * fymu(mu, cens, sigma) + alpha * fy))
lm_naive <- lm(log(yperm[!cens]) ~ X[!cens,] - 1, )
betacur <- coef(lm_naive)
mu <- X %*% betacur
alpha <- 0.5
sigma <- .25 #ASSUME KNOWN
iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
#A <- M1(mu, sigma) * (alpha/(n-1))
#diag(A) <- diag(A) * (1 - alpha)/(alpha/(n-1))
#B <- fyM * (1 - alpha/(n-1))
#diag(B) <- fyM * alpha / (1 - alpha/(n-1))
#objs_pairwise[iter] <- nloglik_pairwise(A,B)

# marginal approach
while(iter < maxiter){

    num <- (1-alpha) * fymu(mu, cens, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)
    #deriv <- (dlnorm(yperm, mu, sigma)^2)/(1 - plnorm(yperm, mu, sigma)) * (-(log(yperm) - mu))/(sigma^2)
    #deriv_num <- (-log(1 - plnorm(yperm, mu + 0.001, sigma)) - -log(1 - plnorm(yperm, mu -0.001, sigma)))/(2*0.001)
  #  deriv_num <- (-log(dlnorm(yperm, mu + 0.00001, sigma))  - -log(dlnorm(yperm, mu -0.00001, sigma)))/(2 * 0.00001)
    #deriv <- (log(yperm) - mu)/(sigma^2)
    #cens_weights <- rep(1,n)
    #cw <- dlnorm(yperm, mu, sigma)^2/(1 - plnorm(yperm, mu, sigma))
    #cw[!is.finite(cw)] <- 0
    #cens_weights[cens] <- cw[cens]
    lm_pcur <- survreg(Surv(yperm, event = 1 - cens) ~ X - 1, scale = 0.25, dist = "weibull", weights = pmax(pcur, 1E-6))
    #lm_pcur  <- lm(log(yperm) ~ X - 1, weights = pcur * cens_weights)
    betacur <- coef(lm_pcur)
    mu <- X %*% betacur
    #sigma <- sqrt(weighted.mean(residuals(lm_pcur)^2, w = pcur))
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
    if(objs[iter] + tol > objs[iter-1])
        break
}

alpha
sqrt(sum((betacur - beta)^2))

### now start working on the COX PH model

### baseline hazard
ls()
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
tgrid <- seq(from = 0, to = 5, by = 0.01)
eval_lambdahat <- lambdahat(tgrid, 0.1)
plot(tgrid, eval_lambdahat, type = "l")
tgrid2 <- seq(from = 0, to = 20, by = 0.01)
eval_Lambdahat <- Lambdahat(tgrid2, 0.1)
plot(times, survf0$cumhaz)
lines(tgrid2, eval_Lambdahat, type = "l", col = "red", lwd = 2)

# compare with "true" underlying baseline hazard

lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )
Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )
eval_Lambdastar <- Lambdastar(tgrid2)
lines(tgrid2, eval_Lambdastar, type = "l", col = "blue", lwd = 2)
plot(Lambdastar(tgrid), Lambdahat(tgrid, 0.1), cex = 0.2)
abline(0,1, lwd = 2, col = "blue")
plot(tgrid2, lambdastar(tgrid2))
lines(tgrid2, lambdahat(tgrid2, 0.1), col = "blue", lwd = 1.5)

# estimation of Cox model assuming BL functions are known

Xcent <- scale(X, center = TRUE, scale = FALSE)
#Xcent <- X
lambdastar_t <- lambdastar(y)#4*y^3#lambdastar(y)
Lambdastar_t <- Lambdastar(y)#y^4#Lambdastar(y)
t1_grad <- -colSums(sweep(Xcent, MARGIN = 1, FUN = "*", STATS = (1-cens)))
obj_cox_star <- function(beta, X, cens) sum( (1 - cens) * (-X%*%beta -log(lambdastar_t)) + exp(X %*% beta) * Lambdastar_t)
grad_cox_star <- function(beta, X, cens)  t1_grad + t(X) %*% (exp(X %*% beta) * Lambdastar_t) #colSums(sweep(X, MARGIN = 1, FUN = "*", STATS = exp(X %*% beta) * Lambdastar_t))
creg <- coxph(Surv(y, event = 1 - cens) ~ Xcent)
bhaz <- basehaz(creg)
#plot(lambdastar_t * exp(Xcent %*%res_star$par),
# estimate regression parameters
res_star <- optim(par = rep(0,d), fn = obj_cox_star, gr = grad_cox_star, method = "BFGS", X = Xcent, cens = cens)

plot(exp(-log(gamma)/sigma) * dweibull(y, shape = 1/sigma)/(1-pweibull(y, shape = 1/sigma)))
plot(-log(gamma)/sigma, Xcent %*% res_star$par)

plot(-log(gamma)/sigma + log(4*y^3), Xcent %*% res_star$par)#%

plot(-log(gamma)/sigma, Xcent %*% res_star$par)
plot(-log(gamma)/sigma, Xcent %*% creg$coef)
times <- bhaz$time
Delta <- c(bhaz$hazard[1], diff(bhaz))

lambdahat_0 <- lambdastar(y)#4*y^3#lambdastar(y)
Lambdahat_t <- Lambdastar(y)

plot(-log(gamma)/sigma, Xcent %*% creg$coef)
abline(0,1, col = "blue", lwd = 3)
yy <- Xcent %*% res_star$par
xx <- -log(gamma)/sigma
lm(yy ~ xx)
#plot(c(bhaz[1,1], diff(bhaz[,1])) * predict(creg, type = "risk"), lambdastar_t * exp(Xcent %*%res_star$par))

### pseudo-code

fy <- dens_y(yperm, cens)
alpha <- 0.5
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)

# compute naive estimator

creg <- coxph(Surv(yperm, event = 1 - cens) ~ Xcent)
bhaz <- basehaz(creg)
mu <- Xcent %*% creg$coef

# compute associated smoothed baseline and cumulative baseline hazard

times <- bhaz$time
Delta <- c(bhaz$hazard[1], diff(bhaz))

lambdahat_0 <- lambdahat(yperm, 0.1)
Lambdahat_0 <- Lambdahat(yperm, 0.1)

fymu <- function(mu, cens, l0, L0) exp( (1 - cens) * (-X%*%beta -log(l0)) + exp(mu) * L0)
nloglik_marginal <- function(mu, cens, alpha, l0, L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy))
# update weights and alpha; then repeat

iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
#A <- M1(mu, sigma) * (alpha/(n-1))
#diag(A) <- diag(A) * (1 - alpha)/(alpha/(n-1))
#B <- fyM * (1 - alpha/(n-1))
#diag(B) <- fyM * alpha / (1 - alpha/(n-1))
#objs_pairwise[iter] <- nloglik_pairwise(A,B)

# marginal approach
while(iter < maxiter){
    num <- (1-alpha) * fymu(mu, cens, sigma)
    denom <- num + alpha * fy
    pcur <- num/denom
    alpha <- 1 - mean(pcur)
    #deriv <- (dlnorm(yperm, mu, sigma)^2)/(1 - plnorm(yperm, mu, sigma)) * (-(log(yperm) - mu))/(sigma^2)
    #deriv_num <- (-log(1 - plnorm(yperm, mu + 0.001, sigma)) - -log(1 - plnorm(yperm, mu -0.001, sigma)))/(2*0.001)
    #deriv_num <- (-log(dlnorm(yperm, mu + 0.00001, sigma))  - -log(dlnorm(yperm, mu -0.00001, sigma)))/(2 * 0.00001)
    #deriv <- (log(yperm) - mu)/(sigma^2)
    #cens_weights <- rep(1,n)
    #cw <- dlnorm(yperm, mu, sigma)^2/(1 - plnorm(yperm, mu, sigma))
    #cw[!is.finite(cw)] <- 0
    #cens_weights[cens] <- cw[cens]
    lm_pcur <- survreg(Surv(yperm, event = 1 - cens) ~ X - 1, scale = 0.25, dist = "weibull", weights = pmax(pcur, 1E-6))
    #lm_pcur  <- lm(log(yperm) ~ X - 1, weights = pcur * cens_weights)
    betacur <- coef(lm_pcur)
    mu <- X %*% betacur
    #sigma <- sqrt(weighted.mean(residuals(lm_pcur)^2, w = pcur))
    iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, cens, sigma, alpha)
    if(objs[iter] + tol > objs[iter-1])
        break
}


alpha
sqrt(sum((betacur - beta)^2))

### now start working on the COX PH model

### baseline hazard
ls()
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


tgrid <- seq(from = 0, to = 5, by = 0.01)
eval_lambdahat <- lambdahat(tgrid, 0.1)

plot(tgrid, eval_lambdahat, type = "l")

tgrid2 <- seq(from = 0, to = 20, by = 0.01)
eval_Lambdahat <- Lambdahat(tgrid2, 0.1)
plot(times, survf0$cumhaz)
lines(tgrid2, eval_Lambdahat, type = "l", col = "red", lwd = 2)

# compare with "true" underlying baseline hazard

lambdastar <- function(t) rowMeans(sweep(dweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma), MARGIN = 2, FUN = "/", STATS = gamma)) /( (1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )

Lambdastar <- function(t) -log((1 - rowMeans(pweibull(outer(t, gamma, FUN = "/"), shape = 1/sigma))) )
zz
eval_Lambdastar <- Lambdastar(tgrid2)

lines(tgrid2, eval_Lambdastar, type = "l", col = "blue", lwd = 2)

plot(Lambdastar(tgrid), Lambdahat(tgrid, 0.1), cex = 0.2)
abline(0,1, lwd = 2, col = "blue")

plot(tgrid2, lambdastar(tgrid2))
lines(tgrid2, lambdahat(tgrid2, 0.1), col = "blue", lwd = 1.5)

# estimation of Cox model assuming BL functions are known

Xcent <- scale(X, center = TRUE, scale = FALSE)
#Xcent <- X
lambdastar_t <- lambdastar(y)#4*y^3#lambdastar(y)
Lambdastar_t <- Lambdastar(y)#y^4#Lambdastar(y)
t1_grad <- -colSums(sweep(Xcent, MARGIN = 1, FUN = "*", STATS = (1-cens)))

obj_cox_star <- function(beta, X, cens) sum( (1 - cens) * (-X%*%beta -log(lambdastar_t)) + exp(X %*% beta) * Lambdastar_t)
grad_cox_star <- function(beta, X, cens)  t1_grad + t(X) %*% (exp(X %*% beta) * Lambdastar_t) #colSums(sweep(X, MARGIN = 1, FUN = "*", STATS = exp(X %*% beta) * Lambdastar_t))

creg <- coxph(Surv(y, event = 1 - cens) ~ Xcent)
bhaz <- basehaz(creg)
#plot(lambdastar_t * exp(Xcent %*%res_star$par),
# estimate regression parameters
res_star <- optim(par = rep(0,d), fn = obj_cox_star, gr = grad_cox_star, method = "BFGS", X = Xcent, cens = cens)

plot(exp(-log(gamma)/sigma) * dweibull(y, shape = 1/sigma)/(1-pweibull(y, shape = 1/sigma)))
plot(-log(gamma)/sigma, Xcent %*% res_star$par)

plot(-log(gamma)/sigma + log(4*y^3), Xcent %*% res_star$par)#%

plot(-log(gamma)/sigma, Xcent %*% res_star$par)
plot(-log(gamma)/sigma, Xcent %*% creg$coef)
times <- bhaz$time
Delta <- c(bhaz$hazard[1], diff(bhaz))

lambdahat_0 <- lambdastar(y)#4*y^3#lambdastar(y)
Lambdahat_t <- Lambdastar(y)

plot(-log(gamma)/sigma, Xcent %*% creg$coef)
abline(0,1, col = "blue", lwd = 3)
yy <- Xcent %*% res_star$par
xx <- -log(gamma)/sigma
lm(yy ~ xx)
#plot(c(bhaz[1,1], diff(bhaz[,1])) * predict(creg, type = "risk"), lambdastar_t * exp(Xcent %*%res_star$par))

### pseudo-code

fy <- dens_y(yperm, cens)
alpha <- 0.5
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)

# compute naive estimator

creg <- coxph(Surv(yperm, event = 1 - cens) ~ Xcent)
bhaz <- basehaz(creg)
mu <- Xcent %*% creg$coef

# compute associated smoothed baseline and cumulative baseline hazard

times <- bhaz$time
Delta <- c(bhaz$hazard[1], diff(bhaz$hazard))

lambdahat_0 <- 4*yperm^3#lambdahat(yperm, 0.1)
Lambdahat_0 <- yperm^4#Lambdahat(yperm, 0.1)

fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *  exp(-exp(mu) * L0)
#obj_cox_star <- function(beta, X, cens) sum( (1 - cens) * (-X%*%beta -log(lambdastar_t)) + exp(X %*% beta) * Lambdastar_t)
nloglik_marginal <- function(mu, cens, alpha, l0, L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy))
# update weights and alpha; then repeat

iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)

# marginal approach
while(iter < maxiter){



num <- (1-alpha) * fymu(mu, cens, lambdahat_0, Lambdahat_0)
denom <- num + alpha * fy
pcur <- num/denom
alpha <- 1 - mean(pcur)

t1_grad <- -colSums(sweep(Xcent, MARGIN = 1, FUN = "*", STATS = (1-cens)*pcur))
obj_cox_star <- function(beta, X, cens) sum(pcur* ((1 - cens) * (-X%*%beta -log(lambdahat_0)) + exp(X %*% beta) * Lambdahat_0))
grad_cox_star <- function(beta, X, cens)  t1_grad + t(X) %*% (exp(X %*% beta) * Lambdahat_0 * pcur) #colSums(sweep(X, MARGIN = 1, FUN = "*", STATS = exp(X %*% beta) * Lambdastar_t))

creg <- coxph(Surv(y, event = 1 - cens) ~ Xcent)
bhaz <- basehaz(creg)
#plot(lambdastar_t * exp(Xcent %*%res_star$par),
# estimate regression parameters
res_star <- optim(par = rep(0,d), fn = obj_cox_star, gr = grad_cox_star, method = "BFGS", X = Xcent, cens = cens)
#creg <- coxph(Surv(yperm, event = 1 - cens) ~ Xcent, weights = pmax(drop(pcur),1E-6))
#bhaz <- basehaz(creg)
#mu <- Xcent %*% creg$coef
mu <- Xcent %*% res_star$par
times <- bhaz$time
Delta <- c(bhaz$hazard[1], diff(bhaz$hazard))

lambdahat_0 <- 4*yperm^3#lambdahat(yperm, 0.1)
Lambdahat_0 <- yperm^4#Lambdahat(yperm, 0.1)

#nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)

   iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

alpha
plot(-log(gamma)/sigma, Xcent %*% creg$coef)
abline(0,1, col = "blue", lwd = 3)

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

fy2 <- lambdahat(yperm, 0.1)^(1-cens) * exp(-Lambdahat(yperm, 0.1))^cens



fy <- fy2 #dens_y(yperm, cens)
alpha <- 0.5
maxiter <- 1000
tol <- 1E-4
objs <- numeric(maxiter)

# compute naive estimator

creg <- coxph(Surv(yperm, event = 1 - cens) ~ Xcent)
bhaz <- basehaz(creg)
mu <- Xcent %*% creg$coef

# compute associated smoothed baseline and cumulative baseline hazard

times <- bhaz$time
Delta <- c(bhaz$hazard[1], diff(bhaz$hazard))

lambdahat_0 <- 4*yperm^3#lambdahat(yperm, 0.1)
Lambdahat_0 <- yperm^4#Lambdahat(yperm, 0.1)

fymu <- function(mu, cens, l0, L0) (exp(mu) * l0)^(1 - cens) *  exp(-exp(mu) * L0)
nloglik_marginal <- function(mu, cens, alpha, l0, L0) sum(-log((1-alpha) * fymu(mu, cens, l0, L0) + alpha * fy))
# update weights and alpha; then repeat

iter <- 1
objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)

# marginal approach
while(iter < maxiter){



num <- (1-alpha) * fymu(mu, cens, lambdahat_0, Lambdahat_0)
denom <- num + alpha * fy
pcur <- num/denom
alpha <- 1 - mean(pcur)

t1_grad <- -colSums(sweep(Xcent, MARGIN = 1, FUN = "*", STATS = (1-cens)*pcur))
obj_cox_star <- function(beta, X, cens) sum(pcur* ((1 - cens) * (-X%*%beta -log(lambdahat_0)) + exp(X %*% beta) * Lambdahat_0))
grad_cox_star <- function(beta, X, cens)  t1_grad + t(X) %*% (exp(X %*% beta) * Lambdahat_0 * pcur) #colSums(sweep(X, MARGIN = 1, FUN = "*", STATS = exp(X %*% beta) * Lambdastar_t))

creg <- coxph(Surv(y, event = 1 - cens) ~ Xcent)
bhaz <- basehaz(creg)
#plot(lambdastar_t * exp(Xcent %*%res_star$par),
# estimate regression parameters
res_star <- optim(par = rep(0,d), fn = obj_cox_star, gr = grad_cox_star, method = "BFGS", X = Xcent, cens = cens)
#creg <- coxph(Surv(yperm, event = 1 - cens) ~ Xcent, weights = pmax(drop(pcur),1E-6))
#bhaz <- basehaz(creg)
#mu <- Xcent %*% creg$coef
mu <- Xcent %*% res_star$par
times <- bhaz$time
Delta <- c(bhaz$hazard[1], diff(bhaz$hazard))

lambdahat_0 <- 4*yperm^3#lambdahat(yperm, 0.1)
Lambdahat_0 <- yperm^4#Lambdahat(yperm, 0.1)

#nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)

   iter <- iter + 1
    objs[iter] <- nloglik_marginal(mu, cens, alpha, lambdahat_0, Lambdahat_0)
    if(objs[iter] + tol > objs[iter-1])
        break
}

alpha

plot(-log(gamma)/sigma, mu)
abline(0,1, col = "blue", lwd = 3)

### Remaining: estimate baseline hazard as well

# IDEA: use kernel estimation of the baseline hazard so long until objective increases (line search?)

