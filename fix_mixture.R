###--------------------------------------------------------------------------------------------
### X:              --- n x p design matrix
### y:              --- n x 1 observed response
### alph:           --- the level of mismatch parameter, default is 0.5
### g:              --- density function of the marginal distribution of y
### control:        --- additional parameters
### init:           --- initial estimator of Expectation-Maximization (EM) algorithm, default is "LS" (naive estimator) ignoring that (xi,yi) - pairs are subject to mismatches
### tol:            --- tolerance level, default is 1E-4
### maxit:          --- maximum iteration of Expectation-Maximization (EM) algorithm,  default is 100
### family:         --- model to be fit for mixture model, default is "Gaussian"
### weights:        --- the observation weights, default is 1 for each observation

#y <- yperm ; alph <- alpha; family <- "binomial"; init <- "IRLS"; tol= 1E-4; maxit= 100; family = "poisson"; family = Gamma(link = "log")
fit_mixture <- function(X, y, alph, init, family, tol= 1E-4, maxit= 100, shape = 5){
  ########################################
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  if(!is.null(init)){
    init <- init
  } else
    init <- "IRLS"
  
  if(!is.null(tol)){
    tol <- tol
  } else
    tol <- 1E-4
  
  if(!is.null(maxit))
  {
    maxit <- maxit
  } else
    maxit <- 100
  
  if(!is.null(family)){
    family <- family
  }else
    family <- "Gaussian"
  
  y <- as.numeric(y)
  
  if(length(y) != n)
    stop("Error: Length of 'y' must be equal to the number of rows of 'X' \n")
  
  if(!(init %in% c("IRLS", "robust", "LL")))
    stop("Error: 'init' must be either 'IRLS', 'robust' or 'LL' \n")
  
  if(init == "IRLS"){
    if(identical(family, "Gamma(link = log)")){
      glm0 <- glm(yperm ~ X - 1, family = Gamma(link = "log"))
      betahat <- coef(glm0)
    }
    else{
      glm0 <- glm(yperm ~ X - 1, family = family)
      betahat <- coef(glm0)
    }
    alphahat <- alph    
  }
  
  if(init == "robust"){
    require(MASS)
    rlmsol <- rlm(x=X, y=y)
    betahat <- coef(rlmsol)
    betahat_LS <- solve(crossprod(X), crossprod(X, y))
    alphahat <- max(min(1, sqrt(sum(betahat_LS^2))/sqrt(sum(betahat^2))),0.01)   
    
    require(robustbase)
    glm_robust <- glmrob(yperm ~ X - 1, family, method="Mqle", control = glmrobMqle.control(tcc=3.5))
    betahat_robust <- coef(glm_robust)
    alphahat <- max(min(1, sqrt(sum(bbetahat_robust ^2))/sqrt(sum(betahat^2))),0.01)   
  }
  
  if(init == "LL"){
    Q = (1 - alph - alph/(n-1)) * diag(n) + (alph/(n-1)) * matrix(1, n, n)
    betahat <- lahiri_larsen(X, y, Q, family)
    alphahat <- alph    
  }
  
  ########################## Main ###################################
  betacur <- betahat
  
  if(identical(family ,"binomial")){
    p <- mean(yperm)
    fymu <- function(mu) (mu^yperm) * ((1 - mu)^(1 - yperm))   
    fy <- (p^yperm) * (1-p)^(1-yperm)
    mu_fun <- function(beta) plogis(X%*%beta)
    mucur <- mu_fun(betacur)
  }
  
  if(identical(family,"poisson")){
    taby <- table(yperm)/n
    fymu <- function(mu) dpois(yperm, mu)
    fy <- as.numeric(taby[match(yperm, as.numeric(names(taby)))])
    mu_fun <- function(beta) exp(X%*%beta)
    mucur <- mu_fun(betacur)
  }
  
  if(identical(family, "Gamma(link = log)")){
    dens <- density(yperm)
    ############################################
    fymu <- function(mu) dgamma(yperm, shape = shape, scale = mu/shape)
    fy <- approx(dens$x,dens$y,xout=yperm)$y
    mu_fun <- function(beta) exp(X%*%beta)
    mucur <- mu_fun(betacur)
    family = Gamma(link = "log")
  }
  
  nloglik <- function(mu, alpha) sum(-log((1-alpha) * fymu(mu) + alpha * fy))
  it <- 1
  objs <- numeric(maxit)           
  nloglik_cur <- nloglik(mucur, alphahat)
  objs[it] <- nloglik_cur
  
  ### iterations
  while(it < maxit){
    num <- (1-alphahat) * fymu(mu)
    denom <- num + alphahat * fy
    pcur <- num/denom
    alphahat <- 1 - mean(pcur)
    wglmfit <- glm(yperm ~ X - 1, family = family, weights = pcur)
    betacur <- coef(wglmfit)
    mu <- mu_fun(betacur)
    it <- it + 1
    objs[it] <- nloglik(mu, alphahat)
    if(objs[it] + tol > objs[it-1])
      break
  }
  ########################## Output ###################################
  res_hat = y - X%*%betahat
  fitted_values = X%*%betahat
  linear_predictors = X%*%betahat
  coverged = (objs[it-1] - objs[it]) < tol
  
  return(list(beta = betahat, alpha = alphahat,
              residuals = res_hat, fitted.values = fitted_values, linear.predictors = linear_predictors, 
              family = family, coverge = coverged, iteration = it, objective = objs[1:it], shape = shape))
}

### X:       --- n x p design matrix
### y:       --- n x 1 observed response
### f0:      --- density function of the marginal distribution of y
### control: --- additional parameters

fit_mixture_guassian <- function(X, y, f0, control = list(init = c("LS", "robust"), tol = 1E-4, maxit = 100, gamma = 0.95)){
  
  if(!is.null(control$init))
    init <- control$init[1]
  else
    init <- "LS"
  if(!is.null(control$tol))
    tol <- control$tol
  else
    tol <- 1E-4
  if(!is.null(control$maxit))
    maxit <- control$maxit
  else
    maxit <- 100
  if(!is.null(control$gamma))
    gamma <- control$gamma
  else
    gamma <- 0.95
  
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  y <- as.numeric(y)
  
  if(length(y) != n)
    stop("Error: Length of 'y' must be equal to the number of rows of 'X' \n")
  
  if(!is.function(f0))
    stop("Error: 'f0' must be a function (density) \n")
  
  
  if(!(init %in% c("LS", "robust")))
    stop("Error: 'init' must be either 'LS' or 'robust' \n")
  
  if(init == "LS"){
    betahat <- solve(crossprod(X), crossprod(X, y))
    sigmahat <- sqrt(sum((y - X %*% betahat)^2)/(n-p))
    ##  sigmahat <- sqrt(sum((y - X %*% betahat)^2)/n)
    alphahat <- 0.5    
  }
  
  if(init == "robust"){
    
    require(MASS)
    rlmsol <- rlm(x=X, y=y)
    betahat <- coef(rlmsol)
    sigmahat <- rlmsol$s
    betahat_LS <- solve(crossprod(X), crossprod(X, y))
    alphahat <- max(min(1, sqrt(sum(betahat_LS^2))/sqrt(sum(betahat^2))),0.01)     
  }
  
  
  nloglik <- function(beta, sigma, alpha){
    mean(-log(alpha * (1/sqrt(2 * pi * sigma^2)) * exp(-(y - (X %*% beta))^2/(2 * sigma^2)) 
              + (1-alpha) * f0(y)))
  }
  
  
  it <- 0
  conv <- numeric(maxit)           
  nloglik_cur <- nloglik(betahat, sigmahat, alphahat)
  conv[it + 1] <- nloglik_cur
  
  ### iterations
  
  while(it < maxit){
    r <- X %*% betahat - y
    
    pnum <-  alphahat * (1/sqrt(2 * pi * sigmahat^2)) * exp(-r^2/(2*sigmahat^2))
    pdenom <- pnum + (1 - alphahat) * f0(y)
    pcur <- pnum / pdenom
    w <- pcur
    
    Xw <- X
    Xw <- sweep(X, 1, STATS = sqrt(w), FUN="*") 
    yw <- sqrt(w) * y
    
    betahatnew <- solve(crossprod(Xw), crossprod(Xw, yw))        
    sigmahatnew <- sqrt(sum((yw - Xw %*% betahatnew)^2) / sum(w))
    
    alphahatnew = mean(pcur);
    nloglik_new = nloglik(betahatnew, sigmahatnew, alphahatnew)
    
    if((nloglik_cur - nloglik_new) < -tol^2)
      break
    else{
      betahat = betahatnew;
      alphahat = alphahatnew;
      sigmahat = sigmahatnew;
      it <- it+1
      conv[it+1] = nloglik_new
      if((nloglik_cur - nloglik_new) < tol)
        break
      else{
        nloglik_cur <- nloglik_new
      }
    }
    
  }
  
  
  return(list(betahat = betahat, sigmahat = sigmahat, alphahat = alphahat, conv = conv[1:(it+1)]))
  
}