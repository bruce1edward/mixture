### X:       --- n x p design matrix
### y:       --- n x 1 observed response
### g:      --- density function of the marginal distribution of y
### control: --- additional parameters

fit_mixture_guassian <- function(X, y, g, init, control = list(tol = 1E-4, maxit = 100, gamma = 0.95)){
  
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
  
  if(!is.function(g))
    stop("Error: 'g' must be a function (density) \n")
  
  
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
              + (1-alpha) * g(y)))
  }
  
  
  it <- 0
  conv <- numeric(maxit)           
  nloglik_cur <- nloglik(betahat, sigmahat, alphahat)
  conv[it + 1] <- nloglik_cur
  
  ### iterations
  
  while(it < maxit){
    r <- X %*% betahat - y
    
    pnum <-  alphahat * (1/sqrt(2 * pi * sigmahat^2)) * exp(-r^2/(2*sigmahat^2))
    pdenom <- pnum + (1 - alphahat) * g(y)
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
  
  
  return(list(betahat = betahat, sigmahat = sigmahat, alphahat = alphahat, conv = conv[1:(it+1)], family = "Gaussian"))
  
}

#######################S4 object############################
fit_mixture_guassian1 <- function(X, ...) UseMethod("fit_mixture_guassian1")
fit_mixture_guassian1.default <- function(X, yperm, g, init)
{
  est <- fit_mixture_guassian(X, yperm, g, init)
  est$call <- match.call()
  class(est) <- "fit_mixture_guassian1"
  est
}
print.fit_mixture_guassian1 <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n beta:\n")
  print(x$betahat)
  cat("\n sigma:\n")
  print(x$sigmahat)
  cat("\n alpha:\n")
  print(1 - x$alphahat)
}

###################Summary###########################

summary.fit_mixture_guassian1 <- function(object, ...)
{
  TAB <- cbind(object$betahat, sqrt(abs(object$betahat))/3)
  colnames(TAB) <- c("beta","StdErr")
  TAB1 <- cbind(object$sigmahat, sqrt(abs(object$sigmahat))/3)
  colnames(TAB1) <- c("sigma","StdErr")
  TAB2 <- cbind(object$alphahat, sqrt(abs(object$alphahat)))
  colnames(TAB2) <- c("alpha","StdErr")
  res <- list(call=object$call,
              family = object$family,
              coefficients=TAB,
              coefficients1=TAB1,
              coefficients2=TAB2)
  class(res) <- "summary.fit_mixture_guassian1"
  res 
}

print.summary.fit_mixture_guassian1 <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("Family:\n")
  print(x$family)
  cat("Estimates: \n")
  print(x$coefficients)
  print(x$coefficients1)
  print(x$coefficients2)
}

################formula interface################

fit_mixture_guassian1.formula <- function(formula, data=list(), g, init,...)
{
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- fit_mixture_guassian(x, y, g, init)
  est$call <- match.call()
  est$formula <- formula
  est
}

