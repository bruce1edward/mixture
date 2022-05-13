set.seed(1234)
require(tictoc)
#######################################
n <- 100000
d <- 3
k_n <- 0.1
b <- 1
k <- round(n * k_n)
X <- matrix(rnorm(n*d), nrow = n, ncol = d)
beta <- rnorm(d)
beta <- beta/ sqrt(sum(beta^2))
beta <- beta*b

###Linear Regression###################
sigma <- 1
y <- rnorm(n, X %*% beta, sigma)
yperm <- y
yperm[1:k] <- y[sample(1:k)]
alpha <-0.5
tausq <- mean(yperm^2)
f0 <- function(z) dnorm(z, mean = 0, sd = sqrt(tausq))
tic()
res <- fit_mixture_guassian(X, y, f0, control = list(init = "robust"))
toc()

###Logistic Regression#################

mu <- plogis(X %*% beta)
y <- rbinom(n, size = 1, prob = mu)
yperm <- y
yperm[1:k] <- y[sample(1:k)]
alpha <-0.5

tic()
mod1 <- fit_mixture(X, yperm, alpha, init = "IRLS", family = "binomial")
toc()

###Poisson Regression##################
mu <- exp(X %*% beta)
y <- rpois(n, lambda = mu)
yperm <- y
yperm[1:k] <- y[sample(1:k)]
alpha <-0.5

tic()
mod2 <- fit_mixture(X, yperm, alpha, init = "IRLS", family = "poisson")
toc()

###Gamma Regression (log link)##########
shape <- 5
mu <- exp(X%*%beta)  # use log-link for convenience
y <- rgamma(n, shape = shape, rate = shape/mu)
dens <- density(y)
yperm <- y
yperm[1:k] <- y[sample(1:k)]
g0 <- function(z) dnorm(z, mean = 0, sd = sqrt(tausq))
alpha <-0.5

tic()
mod3 <- fit_mixture(X, yperm, alpha, init = "IRLS", family = "Gamma(link = log)", shape = 5)
toc()