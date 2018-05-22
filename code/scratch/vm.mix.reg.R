# simulate a case of two von mises distribution

library(circular)
library(CircStats)

pis <- c(.2, .8)
betas <- list(c(1,1), c(1,1))
mus <- list(0, pi)
kapps <- c(50, 50)
n <- 20
set.seed(17)
X <- cbind(rep(1, n), rnorm(n))
mus.g <- lapply(1:length(beta), function(i) X%*%beta[[i]])

M <- length(pis)
set.seed(17)
Z <- sample(1:M, n, pis, replace = TRUE)

Y <- sapply(1:n, function(i) {
  j <- Z[i]
  y <- circular(mu[[j]] + 2*atan(c(X[i,]%*%beta[[j]])))+rvonmises(1, mu=circular(0), kappa=kapp[[j]])
  return(as.numeric(y))
})

par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(circular(Y), stack=TRUE)
plot(circular(Y)[Z==1], stack=TRUE) # need to scale by density
plot(circular(Y)[Z==2], stack=TRUE) # need to scale by density


# initialize parameters
# E-step
# compute responsibilities based on the parameter values

vm.mix.reg <- function(Y, X, M, tol = .01, maxit = 10, verbose = FALSE) 
  
  M <- 2
P <- dim(X)[2]
N <- dim(X)[1]
# initialize mu, beta, kappa
Mu <- seq(0, 2*pi, length.out=M)

# beta is a 3-dimensional array
# beta[j,,] is a 1xd dimensional vector of regression coefficients
Beta <- array(0, dim = c(M,P,1))
G <- array(0, dim = c(M,N,1))

Beta[1,,] <- c(1,1)
Beta[2,,] <- c(1,1)

for (j in 1:M) {
  G[j,,] <- 2*atan(X %*% Beta[j,,])
}

# kappa is a 2-dimensional array
# kappa[j,] is the 1/variance of the j-th class
Kappa <- array(0, dim=c(M, P))
Kappa[1,] <- 5
Kappa[2,] <- 5

Lambda <- array(0, dim=c(M))
Lambda <- c(.5,.5)

# compute responsibilites
library(CircStats)

W.j <- t(sapply(1:M, function(j) {
  Lambda[j]*dvm(theta = Y[j], mu = Mu[j] + G[j,,], kappa = Kappa[j,]) }))
W <- t(t(W.j)/colSums(W.j))


# M-step
# given the responsibilities,
# compute estimates 

y <- Y
x <- X
w <- W[j,]
init <- Beta[j,,]
lambda <- Lambda[j]

# 
M.step <- lapply(1:M, function(j) {
  get.LmCircularclRad(y=Y, x=X, w=W[j,], init=c(0,0))$log.lik
})
M.step.log.lik.prev <- sapply(M.step,"[[","log.lik")

log.lik <- sapply(1:M, function(j) {
  res <- M.step.[[j]]
  (-1)*log(besselI(x = res$kappa, nu = 0, expon.scaled = FALSE)) + res$kappa*(cos(y-res$mu-2*atan(x%*%res$beta))))
return(res)
})

W <- t(t(log.lik)/colSums(log.lik))

M.step <- lapply(1:M, function(j) {
  get.LmCircularclRad(y=Y, x=X, w=W[j,], init=c(0,0))$log.lik
})
M.step.log.lik.new <- sapply(M.step,"[[","log.lik")

diff <- M.step.log.lik.new - M.step.log.lik.prev
while(diff > tol) {
  M.step <- lapply(1:M, function(j) {
    get.LmCircularclRad(y=Y, x=X, w=W[j,], init=c(0,0))$log.lik
  })
  M.step.log.lik.new <- sapply(M.step,"[[","log.lik")
  M.step.log.lik.prev <- M.step.log.lik.prev  
}

library(CircStats)
#lm.circular(circular(y),x,init=c(0,0), type = "c-l")

get.LmCircularclRad <- function(y, x, w, init, verbose=TRUE, tol=1e-1) {
  n <- length(y)
  y <- y%%(2*pi)
  betaPrev <- init  
  S <- sum(w*sin(y-2*atan(x%*%betaPrev)))/n
  C <- sum(w*cos(y-2*atan(x%*%betaPrev)))/n
  R <- sqrt(S^2 + C^2)
  mu <- atan2(S,C)
  k  <- A1inv(R)
  #  if (k < 0) {k <- 1}
  diff <- tol+1
  iter <- 0
  while (diff > tol){
    iter <- iter + 1
    u <- k*sin(y - mu - 2*atan(x%*%betaPrev))
    A <- diag(k*A1(k), nrow=n)
    g.p <- diag(apply(x, 1, function(row, betaPrev) 2/(1+(t(betaPrev)%*%row)^2), 
                      betaPrev=betaPrev), nrow=n)
    g.p.w <- g.p %*% diag(sqrt(w))
    D <- g.p.w%*%x
    betaNew <- lm(t(D)%*%(u+A%*%D%*%betaPrev) ~ t(D)%*%A%*%D - 1)$coefficients
    diff <- max(abs(betaNew - betaPrev))
    
    cat(diff, ";", betaNew, ";",betaPrev, "\n")
    
    betaPrev <- betaNew
    
    S <- sum(w*sin(y-2*atan(x%*%betaPrev)))/n
    C <- sum(w*cos(y-2*atan(x%*%betaPrev)))/n
    R <- sqrt(S^2 + C^2)
    mu <- atan2(S,C)
    k  <- A1inv(R)
    
    if (verbose){
      log.lik <- sum(w*(-1)*log(besselI(x = k, nu = 0, expon.scaled = FALSE)) + 
                       w*k*(cos(y-mu-2*atan(x%*%betaNew))) )
      cat("Iteration ", iter, ":    Log-Likelihood = ", log.lik, "\n")
    }
  }
  
  log.lik <- sum(w*(-1)*log(besselI(x = k, nu = 0, expon.scaled = FALSE)) + 
                   w*k*(cos(y-mu-2*atan(x%*%betaNew))))
  
  cov.beta <- solve(t(D)%*%A%*%D)
  se.beta <- sqrt(diag(cov.beta))
  se.kappa <- sqrt(1/(n*(1-A1(k)^2-A1(k)/k)))
  se.mu <- 1/sqrt((n-ncol(x))*k*A1(k))
  t.values <- abs(betaNew/se.beta)
  p.values <- 1-pnorm(t.values)
  betaNew <- as.vector(betaNew)
  output <- list()
  output$x <- x
  output$y <- y
  output$beta <- betaNew
  output$cov.beta <- cof.beta
  output$se.beta <- se.beta
  output$mu <- mu
  output$se.mu <- se.mu
  output$kappa <- k
  output$se.kappa <- se.kappa
  output$log.lik <- log.lik 
  return(output)
}
