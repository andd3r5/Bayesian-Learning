#2.
#a).
data <- read.table("~/Documents/Bayesian Learning/Lab 3/eBayNumberOfBidderData.dat", header = TRUE)


# beta <- glm(formula = nBids ~ . - Const, data = data)
# beta_hat <- as.matrix(beta$coefficients)
# 
# 
# ml <- 0
# for (i in 1:nrow(data)) {
#   lambda <- exp(as.matrix(x[i,]) %*% beta_hat)
#   ml <- ml + data$nBids[i] * log(lambda)-lambda
# }

poisson_beta <- glm(formula = nBids ~ . - Const, data = data, family = 'poisson')
summary(poisson_beta)

#looking at the pvalue of the z scores it appears that VerifyID, Sealed, MajBlem, LogBook and MinBidShare are significant since they are all < 0.05.

#b).
library(mvtnorm)
X <- as.matrix(data[,2:ncol(data)])
y <- as.matrix(data[,1])
nPara <- dim(X)[2]


# Setting up the prior
mu <- as.vector(rep(0,nPara))
# tau <- 10
Sigma <- 100*solve(t(X)%*%X)
betaVect <- rmvnorm(n = 1,mean = mu,sigma = 100*solve(t(X)%*%X))


PostPoisson <- function(betaVect,y,X,Sigma){
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  # evaluating the log-likelihood
  #logLik <- sum( linPred*y -log(1 + exp(linPred)));
  logLik <- sum(y * linPred - exp(linPred)) 

  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
      # evaluating the prior
      logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log = TRUE);
      # add the log prior and log-likelihood together to get log posterior
      return(logLik + logPrior)
}

#initVal <- as.vector(rmvnorm(n = 1, mean = rep(0,dim(x)[2]), sigma = tau^2*diag(dim(x)[2])))
initVal <- as.vector(rep(0,dim(X)[2]))
posterior <- optim(initVal,PostPoisson,gr=NULL,y,X,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
posterior_mode <- posterior$par
hessian_matrix <- posterior$hessian
hessian_sigma <- -solve(posterior$hessian)
print("Optimal beta:")
print(posterior_mode)
print(hessian_sigma)

#c).


RWMSampler <- function(logPostFunc, theta_0, hessian_sigma, c, niter, ...){
  
  theta <- c()
  i <- 1
  reject <- 0
  while (i <= niter) {

    if (i == 1) {
      theta_p <- as.vector(rmvnorm(n = 1, mean = theta_0, sigma = c * hessian_sigma))
      alpha <- exp(logPostFunc(theta_p, ...)-logPostFunc(mu, ...))
    } else{
      theta_p <- as.vector(rmvnorm(n = 1, mean = as.vector(theta[i-1,]), sigma = c * hessian_sigma)) #maybe t(theta)
      alpha <- exp(logPostFunc(theta_p, ...)-logPostFunc(as.vector(theta[i-1,]), ...))
    }
    
    if (runif(n = 1, min = 0,max = 1) <= alpha) {
      theta <- rbind(theta, theta_p)
      i <- i + 1
    }else{reject <- reject + 1}
  }
  # theta_p <- as.vector(rmvnorm(n = 1, mean = theta_0, sigma = c * hessian_sigma))
  # logPostFunc(theta_p, ...)
  return(list("theta"=theta, "reject"=reject, "accept"=niter,"acceptance_rate"=(niter)/(reject+niter)))
}

X <- as.matrix(data[,2:ncol(data)])
y <- as.matrix(data[,1])
nPara <- dim(X)[2]
# Setting up the prior
mu <- as.vector(rep(0,nPara))
Sigma <- 100*solve(t(X)%*%X)
betaVect <- rmvnorm(n = 1,mean = mu,sigma = 100*solve(t(X)%*%X))

result_poisson <- RWMSampler(PostPoisson, mu, hessian_sigma, 0.65, 1000, y, X, Sigma)
result_poisson$theta[nrow(result_poisson$theta),]

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(as.vector(result_poisson$theta[,1]), type = "l", ylim = c(min(result_poisson$theta),max(result_poisson$theta)), col = 1)
for (i in 2:ncol(result_poisson$theta)) {
  lines(result_poisson$theta[,i], col = i)
}
legend("topright", inset=c(-0.22,0), legend=c("Intercept", sprintf("beta[%s]",seq(1:(ncol(result_poisson$theta)-1)))),
       col=1:ncol(result_poisson$theta), lty = 1)
par(xpd=FALSE)

#d).
x_d <- c(1,1,1,1,0,0,0,1,0.5)
x_d %*% result_poisson$theta[nrow(result_poisson$theta),]
