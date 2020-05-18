#1.
#a).
n <- 1
data <- read.csv("/home/erik/Documents/SML/Semester 2/Bayesian Learning/Labs/Lab2/TempLinkoping.txt", sep="")
X <- cbind(1,data$time,data$time^2)

mu_0 <- c(-10,130,-130)
sigma2_0 <- 1
omega_0 <- 0.1*diag(3)
vu_0 <- 4

library(MASS)
library(mvtnorm)

for (i in 1:100) {
  sigma2 <- (sigma2_0*vu_0)/rchisq(n = 1, df = vu_0)
  #beta_sigma2 <- mvrnorm(n = n, mu = mu_0, Sigma = sigma2*solve(omega_0))
  beta_sigma2 <- rmvnorm(n = 1, mean = mu_0, sigma = sigma2*solve(omega_0))
  
  #temp <- X%*%beta_sigma2
  E <- rnorm(n = 1, mean = 0,sd = sqrt(sigma2))
  temp <- X%*%t(beta_sigma2) + E
  
  if (i == 1) {
    plot(temp, type = "l", ylim = c(-50,50))
  }
  else{
    lines(temp, type = "l")
  }
}

#b).
X <- cbind(1,data$time,data$time^2)
y <- data$temp
n <- length(data$time)

beta_hat <- solve(t(X) %*% X) %*% (t(X) %*% y)
mu_n <- solve(t(X)%*%X+omega_0)%*%(t(X)%*%X%*%beta_hat+omega_0%*%mu_0)
omega_n <- t(X)%*%X+omega_0
vu_n <- vu_0 + n #what n?
vu_nsigma2 <- vu_0*sigma2_0 + (t(y)%*%y+t(mu_0)%*%omega_0%*%mu_0-t(mu_n)%*%omega_n%*%mu_n)

sigma2_n <- var(data$temp)
sigma2_i <- c()
beta_0 <- c()
beta_1 <- c()
beta_2 <- c()
f_time <- c()
f_median <- c()

for (i in 1:100) {
  sigma2 <- (vu_nsigma2)/rchisq(n = 1, df = vu_n)
  beta_sigma2 <- rmvnorm(n = 1, mean = mu_n, sigma = sigma2[1,1]*solve(omega_n))
  
  #temp <- X%*%beta_sigma2
  E <- rnorm(n = 1, mean = 0,sd = sqrt(sigma2))
  temp <- X%*%t(beta_sigma2) + E
  f_median[i] <- median(temp)
  f_time <- cbind(f_time, as.vector(temp))
  
  if (i == 1) {
    plot(temp, type = "l", ylim = c(-50,50))
  }
  else{
    lines(temp, type = "l")
  }
  beta_0[i] <- beta_sigma2[1,1]
  beta_1[i] <- beta_sigma2[1,2]
  beta_2[i] <- beta_sigma2[1,3]
  sigma2_i[i] <- sigma2 
}

hist(beta_0, main = "marginal posterior beta_0")
hist(beta_1, main = "marginal posterior beta_1")
hist(beta_2, main = "marginal posterior beta_2")
hist(sigma2_i, main = "marginal posterior sigma²")



#95% equal tail credible interval
G <- f_time[,1]
G_sort <- sort(G)
#G_cut <- G_sort[(length(G_sort)*0.025):(length(G_sort)*0.975-1)]
plot(G)
abline(v=(length(G_sort)*0.025), col = "red")
abline(v=(length(G_sort)*0.975-1), col = "red")

#c).
x_max <- -beta_1/(2*beta_2)
plot(density(x_max))
abline(v=X[which.max(f_median)], col = "red")


#d).


#2.
#a).
data <- read.table("/home/erik/Documents/SML/Semester 2/Bayesian Learning/Labs/Lab2/WomenWork.dat", header = TRUE)
n <- dim(data)[1]

x <- as.matrix(data[,2:9])
y <- as.vector(data[,1])

nPara <- dim(x)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara))
tau <- 10
sigma <- tau^2*diag(nPara)


LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
  
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  # evaluating the log-likelihood                                    
  logLik <- sum( linPred*y -log(1 + exp(linPred)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  
  # evaluating the prior
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLik + logPrior)
}


#initVal <- as.vector(rmvnorm(n = 1, mean = rep(0,dim(x)[2]), sigma = tau^2*diag(dim(x)[2])))
initVal <- as.vector(rep(0,dim(x)[2]))
posterior <- optim(initVal,LogPostLogistic,gr=NULL,y,x,mu,sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
posterior_mode <- posterior$par
hessian_matrix <- posterior$hessian
hessian_sigma <- -solve(posterior$hessian)

print("Optimal beta:")
print(posterior_mode)

#
print(hessian_sigma)

#95% equal tail credible interval
for (i in 1:1000) {
  G[i] <- rnorm(1, posterior_mode[7], hessian_sigma[7,7])
}

G_sort <- sort(G)
#G_cut <- G_sort[(length(G_sort)*0.05):(length(G_sort)*0.95-1)]
plot(G_sort)
abline(v=(length(G_sort)*0.05), col = "red")
abline(v=length(G_sort)*0.95-1, col = "red")
#Yes it is an important feature, it's absolute beta value is the highest of the entire model.
#This means the value of the feature has a strong impact on the result.

glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)
glmModel$coefficients

#Looking at the estimated parameters using maximum likelihood we can see that tey are very close the the ones we predicted.

#b).
X <- c(1,10,8,10,1,40,1,1)
result <- c()
for (i in 1:100) {
  beta_posterior <- rmvnorm(1, posterior_mode, hessian_sigma)
  y <- exp(t(X) %*% t(beta_posterior))/(1+(exp(t(X) %*% t(beta_posterior))))
  result[i] <- rbinom(1,1, y)
}
plot(density(result), main = "p. distribution for Work")

sum(result==1)

#c).
X <- c(1,10,8,10,1,40,1,1)
result <- c()
for (i in 1:100) {
  beta_posterior <- rmvnorm(1, posterior_mode, hessian_sigma)
  y <- exp(t(X) %*% t(beta_posterior))/(1+(exp(t(X) %*% t(beta_posterior))))
  result[i] <- rbinom(1,10, y)
}
plot(density(result), main = "working women" , xlab = "working women")
