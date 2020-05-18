#1.
#a).

avg <- c()
sdv <- c()
for (n in c(10,100,1000,10000)) {
  dist <- rbeta(n = n,shape1 = 2+5,shape2 = 2+15)
  plot(dist)
  avg <- c(avg,mean(dist))
  sdv <- c(sdv, sd(dist))
}

#E = (alpha + s)/(alpha + s +beta+ f)
# = (alpha + 5)/(alpha + 5+beta + 15)
# = (alpha + 5)/(alpha + beta + 20)
# = 7/24 = 0.2916667

plot(c(10,100,1000,10000),avg, type = "l", ylim = c(0,0.4), ylab = "sd, mean")
lines(c(10,100,1000,10000),sdv, col="red")

#b).
dist <-rbeta(n = 10000,shape1 = 2+5,shape2 = 2+15)
sum(dist>0.3)/length(dist)

1-pbeta(0.3,shape1 = 2+5,shape2 = 2+15)

#c).
log_odds<-log(dist/(1-dist))
hist(log_odds)
density(log_odds)


#2.
#a).
obs <- c(44,25,45,52,30,63,19,50,34,67)

n <- 10000
tau2 <- sum((log(obs)-3.7)^2)/length(obs)
posterior <- (length(obs)*tau2)/rchisq(n = n, df = length(obs))
hist(posterior)

#We compare mean of simulated posterior and theoretical value
#simulated:
mean(posterior)
#theoretical:
(length(obs)*tau2)/(length(obs)-2)

#b).
G <- 2* pnorm(q = sqrt(posterior)/sqrt(2), mean = 0, sd = 1)-1
hist(G, breaks = 30)

#c).
#90% equal tail credible interval
G_sort <- sort(G)
G_cut <- G_sort[(length(G_sort)*0.05):(length(G_sort)*0.95-1)]
#hist(G_cut)
hist(G, breaks = 30)
abline(v=min(G_cut), col = "red")
abline(v=max(G_cut), col = "red")

#90% Highest Posterior Density interval
density(G)

#3.
#a).
wind_dir <- c(40,303,326,285,296,314,20,308,299,296)
radians <- c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)

mu <- 2.39
k <- seq(0.001, 10.001, by = 0.001)
n <- length(radians)
ml <- exp(k*sum(cos(radians-mu)))/(besselI(k, 0))^n
prior <- exp(-k)
#posterior <- exp(k*sum(cos(radians-mu))-k)/(besselI(k, 0))^n
posterior <- ml * prior #distribution the same but values very small
plot(k,posterior, type = "l")

#b).
abline(v=k[which(posterior==max(posterior))], col = "red")
