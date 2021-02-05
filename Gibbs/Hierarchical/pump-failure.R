# adapt from mcsm's Chapter.7.R (https://github.com/cran/mcsm/blob/master/demo/Chapter.7.R)
rm(list = ls())
xdata = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
time = c(94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.05, 1.05, 2.10, 10.48)
nx = length(xdata)
nsim = 10^4; alpha = 1.8; gamma= 0.01; delta=1
set.seed(1234)
lambda = jitter(xdata / time)
beta = rgamma(1, shape = gamma + 10*alpha, rate = delta + sum(lambda))
for (i in 2:nsim) {
  lambda = rbind(lambda, rgamma(nx, shape = xdata + alpha, rate = time + beta[i-1]))
  beta = c(beta, rgamma(1, shape = gamma+nx*alpha, rate = delta+sum(lambda[i,])))
}
par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1.5)); st=.5*nsim
hist(lambda[,1][st:nsim],breaks=25,col="grey",xlab="",main=expression(lambda[1]))
hist(lambda[,2][st:nsim],breaks=25,col="gold3",xlab="",main=expression(lambda[2]))
hist(beta[st:nsim],breaks=25,col="sienna",xlab="",main=expression(beta))
acf(lambda[,1][st:nsim],lwd=3,xlab="")
acf(lambda[,2][st:nsim],lwd=3,xlab="",col="gold3")
acf(beta[st:nsim],lwd=3,col="sienna",xlab="")

## single chain KS
ks = NULL
G = 10
for (t in seq(nsim/10, nsim, length.out = 100)) {
  beta1 = beta[1:(t/2)]
  beta2 = beta[(t/2):(1:t/2)]
  beta1 = beta1[seq(1, t/2, by = G)]
  beta2 = beta2[seq(1, t/2, by = G)]
  ks = c(ks, ks.test(beta1, beta2)$p)
}

## dual chain KS
oldbeta = beta[seq(1, nsim, by = G)]
olks = ks

lambda = jitter(xdata / time)
beta = rgamma(1, gamma + 10*alpha) / (delta + sum(lambda))
for (t in 2:nsim) {
  lambda = rbind(lambda, rgamma(10, xdata + alpha) / (time + beta[t-1]))
  beta = c(beta, rgamma(1, gamma + 10*alpha) / (delta + sum(lambda[t,])))
}

beta = beta[seq(1, nsim, by = G)]
ks = NULL
for (t in seq(nsim/(10*G), nsim/G, length.out = 100))
  ks = c(ks, ks.test(beta[1:t], oldbeta[1:t])$p)

par(mfrow = c(2, 1))
plot(seq(1, nsim, le = 100), olks, xlab = "iter", ylab = "p-value")
plot(seq(1, nsim, le = 100), ks, xlab = "iter", ylab = "p-value")