# adapt from mcsm's Chapter.7.R (https://github.com/cran/mcsm/blob/master/demo/Chapter.7.R)

xdata = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
Time = c(94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.05, 1.05, 2.10, 10.48)
nx = length(xdata)
nsim = 10^4; alpha = 1.8; gamma= 0.01; delta=1
lambda = array(xdata*Time/sum(Time), dim=c(nsim,nx))
beta = rep(gamma*delta,nsim)
for (i in 2:nsim) {
  for (j in 1:nx){
    lambda[i,j] = rgamma(1,shape=xdata[j]+alpha, rate=Time[j]+beta[i-1])
  }
  beta[i] = rgamma(1,shape=gamma+nx*alpha,rate=delta+sum(lambda[i,]))
}
par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1.5)); st=.5*nsim
hist(lambda[,1][st:nsim],breaks=25,col="grey",xlab="",main=expression(lambda[1]))
hist(lambda[,2][st:nsim],breaks=25,col="gold3",xlab="",main=expression(lambda[2]))
hist(beta[st:nsim],breaks=25,col="sienna",xlab="",main=expression(beta))
acf(lambda[,1][st:nsim],lwd=3,xlab="")
acf(lambda[,2][st:nsim],lwd=3,xlab="",col="gold3")
acf(beta[st:nsim],lwd=3,col="sienna",xlab="")