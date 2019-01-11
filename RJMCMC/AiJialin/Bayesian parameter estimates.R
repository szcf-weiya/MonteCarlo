# Copyright (C) 2012  Ai Jialin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# See http://www1.maths.leeds.ac.uk/~voss/projects/2011-RJMCMC/ for
# information about this R code.

# value of mu
mu<-runif(2,min=-10,max=10)
mu

# value from x0 to xm
m=100
x_axis <- rnorm(m,mean=mu[1],sd=1)
y_axis <- rnorm(m,mean=mu[2],sd=1)
x <- matrix (c(x_axis,y_axis),nrow=m,ncol=2,byrow=FALSE)
mean(x[,1])
mean(x[,2])

# value of log(p(x|mu))
log_p_condition <- function(x,mu){
    sum=0
    for (i in 1:m){
        sum=sum+log(1/(2*pi))-0.5*t(x[i,]-mu)%*%(x[i,]-mu)
    }
    return(sum)
}

#value of log(p(mu))
log_p_prior <- function(mu){
    if(mu[1]>=-10 && mu[1]<=10 && mu[2]>=-10 && mu[2]<=10){return (log(0.0025))}
    else return (-Inf)
}


# the log value of non.normalised target density
log_pi.Z <- function(mu){
    return(log_p_condition(x,mu)+log_p_prior(mu))
}


# the acceptance probability
alpha <- function(mu,mu_tilde){
  return(exp(log_pi.Z(mu_tilde)-log_pi.Z(mu)))
}


#generate paths from the MH process

MH <- function(n,mu,sigma){
    path <- matrix(0,nrow=n,ncol=2)
    for (i in 1:n){
        mu_tilde <- rnorm (2,mu,sigma)
        U <- runif (1)
        if (U < alpha(mu,mu_tilde)){
           mu <- mu_tilde
        }
        path[i,] <- mu
    }
    return(path)
}


# M-H Markov chain with burn-in period
par(mfrow=c(3,2))
t <- c()
for(i in 1:3000){t<-c(t,i)}

y <- MH(3000,c(0,0),1)
plot(t,y[t,1],xlab="estimating time",ylab="x axis of mu",type="l")
plot(t,y[t,2],xlab="estimating time",ylab="y axis of mu",type="l")

y <- MH(3000,c(0,0),0.1)
plot(t,y[t,1],xlab="estimating time",ylab="x axis of mu",type="l")
plot(t,y[t,2],xlab="estimating time",ylab="y axis of mu",type="l")

y <- MH(3000,c(0,0),0.01)
plot(t,y[t,1],xlab="estimating time",ylab="x axis of mu",type="l")
plot(t,y[t,2],xlab="estimating time",ylab="y axis of mu",type="l")


# M-H Markov chain without burn-in period
par(mfrow=c(3,2))
t <- c()
for(i in 1:3000){t<-c(t,i)}

y <- MH(3000,c(mean(x[,1]),mean(x[,2])),1)
plot(t,y[t,1],xlab="estimating time",ylab="x axis of mu",type="l")
plot(t,y[t,2],xlab="estimating time",ylab="y axis of mu",type="l")

y <- MH(3000,c(mean(x[,1]),mean(x[,2])),0.1)
plot(t,y[t,1],xlab="estimating time",ylab="x axis of mu",type="l")
plot(t,y[t,2],xlab="estimating time",ylab="y axis of mu",type="l")

y <- MH(3000,c(mean(x[,1]),mean(x[,2])),0.01)
plot(t,y[t,1],xlab="estimating time",ylab="x axis of mu",type="l")
plot(t,y[t,2],xlab="estimating time",ylab="y axis of mu",type="l")


# The paths of the M-H Markov chain
par(mfrow=c(3,2))
t <- c()
for(i in 1:3000){t<-c(t,i)}

y <- MH(3000,c(0,0),1)
plot(y[t,1],y[t,2],xlab="x axis of mu",ylab="y axis of mu",type="l",asp=1)
y <- MH(3000,c(mean(x[,1]),mean(x[,2])),1)
plot(y[t,1],y[t,2],xlab="x axis of mu",ylab="y axis of mu",type="l",asp=1)

y <- MH(3000,c(0,0),0.1)
plot(y[t,1],y[t,2],xlab="x axis of mu",ylab="y axis of mu",type="l",asp=1)
y <- MH(3000,c(mean(x[,1]),mean(x[,2])),1)
plot(y[t,1],y[t,2],xlab="x axis of mu",ylab="y axis of mu",type="l",asp=1)

y <- MH(3000,c(0,0),0.01)
plot(y[t,1],y[t,2],xlab="x axis of mu",ylab="y axis of mu",type="l",asp=1)
y <- MH(3000,c(mean(x[,1]),mean(x[,2])),1)
plot(y[t,1],y[t,2],xlab="x axis of mu",ylab="y axis of mu",type="l",asp=1)
