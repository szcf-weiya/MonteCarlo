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

k <- 3
s0 <- c(1851, 1890.5, 1927, 1963)
h0 <- c(1.7,1.7,1.7)
X<-Move(1000,h0,s0)
par(mfrow=c(2,1))
plot(X[,1],ylim=c(0,4),xlab="Estimating Time",ylab="Height",type="l",col="red",main="Paths of Simulation for Heights")
lines(X[,2],xlab="Estimating Time",ylab="Height",type="l",col="blue")
lines(X[,3],xlab="Estimating Time",ylab="Height",type="l",col="green")
plot(density(X[,1]),xlim=c(0,4),ylim=c(0,3),col="red",main="Densities for Heights")
lines(density(X[,2]),col="blue")
lines(density(X[,3]),col="green")

par(mfrow=c(2,1))
hist(X[,5])
plot(density(X[,5]),xlim=c(1890,1893),col="red")
par(mfrow=c(2,1))
hist(X[,6])
plot(density(X[,6]),xlim=c(1945,1955),col="blue")
pdf("10000position3.pdf",height=10,width=6)
par(mfrow=c(2,1))
plot(X[,5],ylim=c(1851,1963),xlab="Estimating Time",ylab="Position",type="l",col="red",main="Paths of Simulation for Positions")
lines(X[,6],xlab="Estimating Time",ylab="Position",type="l",col="blue")
plot(density(X[,5]),xlim=c(1851,1963),col="red",main="Densities for Positions")
lines(density(X[,6]),col="blue")

X<-Move(10000,h0,s0)
...
X<-Move(100000,h0,s0)
...


k <- 2
s0 <- c(1851, 1907, 1963)
h0 <- c(1.7,1.7)
X<-Move(100000,h0,s0)
par(mfrow=c(2,1))
plot(X[,1],ylim=c(0,4),xlab="Estimating Time",ylab="Height",type="l",col="red",main="Paths of Simulation for Heights")
lines(X[,2],xlab="Estimating Time",ylab="Height",type="l",col="blue")
plot(density(X[,1]),xlim=c(0,4),ylim=c(0,5),col="red",main="Densities for Heights")
lines(density(X[,2]),col="blue")

par(mfrow=c(2,1))
plot(X[,4],ylim=c(1851,1963),xlab="Estimating Time",ylab="Position",type="l",col="red",main="Paths of Simulation for Positions")
plot(density(X[,4]),xlim=c(1851,1963),col="red",main="Density for Position")


k <- 4
s0 <- c(1851, 1892, 1910, 1948, 1963)
h0 <- c(1.7,1.7,1.7,1.7)
X<-Move(100000,h0,s0)
par(mfrow=c(2,1))
plot(X[,1],ylim=c(0,4),xlab="Estimating Time",ylab="Height",type="l",col="red",main="Paths of Simulation for Heights")
lines(X[,2],xlab="Estimating Time",ylab="Height",type="l",col="blue")
lines(X[,3],xlab="Estimating Time",ylab="Height",type="l",col="green")
lines(X[,4],xlab="Estimating Time",ylab="Height",type="l",col="yellow")
plot(density(X[,1]),xlim=c(0,4),ylim=c(0,3),col="red",main="Densities for Heights")
lines(density(X[,2]),col="blue")
lines(density(X[,3]),col="green")
lines(density(X[,4]),col="yellow")

par(mfrow=c(2,1))
plot(X[,6],ylim=c(1851,1963),xlab="Estimating Time",ylab="Position",type="l",col="red",main="Paths of Simulation for Positions")
lines(X[,7],xlab="Estimating Time",ylab="Position",type="l",col="blue")
lines(X[,8],xlab="Estimating Time",ylab="Position",type="l",col="green")
plot(density(X[,6]),xlim=c(1851,1963),col="red",main="Densities for Positions")
lines(density(X[,7]),col="blue")
lines(density(X[,8]),col="green")
