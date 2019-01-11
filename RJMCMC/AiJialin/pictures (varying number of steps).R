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

X <- Move(100000, h0, s0, k0)
H <- X$h
K <- X$k
S <- X$s

# the distribution of K (K>=1)
hist(K,breaks=seq(0.5, max(K)+0.5, by=0.5)+0.20)

# grouping data in different dimensions
T2<-list()
T3<-list()
T4<-list()
R2<-list()
R3<-list()
R4<-list()
j2<-0
j3<-0
j4<-0
for (i in 1:length(K)){
if (K[i]==3){
    j3=j3+1
    T3[j3]<-X$h[i]
        R3[j3]<-X$s[i]
        }
else if (K[i]==4){
    j4=j4+1
    T4[j4]<-X$h[i]
        R4[j4]<-X$s[i]
        }
else if (K[i]==2){
    j2=j2+1
    T2[j2]<-X$h[i]
        R2[j2]<-X$s[i]
        }
}
H2 <- matrix(unlist(T2), ncol=2, byrow=TRUE)
H3 <- matrix(unlist(T3), ncol=3, byrow=TRUE)
H4 <- matrix(unlist(T4), ncol=4, byrow=TRUE)
S2 <- matrix(unlist(R2), ncol=3, byrow=TRUE)
S3 <- matrix(unlist(R3), ncol=4, byrow=TRUE)
S4 <- matrix(unlist(R4), ncol=5, byrow=TRUE)


# conditioned on k=3
par(mfrow=c(2,1))
plot(H3[,1],ylim=c(0,4),xlab="Estimating Time",ylab="Height",type="l",col="red",main="Paths of Simulation for Heights")
lines(H3[,2],xlab="Estimating Time",ylab="Height",type="l",col="blue")
lines(H3[,3],xlab="Estimating Time",ylab="Height",type="l",col="green")
plot(density(H3[,1]),xlim=c(0,4),ylim=c(0,3),col="red",main="Densities for Heights")
lines(density(H3[,2]),col="blue")
lines(density(H3[,3]),col="green")

par(mfrow=c(2,1))
plot(S3[,2],ylim=c(1851,1963),xlab="Estimating Time",ylab="Position",type="l",col="red",main="Paths of Simulation for Positions")
lines(S3[,3],xlab="Estimating Time",ylab="Position",type="l",col="blue")
plot(density(S3[,2]),xlim=c(1851,1963),col="red",main="Densities for Positions")
lines(density(S3[,3]),col="blue")


#conditioned on k=2
par(mfrow=c(2,1))
plot(H2[,1],ylim=c(0,4),xlab="Estimating Time",ylab="Height",type="l",col="red",main="Paths of Simulation for Heights")
lines(H2[,2],xlab="Estimating Time",ylab="Height",type="l",col="blue")
plot(density(H2[,1]),xlim=c(0,4),ylim=c(0,3),col="red",main="Densities for Heights")
lines(density(H2[,2]),col="blue")

par(mfrow=c(2,1))
plot(S2[,2],ylim=c(1851,1963),xlab="Estimating Time",ylab="Position",type="l",col="red",main="Paths of Simulation for Positions")
plot(density(S2[,2]),xlim=c(1851,1963),col="red",main="Densities for Positions")


#conditioned on k=4
par(mfrow=c(2,1))
plot(H4[,1],ylim=c(0,4),xlab="Estimating Time",ylab="Height",type="l",col="red",main="Paths of Simulation for Heights")
lines(H4[,2],xlab="Estimating Time",ylab="Height",type="l",col="blue")
lines(H4[,3],xlab="Estimating Time",ylab="Height",type="l",col="green")
lines(H4[,4],xlab="Estimating Time",ylab="Height",type="l",col="yellow")
plot(density(H4[,1]),xlim=c(0,4),ylim=c(0,3),col="red",main="Densities for Heights")
lines(density(H4[,2]),col="blue")
lines(density(H4[,3]),col="green")
lines(density(H4[,4]),col="yellow")

par(mfrow=c(2,1))
plot(S4[,2],ylim=c(1851,1963),xlab="Estimating Time",ylab="Position",type="l",col="red",main="Paths of Simulation for Positions")
lines(S4[,3],xlab="Estimating Time",ylab="Position",type="l",col="blue")
lines(S4[,4],xlab="Estimating Time",ylab="Position",type="l",col="green")
plot(density(S4[,2]),xlim=c(1851,1963),col="red",main="Densities for Positions")
lines(density(S4[,3]),col="blue")
lines(density(S4[,4]),col="green")
