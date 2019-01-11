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

library("boot")

# disaster number function m(.)
count.disasters <- function(s1, s2) {
        return(sum(coal >= s1 & coal < s2))
}

# value of Gamma(a,b)
a=1
b=200/365.24

# accept probability function for Height move
alpha.hmove <- function(s,h,h_tilde,j) {
        m <- count.disasters(s[j],s[j+1])
        log_likelihood_ratio <- (h[j]-h_tilde)*(s[j+1]-s[j])+m*(log(h_tilde)-log(h[j]))
        log_prior_ratio <- a*(log(h_tilde)-log(h[j]))-b*(h_tilde-h[j])
        log_pi_ratio <- log_likelihood_ratio + log_prior_ratio
        return(exp(log_pi_ratio))
}

# accept probability function for Position move
alpha.smove <- function(s,h,s_tilde,j) {
              m1 <- count.disasters(s[j-1],s[j])
        m1_tilde <- count.disasters(s[j-1],s_tilde)
              m2 <- count.disasters(s[j],s[j+1])
        m2_tilde <- count.disasters(s_tilde,s[j+1])
        log_likelihood_ratio <- (h[j]-h[j-1])*(s_tilde-s[j])+(m1_tilde - m1)*log(h[j-1]) + (m2_tilde-m2) * log(h[j])
        log_prior_ratio <- log(s_tilde-s[j-1])+log(s[j+1]-s_tilde)-log(s[j]-s[j-1])-log(s[j+1] - s[j])
        log_pi_ratio <- log_likelihood_ratio + log_prior_ratio
        return(exp(log_pi_ratio))
}

# main program
Move <- function(n,h_initial,s_initial){
    H <-matrix(NA,nrow=n,ncol=k)
    S <- matrix(NA,nrow=n,ncol=k+1)
    h <- h_initial
        s <- s_initial

    for (i in 1:n){
        j <- sample(1:(2*k-1),size=1,replace=TRUE)

                if (j <= k) {                       #from 1 to k is for Height
                    u <- runif(1,-0.5,0.5)
            h_tilde <- h[j] * exp(u)
            U <- runif (1)
            if (U < alpha.hmove(s,h,h_tilde,j)){
                h[j] <- h_tilde
            }
            }

                if (j > k) {                        #from k+1 to 2k-1 is for Position
                    j = j - k + 1
            s_tilde <- runif(1,s[j-1],s[j+1])
            U <- runif (1)
            if (U < alpha.smove(s,h,s_tilde,j)){
                s[j] <- s_tilde
            }
                }

                H[i,] <- h
        S[i,] <- s
        }
        return(cbind(H,S))
}
