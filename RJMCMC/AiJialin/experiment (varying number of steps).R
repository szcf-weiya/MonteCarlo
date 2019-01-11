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

k0 <- 1
s0 <- c(1851,1963)
h0 <- c(1.7)

lambda = 3

# value of Gamma(a,b)
a=1
b=200/365.24


# accept probability function for Height move

alpha.hmove <- function(s,h,h_tilde,j) {
        m <- count.disasters(s[j],s[j+1])
        log_likelihood_ratio <- (h[j] - h_tilde) * (s[j+1] - s[j]) + m * (log(h_tilde) - log(h[j]))
        log_prior_ratio <- a * (log(h_tilde) - log(h[j])) - b * (h_tilde - h[j])
        log_pi_ratio <- log_likelihood_ratio + log_prior_ratio
        return(exp(log_pi_ratio))

}


# accept probability function for Position move

alpha.smove <- function(s,h,s_tilde,j) {
              m1 <- count.disasters(s[j-1],s[j])
        m1_tilde <- count.disasters(s[j-1],s_tilde)
              m2 <- count.disasters(s[j],s[j+1])
        m2_tilde <- count.disasters(s_tilde,s[j+1])
        log_likelihood_ratio <- (h[j] - h[j-1]) * (s_tilde - s[j]) + (m1_tilde - m1) * log(h[j-1]) + (m2_tilde - m2) * log(h[j])
        log_prior_ratio <- log(s_tilde - s[j-1]) + log(s[j+1] - s_tilde) - log(s[j] - s[j-1]) - log(s[j+1] - s[j])
        log_pi_ratio <- log_likelihood_ratio + log_prior_ratio
        return(exp(log_pi_ratio))

}


# accept probability function for birth of step

alpha.birth <- function(s,h,s_star,h1_prime,h2_prime,j,k) {
        m1_prime <- count.disasters(s[j],s_star)
        m2_prime <- count.disasters(s_star,s[j+1])
                mj <- m1_prime + m2_prime
        log_likelihood_ratio <- (
           - h1_prime * (s_star - s[j]) - h2_prime * (s[j+1] - s_star)
           + m1_prime * log(h1_prime) + m2_prime * log(h2_prime)
           + h[j] * (s[j+1] - s[j]) - mj * log(h[j])
                )
        log_prior_ratio <- (
                   log(2) + log(lambda) + log(2*k+1)
                   - log(s[k+1] - s[1]) + log(s_star - s[j]) + log (s[j+1] - s_star) - log(s[j+1] - s[j])
                   + a * log(b) - log(gamma(a)) + (a-1) * (log(h1_prime) + log(h2_prime)) - a * log(h[j])
                   - b * (h1_prime + h2_prime - h[j])
                )
                log_proposal_ratio <- log(death(k)) - log(birth(k-1)) - log(k)
                log_Jacobian <- 2 * log(h1_prime + h2_prime)
        log_pi_ratio <- log_likelihood_ratio + log_prior_ratio + log_proposal_ratio + log_Jacobian
        return(exp(log_pi_ratio))
}



# accept probability function for death of step

alpha.death <- function(s,h,hj_prime,j,k) {
        m1 <- count.disasters(s[j],s[j+1])
        m2 <- count.disasters(s[j+1],s[j+2])
        mj_prime <- m1 + m2
        log_likelihood_ratio <- (
                   - hj_prime * (s[j+2] - s[j]) + mj_prime * log(hj_prime)
                   + h[j] * (s[j+1] - s[j]) + h[j+1] * (s[j+2] - s[j+1])
           - m1 * log(h[j]) - m2 * log(h[j+1])
                )
        log_prior_ratio <- (
                - log(2) - log (lambda) - log(2*k-1)
                + log(s[k+1] - s[1]) - log(s[j+1] - s[j]) - log (s[j+2] - s[j+1]) + log(s[j+2] - s[j])
                - a * log(b) + log(gamma(a)) - (a-1) * (log(h[j]) + log(h[j+1])) + a * log(hj_prime)
                + b * (h[j] + h[j+1] - hj_prime)
                )
                log_proposal_ratio <- + log(birth(k-2)) + log(k-1) -log(death(k-1))
                log_Jacobian <- - 2 * log(h[j] + h[j+1])
        log_pi_ratio <- log_likelihood_ratio + log_prior_ratio + log_proposal_ratio + log_Jacobian
        return(exp(log_pi_ratio))

}


birth <- function(change_points){return(3.6/7 * min(1,lambda/(change_points+1)))}
death <- function(change_points){return(3.6/7 * min(1,change_points/lambda))}

Move <- function(n, h, s, k) {
    H <- list()
    S <- list()
        K <- vector(length=n)
        K[1] <- k
    H[[1]] <- h
        S[[1]] <- s

    for (i in 2:n){
        position_prob <- ifelse(k<=1, 0, 0.5 * (1 - birth(k-1) - death(k-1)))
                #if (k<=1) {position_prob <- 0} else {position_prob <- 0.5 * (1 - birth(k-1) - death(k-1))}

                height_prob <- 1 - birth(k-1) - death(k-1) - position_prob

        type <- runif(1)

                if (type>1-height_prob) {                       #from 1 to k is for Height h_1 to h_k
                j <- sample(1:k,size=1)
                u <- runif(1,-0.5,0.5)
            h_tilde <- h[j] * exp(u)
            U <- runif (1)
            if (U < alpha.hmove(s,h,h_tilde,j)) {
                h[j] <- h_tilde
            }
            }

                if (type<=1-height_prob && type>1-height_prob-position_prob) {                     #from 2 to k is for Position s_2 to s_k
                    j <- sample(1:(k-1), size=1) + 1
            s_tilde <- runif(1,s[j-1],s[j+1])
                        U <- runif (1)
            if (U < alpha.smove(s,h,s_tilde,j)) {
                s[j] <- s_tilde
            }
                }

                if (type>=birth(k-1) && type<=birth(k-1)+death(k-1)) {                    #from 1 to k-1 is for death of steps d_2 to d_k
                    j <- sample(1:(k-1), size=1)
                        r <- (s[j+2] - s[j+1]) / (s[j+2] - s[j])
                        hj_prime <- h[j]^(1-r) * h[j+1]^r                    #exp((1-r) * log(h[j]) + r * log(h[j+1]))
                        U <- runif(1)
                        if (U < alpha.death(s,h,hj_prime,j,k)){
                            k <- k - 1
                            h[j] <- hj_prime
                h <- h[-(j+1)]
                                s <- s[-(j+1)]
                        }
                }

                if (type<=birth(k-1)){                                             #from 1 to k is for birth of steps b_1 to b_k
                    j <- sample(1:k,size=1)
                        s_star <- runif(1,s[j],s[j+1])
                        u <- runif(1)
                        r <- exp(log(s[j+1] - s_star) - log(s[j+1] - s[j]))
                        h1_prime = h[j] * exp(r * (log(u) - log(1-u)))
                        h2_prime = h[j] * exp((1-r) * (log(1-u) - log(u)))
                        U <- runif(1)
                        if (U < alpha.birth(s,h,s_star,h1_prime,h2_prime,j,k)){
                                s <- c(s[1:j], s_star, s[(j+1):(k+1)])
                                if (j > 1) {
                                        left <- h[1:(j-1)]
                                } else {
                                        left <- c()
                                }
                                if (j < k) {
                                        right <- h[(j+1):k]
                                } else {
                                        right <- c()
                                }
                                h <- c(left, h1_prime, h2_prime, right)
                                k <- k + 1
                        }
                }

                K[i] <- k
                H[[i]] <- h
                S[[i]] <- s
        }

        return(list(h=H,k=K,s=S))
}
