logdnorm <- function(y, mu, sd){
    return( -0.5 * log(2*pi) - 1.0 * log(sd) - 1/(2*sd^2) * (y - mu)^2 )
}

toy.mh <- function(y, S)
{
    s2 = 1; t2 = 10; mu = 5;
    theta = 0; delta2 = 2; S = 10000; THETA = NULL;
    for (s in 1:S)
    {
        theta_star = rnorm(1) * sqrt(delta2) + theta
        log_r = ( sum(logdnorm(y, theta_star, sqrt(s2))) + logdnorm(theta_star, mu, sqrt(t2)) ) - 
                     ( sum(logdnorm(y, theta, sqrt(s2))) + logdnorm(theta, mu, sqrt(t2)) )
        if (log(runif(1)) < log_r){
            theta = theta_star
        }
        THETA = c(THETA, theta)
    }
    return(THETA)
}

# run
#y = c(9.37, 10.18, 9.16, 11.60, 10.33)
y = rnorm(100) * sqrt(10)
#res = toy.mh(y, 1e10)
png("hist-r.png")
hist(toy.mh(y, 1e10))
dev.off()
