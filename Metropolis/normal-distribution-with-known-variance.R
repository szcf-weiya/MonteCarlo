s2 = 1; t2 = 10; mu = 5;
y = c(9.37, 10.18, 9.16, 11.60, 10.33)
theta = 0; delta2 = 2; S = 10000; THETA = NULL;
set.seed(1)

for (s in 1:S)
{
  theta.star = rnorm(1, theta, sqrt(delta2))
  log.r = (sum(dnorm(y, theta.star, sqrt(s2), log = TRUE)) +
             dnorm(theta.star, mu, sqrt(t2), log = TRUE)) -
          (sum(dnorm(y, theta, sqrt(s2), log = TRUE)) + 
             dnorm(theta, mu, sqrt(t2), log = TRUE))
  if (log(runif(1)) < log.r){
    theta = theta.star
  }
  THETA = c(THETA, theta)
}

# figure 10.3 
## left
plot(THETA)
## right
hist(THETA, probability = TRUE, breaks = 70, xlim = c(8, 12))
lines(sort(THETA), dnorm(sort(THETA), 10.03, 0.44))