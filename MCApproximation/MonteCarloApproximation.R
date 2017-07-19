## ########################################
## prior distribution: gamma(a, b)
## Y_1, ..., Y_n\mid \theta: iid Poisson(\theta)
## posterior distribution: gamma(a + \sum y_i, b+n)
## ########################################

## ########################################
## Expectation
## posterior mean: (a+\sum y_i)/(b+n)=1.51
## ########################################

a = 2; b = 1
sy = 66; n = 44

theta.mc10 = rgamma(10, a+sy, b+n)
theta.mc100 = rgamma(100, a+sy, b+n)
theta.mc1000 = rgamma(1000, a+sy, b+n)

mean(theta.mc10)
mean(theta.mc100)
mean(theta.mc1000)

## ########################################
## Probabilities
## 
## ########################################

## posterior Probabilities

pgamma(1.75, a+sy, b+n)

## MC approximations
mean(theta.mc10 < 1.75)
mean(theta.mc100 < 1.75)
mean(theta.mc1000 < 1.75)

## ########################################
## quantiles
## 
## ########################################

## posterior quantiles
qgamma(c(.025, .975), a+sy, b+n)

## MC approximations

quantiles(theta.mc10, c(.025, .975))
quantiles(theta.mc100, c(.025, .975))
quantiles(theta.mc1000, c(.025, .975))


## #######################################
## Log-odds
## #######################################

a = 1; b = 1
theta.prior.mc = rbeta(10000, a, b)
gamma.prior.mc = log(theta.prior.mc/(1-theta.prior.mc))

n0 = 860-441; n1 = 441
theta.post.mc = rbeta(10000, a+n1, b+n0)
gamma.post.mc = log(theta.post.mc/(1-theta.post.mc))

## #######################################
## Functions of two parameters
## #######################################

a = 2; b = 1
sy1 = 217; n1 = 111
sy2 = 66; n2 = 44
theta1.mc = rgamma(10000, a+sy1, b+n1)
theta2.mc = rgamma(10000, a+sy2, b+n2)

mean(theta1.mc > theta2.mc)

## ######################################
## posterior 
##
## ######################################

a = 2; b = 1
sy1 = 217; n1 = 111
sy2 = 66; n2 = 44

theta1.mc = rgamma(10000, a+sy1, b+n1)
theta2.mc = rgamma(10000, a+sy2, b+n2)
y1.mc = rpois(10000, theta1.mc)
y2.mc = rpois(10000, theta2.mc)

mean(y1.mc > y2.mc)

## #####################################
## ratio
##
## #####################################
a = 1; b = 2
t.mc = NULL

for (s in 1:10000){
    theta1 = rgamma(1, a+sy1, b+n1)
    y1.mc = rpois(n1, theta1)
    t.mc = c(t.mc, sum(y1.mc==2)/sum(y1.mc==1))
    }
