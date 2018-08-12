# Monte Carlo Approximation

## For Bayesian

直观上理解，已知后验分布，但是该分布很难进行积分，而我们想计算均值、概率、分位数等。这时我们可以采用MC近似。也就是从后验分布中产生规模为$n$的样本，

- 用样本的均值代替总体均值
- 用样本的经验分布代替该分布
- 用样本的分位数近似总体的分位数

总结起来，就是先从后验分布中产生样本，然后用样本的量近似总体的量。

```r
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
```

## For Integration

为了计算积分
![](https://latex.codecogs.com/gif.latex?I%20%3D%20%5Cint%20_D%20g%28%5Cmathbf%20x%29d%5Cmathbf%20x)
我们通常采用MC模拟:

1. 计算区域的volume：![](https://latex.codecogs.com/gif.latex?V%20%3D%20%5Cint_D%20d%5Cmathbf%20x)
2. 近似：![](https://latex.codecogs.com/gif.latex?%5Chat%20I_m%3DV%5Cfrac%7B1%7D%7Bm%7D%5Csum%5Climits_%7Bi%3D1%7D%5Emg%28%5Cmathbf%20x%5E%7B%28m%29%7D%29)

根据大数律有

![](https://latex.codecogs.com/gif.latex?%5Clim_%7Bm%5Crightarrow%20%5Cinfty%7D%20%5Chat%20I_m%3DI)

并且由中心极限定理有

![](https://latex.codecogs.com/gif.latex?%5Cfrac%7B1%7D%7BV%7D%7B%7D%5Csqrt%7Bm%7D%28%5Chat%20I_m-I%29%5Crightarrow%20N%280%2C%20%5Csigma%5E2%29)

其中![](https://latex.codecogs.com/gif.latex?%5Csigma%5E2%3Dvar%28g%28%5Cmathbf%20x%29%29)

```r
## estimate pi
##
## I = \int H(x, y)dxdy
## where H(x, y) = 1 when x^2+y^2 <= 1;
## otherwise H(x, y) = 0

## volume
V = 4

n = 100000
x = runif(n, -1, 1)
y = runif(n, -1, 1)
H = x^2+y^2
H[H<=1] = 1
H[H>1] = 0
I = V* mean(H)
cat("I = ", I, "\n")

## n = 100, I = 2.96
## n = 1000, I = 3.22
## n = 10000, I = 3.1536
## n = 100000, I = 3.14504
```
