# Monte Carlo Integration

## Toy Example

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
