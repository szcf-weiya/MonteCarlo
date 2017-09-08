# Monte Carlo Approximation

## For Bayesian

直观上理解，已知后验分布，但是该分布很难进行积分，而我们想计算均值、概率、分位数等。这时我们可以采用MC近似。也就是从后验分布中产生规模为$n$的样本，

- 用样本的均值代替总体均值
- 用样本的经验分布代替该分布
- 用样本的分位数近似总体的分位数

总结起来，就是先从后验分布中产生样本，然后用样本的量近似总体的量。

# For Integration

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
