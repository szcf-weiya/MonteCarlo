# State Space Model

Let us take a simple example to illustrate how SMCTC library implement this kind of problem.

The state vector ![](https://latex.codecogs.com/gif.latex?X_n) contains the position and velocity of an object moving in a plane:

![](https://latex.codecogs.com/gif.latex?X_n%3D%28s_n%5Ex%2Cu_n%5Ex%2Cs_n%5Ey%2Cu_n%5Ey%29)

The observation of the position is possible at each time instance. The state and observation equations are linear with additive noise:

![](https://latex.codecogs.com/gif.latex?%5Cbegin%7Barray%7D%7Bccc%7D%20X_n%20%26%3D%26%20AX_%7Bn-1%7D%20&plus;V_n%5C%5C%20Y_n%20%26%3D%26%20BX_n%20&plus;%20%5Calpha%20W_n%20%5Cend%7Barray%7D)

where

![](https://latex.codecogs.com/gif.latex?A%20%3D%20%5Cleft%5B%20%5Cbegin%7Barray%7D%7Bcccc%7D%201%20%26%20%5CDelta%20%26%200%20%26%200%5C%5C%200%20%26%201%20%26%200%20%26%200%5C%5C%200%20%26%200%20%26%201%20%26%20%5CDelta%5C%5C%200%20%26%200%20%26%200%20%26%201%20%5Cend%7Barray%7D%20%5Cright%5D%20%5Cqquad%20B%20%3D%20%5Cleft%5B%20%5Cbegin%7Barray%7D%7Bcccc%7D%201%20%26%200%20%26%200%20%260%5C%5C%200%20%26%200%20%26%201%20%26%200%20%5Cend%7Barray%7D%20%5Cright%20%5D%20%5Cqquad%20%5Calpha%20%3D%200.1)

And assuming that the elements of ![](https://latex.codecogs.com/gif.latex?V_n) are independent normal with variances 0.02 and 0.001 for position and velocity components, respectively. The observation noise, ![](https://latex.codecogs.com/gif.latex?W_n), comprise independent, identically distributed t-distributed random variables with ![](https://latex.codecogs.com/gif.latex?%5Cnu%3D10) degrees of freedom.


The observations are as follows:

![](observation.png)

Apply the SMC algorithm to the SSM, and calculate the mean of position at each iteration, we can obtain the following results:

![](res.png)
