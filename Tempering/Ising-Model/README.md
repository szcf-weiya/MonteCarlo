# Ising Model 

## General Chain-Structured Models

There is an important probability distribution used in applications (Liu, 2008):

$$
\pi(\mathbf x) \propto \exp\Big\{-\sum_{i=1}^dh_i(x_{i-1},x_i)\Big\} \tag{*}
$$

where $$\mathbf x = (x_0,x_1,\ldots,x_d)$$. This type of model exhibits a so-called "Markovian structure" because 

$$
\pi(x_i\mid \mathbf x_{[-i]}) \propto \exp\{-h(x_{i-1},x_i)-h(x_i,x_{i+1})\}\,.
$$


