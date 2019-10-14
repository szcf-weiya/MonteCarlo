# Simulation for a Mixture of Exponential Distributions

Suppose we want to simulate pseudo random numbers from the following distribution, which I came across in [r - Finding a way to simulate random numbers for this distribution - Cross Validated](https://stats.stackexchange.com/questions/391054/finding-a-way-to-simulate-random-numbers-for-this-distribution),

$$
F(x) = 1-\exp\left(-ax - \frac{b}{p+1}x^{p+1}\right)\,\quad x\ge 0
$$

where $$a, b>0$$ and $$p\in (0,1)$$.

[Xi'an](https://stats.stackexchange.com/users/7224/xian) provided a very elegant solution.


Note that

$$
(1-F(x))=\exp\left\{-ax-\frac{b}{p+1}x^{p+1}\right\}=\underbrace{\exp\left\{-ax\right\}}_{1-F_1(x)}\underbrace{\exp\left\{-\frac{b}{p+1}x^{p+1}\right\}}_{1-F_2(x)}
$$

the distribution $$F$$ is the distribution of

$$
X=\min\{X_1,X_2\}\qquad X_1\sim F_1\,,X_2\sim F_2
$$

since

$$
F(x) = P(X\le x) = 1- P(X > x)  = 1-P(X_1 > x)P(X_2 > x) = 1-(1-F_1(X))(1-F_2(X))\,.
$$

Thus, the R code to simulate is simple to be

```r
x = pmin(rexp(n, a), rexp(n, b/(p+1))^(1/(p+1)))
```

## References

- [Finding a way to simulate random numbers for this distribution](https://stats.stackexchange.com/questions/391054/finding-a-way-to-simulate-random-numbers-for-this-distribution)
- [simulation fodder for future exams](https://xianblog.wordpress.com/2019/02/20/simulation-fodder-for-future-exams/)
