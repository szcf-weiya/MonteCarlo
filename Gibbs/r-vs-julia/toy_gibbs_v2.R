bigibbs <- function(T, rho)
{
    x = numeric(T+1)
    y = numeric(T+1)
    for (t in 1:T){
        x[t+1] = rnorm(1) * sqrt(1-rho^2) + rho*y[t]
        y[t+1] = rnorm(1) * sqrt(1-rho^2) + rho*x[t+1]
    }
    return(list(x=x, y=y))
}
## example
res = bigibbs(2e6, 0.5)

## plot
# png("hist-gibbs-r.png")
# par(mfrow = c(2, 1))
# hist(res$x)
# hist(res$y)
# dev.off()