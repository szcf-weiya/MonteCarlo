## R program for pool-adjacent-violators algorithm
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-19

pooladv <- function(f, w)
{
  n = length(f)
  lag = diff(f)
  if (sum(lag < 0) == 0) # f is isotonic
    return(f)
  while (TRUE) {
    idx = which(lag < 0)[1] + 1
    newf = (w[idx]*f[idx] + w[idx-1]*f[idx-1])/(w[idx]+w[idx-1])
    f[idx] = newf
    f[idx-1] = newf
    lag = diff(f)
    if (sum(lag < 0) == 0)
      return(f)
  }
}

## example
pooladv(c(23, 27, 25, 28), rep(1, 4))
# ans: 23, 26, 26, 28