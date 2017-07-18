fmvgamma <- function(p, a)
{
  if (p == 1)
    return(gamma(a))
  else
    return(pi^{(p-1)/2}*gamma(a)*fmvgamma(p-1, a-1/2))
}
  