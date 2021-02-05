# adapt from mcsm's sqar.R (https://github.com/cran/mcsm/blob/master/R/sqar.R)

library(coda)

T = 10^4
set.seed(1234)
rho = 0.85; tau = 0.2
# x minus one; plus one; y current
xm = -0.94; xp = -1.12; yc = 3.17

# target function
ef = function(x) -.5*((xm*rho-x)^2+(x*rho-xp)^2+(yc-x^2)^2/tau^2 )

eef=function(x) exp(ef(x))

# random walk MH
smpl = NULL
for (scale in c(.1, .9)) {
  xmc = sqrt(abs(yc))
  for (t in 2:T) {
    prop = xmc[t-1] + scale*rnorm(1)
    if (log(runif(1)) > ef(prop) - ef(xmc[t-1]))
      prop = xmc[t-1]
    xmc = c(xmc, prop)
  }
  smpl = cbind(smpl, xmc)
}

G = 10
par(mfrow = c(2, 2))
for (i in 1:2) {
  xmc = smpl[, i]
  print(geweke.diag(mcmc(xmc)))
  print(heidel.diag(mcmc(xmc)))
  
  thn = xmc[seq(1, T, by = G)]
  kst = NULL
  for (m in seq(T/(10*G), T/G, le = 100)) kst = c(kst, ks.test(thn[1:(m/2)], thn[(m/2)+(1:(m/2))])$p)
  hist(xmc, pro = T, col = "grey85", nclass = 150, main = "", xlab = "", ylab = "")
  x = seq(min(xmc), max(xmc), le = 200)
  # unnormalized prob.
  uprob = sapply(x, eef)
  # normalized constant
  Z = max(uprob) / max(density(xmc)$y)
  lines(x, uprob/Z, lwd = 2, col = "gold4")
  plot(seq(1, T, le = 100), kst, pch = 19, xlab = "iter", ylab = "p-value")
}