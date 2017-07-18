## simulate wishart distribution
y = matrix(c(1, 1, -1, -1, 2, 2, -2, -2, NA, NA, NA, NA,
             1, -1, 1, -1, NA, NA, NA, NA, 2, 2, -2, -2), byrow = T, ncol = 12)
y4 = y[, 1:4]
S4 = matrix(nrow = 2, ncol = 2)
S4[1, 1] = sum(y4[1, ]^2)
S4[1, 2] = sum(y4[1, ]*y4[2, ])
S4[2, 1] = sum(y4[2, ]*y4[1, ])
S4[2, 2] = sum(y4[2, ]*y4[2, ])
S = S4
t = 4

## ###########################
## |    | yobs1 | ymis2 |
## |    | ymis1 | yobs2 | 
## ###########################

yobs1 = y[1, 5:8]
yobs2 = y[2, 9:12]
m = 5000
Z = rWishart(m, t, S)
w = numeric(m)
pi.fWishart = numeric(m)
rho = numeric(m)
for (i in 1:m)
{
  sigma = solve(Z[,,i])
  rho[i] = sigma[1,2]/(sqrt(sigma[1,1]*sigma[2,2]))
  rs.sigma12 = rho[i]*sqrt(sigma[1,1]/sigma[2,2])
  ymis.sigma = sqrt((1-rho[i]^2)*sigma[1,1])
  ymis1 = rnorm(4, yobs1*rs.sigma12, ymis.sigma)
  ymis2 = rnorm(4, yobs2*rs.sigma12, ymis.sigma)
yy = y
yy[2, 5:8] = ymis1
yy[1, 9:12] = ymis2
S12 = matrix(nrow = 2, ncol = 2)
S12[1, 1] = sum(yy[1, ]^2)
S12[1, 2] = sum(yy[1, ]*yy[2, ])
S12[2, 1] = sum(yy[2, ]*yy[1, ])
S12[2, 2] = sum(yy[2, ]*yy[2, ])
 p.fWishart = fWishart(Z[,,i], 12, S12)
 g.fWishart = fWishart(Z[,,i], 4, S4)
  # g.norm = numeric(8)
  # g.norm[1:4] = sapply(1:4, function(ii) fnorm(ymis1[ii], yobs1[ii]*rs.sigma12, ymis.sigma))
  # g.norm[5:8] = sapply(1:4, function(ii) fnorm(ymis2[ii], yobs2[ii]*rs.sigma12, ymis.sigma))
  # 
  # pi.norm = numeric(8)
  # pi.norm[1:4] = sapply(1:4, function(ii) fnorm(ymis1[ii], 0, sqrt(sigma[1,1])))
  # pi.norm[5:8] = sapply(1:4, function(ii) fnorm(ymis2[ii], 0, sqrt(sigma[2,2])))
  #wei = pi.norm/g.norm
  #wei[wei > 10] = 10
  #wei[wei < 0.01] = 0.01
  #w[i] = mean(pi.norm/g.norm)
  #w[i] = prod(pi.norm)/prod(g.norm)
  w[i] = p.fWishart/g.fWishart
  #pi.fWishart[i] = p.fWishart
}

#sigma.post = sum(pi.fWishart*w)/sum(m)

fWishart <- function(X, t, S)
{
  d = dim(S)[1]
  return(det(X)^{(t-d-1)/2}*exp(-0.5*sum(diag(X%*%S)))/(2^{t*d/2}*det(solve(S))^{t/2}*fmvgamma(d, t/2)))
}

fnorm <- function(x, mu, sigma)
{
  return(1/sigma*exp(-0.5*x^2/sigma))
}

# analytical form
rho = seq(-1, 1, by = 0.001)
p = (1-rho^2)^4.5/(1.25-rho^2)^8
plot(rho, p, type = 'l')