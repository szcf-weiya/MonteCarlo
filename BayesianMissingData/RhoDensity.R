## simulate wishart distribution
y = matrix(c(1, 1, -1, -1, 2, 2, -2, -2, NA, NA, NA, NA,
             1, -1, 1, -1, NA, NA, NA, NA, 2, 2, -2, -2), byrow = T, ncol = 12)
y4 = y[, 1:4]
#y4 = y[, 1:8]
#y4[2, 5:8] = y[2, 9:12]

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
m = 50000
Z = rWishart(m, t, S)
rho = numeric(m)
for (i in 1:m)
{
  sigma = solve(Z[,,i])
  rho[i] = sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])
}

# analytical form
rho = seq(-1, 1, by = 0.001)
p = (1-rho^2)^4.5/(1.25-rho^2)^8
plot(rho, p, type = 'l')