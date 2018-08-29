higibbs <- function(Y, T, nu0 = 1, s20 = 100, eta0 = 1, t20 = 100, mu0 = 50, g20 = 25)
{
    ## starting values
    m = length(unique(Y[, 1]))
    n = sv = ybar = rep(NA, m)

    for (j in 1:m)
    {
    ybar[j] = mean(Y[Y[, 1]==j, 2])
    sv[j] = var(Y[Y[, 1]==j, 2])
    n[j] = sum(Y[, 1]==j)
    }
    theta = ybar
    sigma2 = mean(sv)
    mu = mean(theta)
    tau2 = var(theta)

    ## setup MCMC
    THETA = matrix(nrow = T, ncol = m)
    SMT = matrix(nrow = T, ncol = 3)

    ## MCMC algorithm
    for (s in 1:T)
    {
        # sample new values of the thetas
        for (j in 1:m)
        {
            vtheta = 1/(n[j]/sigma2+1/tau2)
            etheta = vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
            theta[j] = rnorm(1, etheta, sqrt(vtheta))
        }
        # sample new value of sigma2
        nun = nu0 + sum(n)
        ss = nu0*s20
        for (j in 1:m)
        {
            ss = ss + sum(Y[Y[,1]==j, 2]-theta[j])^2
        }
        sigma2 = 1/rgamma(1, nun/2, ss/2)
        
        # sample new value of mu
        vmu = 1/(m/tau2+1/g20)
        emu = vmu*(m*mean(theta)/tau2 + mu0/g20)
        mu = rnorm(1, emu, sqrt(vmu))
        
        # sample a new value of tau2
        etam = eta0 + m
        ss = eta0*t20 + sum((theta-mu)^2)
        tau2 = 1/rgamma(1, etam/2, ss/2)
        
        # store results
        THETA[s, ] = theta
        SMT[s, ] = c(sigma2, mu, tau2)
    }
    return(list(theta = THETA, smt = SMT))
}

# run
Y = read.table("math-score-Y.csv")
res = higibbs(Y, 5000)
png("hist-mu-r.png")
hist(res$smt[,2])
dev.off()