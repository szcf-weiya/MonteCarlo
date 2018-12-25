# #######################################
# Example 9.8 and 10.24 Grouped Multinomial Data of
# Robert, C. P., & Casella, G. (2010). Monte Carlo statistical methods (2. ed., softcover reprint of the hardcover 2. ed. 2004). New York, NY: Springer.
# #######################################

using Distributions
using Plots
using LaTeXStrings

# configure
a = [0.06, 0.14, 0.11, 0.09]
b = [0.17, 0.24, 0.19, 0.20]
x = [9, 15, 12, 7, 8]

# normal gibbs
function gibbs(M = 1000)
    Z = ones(M, 4)
    Params = ones(M, 2)
    Eparams = ones(M, 2)
    # prior
    alpha = [1/2, 1/2, 1/2]
    # initial 
    mu = 0.1
    eta = 0.1
    for m in 1:M
        # sample zi
        zi = ones(4)
        for i = 1:2
            aimu = a[i] * mu
            zi[i] = rand(Binomial(x[i], aimu / (aimu + b[i])))
        end
        for i = 3:4
            aieta = a[i] * eta
            zi[i] = rand(Binomial(x[i], aieta / (aieta + b[i])))
        end
        # sample mu eta 
        mu, eta, = rand(Dirichlet([ zi[1] + zi[2] + alpha[1], zi[3] + zi[4] + alpha[2], x[5] + alpha[3] ]))
        # store
        Z[m, :] = zi
        Params[m, :] = [mu, eta]
        #Eparams[m, :] = mean(Params[1:m, :], dims = 1)
    end
    Eparams[:, 1] = cumsum(Params[:, 1]) ./ Array(1:M)
    Eparams[:, 2] = cumsum(Params[:, 2]) ./ Array(1:M)
    return Z, Params, Eparams
end
Z, Params, Eparams = gibbs(10000)

plot(Eparams[:,1], label = L"\mu")
plot!(Eparams[:,2], label = L"\eta")
savefig("gibbs.png")

# Metropolization
function Metrogibbs(M = 1000)
    Z = ones(Int, M, 4)
    Params = ones(M, 2)
    Eparams = ones(M, 2)
    # prior
    alpha = [1/2, 1/2, 1/2]
    # initial 
    mu = 0.1
    eta = 0.1
    for m in 1:M
        # sample zi
        zi = ones(Int, 4)
        for i = 1:4
            if i == 1 || i == 2
                aimu = a[i] * mu
            else
                aimu = a[i] * eta
            end
            p1 = aimu / (aimu + b[i])
            p2 = 1 - p1 
            while true
                zi[i] = rand(Binomial(x[i], p1))
                if m == 1 || zi[i] != Z[m-1, i] # short-circuiting boolean 
                    break
                end
            end
            # accept or not 
            if m == 1
                rho = 1.1
            else
                num1 = 1 - binomial(Z[m-1, i], x[i]) * p1^Z[m-1,i] * p2^(x[i] - Z[m-1, i])
                num2 = 1 - binomial(zi[i], x[i]) * p1^zi[i] * p2^(x[i] - zi[i])
                rho = num1 / num2
            end
            if rand() > rho
                zi[i] = Z[m-1, i]
            end
            
        end
        # sample mu eta 
        mu, eta, = rand(Dirichlet([ zi[1] + zi[2] + alpha[1], zi[3] + zi[4] + alpha[2], x[5] + alpha[3] ]))
        # store
        Z[m, :] = zi
        Params[m, :] = [mu, eta]
        #Eparams[m, :] = mean(Params[1:m, :], dims = 1)
    end
    Eparams[:, 1] = cumsum(Params[:, 1]) ./ Array(1:M)
    Eparams[:, 2] = cumsum(Params[:, 2]) ./ Array(1:M)
    return Z, Params, Eparams
end

Z, Params, Eparams = Metrogibbs(10000)
plot(Eparams[:,1], label = L"\mu")
plot!(Eparams[:,2], label = L"\eta")
savefig("metrogibbs.png")