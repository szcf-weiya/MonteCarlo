## fit logistic regression

using DataFrames, GLM, Plots
temp = [53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70, 70, 70, 72, 73, 75, 75, 76, 76, 78, 79, 81]
failure = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]

data = DataFrame(temp = temp, failure = failure)

logit_fit = glm(@formula(failure ~ temp), data, Binomial(), LogitLink())
#plot(temp, predict(logit_fit), legend = false, xlabel = "Temperature", ylab = "Probability")
#scatter!(temp, predict(logit_fit))

## metropolis-hastings
using Distributions
γ = 0.57721
function ll(α::Float64, β::Float64)
    a = exp.(α .+ β*temp)
    return prod( (a ./ (1 .+ a) ).^failure .* (1 ./ (1 .+ a)).^(1 .- failure) )
end
function mh_logit(T::Int, α_hat::Float64, β_hat::Float64, σ_hat::Float64)
    φ = Normal(β_hat, σ_hat)
    π = Exponential(exp(α_hat+γ))
    Α = ones(T)
    Β = ones(T)
    for t = 1:T-1
        α = log(rand(π))
        β = rand(φ)
        r = ( ll(α, β) / ll(Α[t], Β[t]) ) * ( pdf(φ, Β[t]) / pdf(φ, β) )
        if rand() < r
            Α[t+1] = α
            Β[t+1] = β
        else
            Α[t+1] = Α[t]
            Β[t+1] = Β[t]
        end
    end
    return Α, Β
end

# trace plot

Α, Β = mh_logit(10000, 15.04, -0.233, 0.108)

p1 = plot(Α, legend = false, xlab = "Intercept")
hline!([15.04])

p2 = plot(Β, legend = false, xlab = "Slope")
hline!([-0.233])

plot(p1, p2, layout = (1,2))

# mean trace plot

Αmean = cumsum(Α) ./ collect(1:length(Α))
Βmean = cumsum(Β) ./ collect(1:length(Β))

p1 = plot(Αmean, legend = false, xlab = "Intercept")
hline!([15.04])

p2 = plot(Βmean, legend = false, xlab = "Slope")
hline!([-0.233])

plot(p1, p2, layout = (1,2))
