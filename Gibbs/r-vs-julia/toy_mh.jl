## log density function of N(mu, sd^2) 
function logdnorm(y::Array, mu::Float64, sd::Float64)
    return( -0.5 * log(2*pi) .- 1.0 * log(sd) .- 1/(2*sd^2) * (y .- mu).^2 )
end

function logdnorm(y::Float64, mu::Float64, sd::Float64)
    return( -0.5 * log(2*pi) .- 1.0 * log(sd) .- 1/(2*sd^2) * (y .- mu).^2 )
end
function toy_mh(y::Array, S::Int64)
    s2 = 1.0; t2 = 10.0; mu = 5.0;
    
    theta = 0.0; delta2 = 2; S = 10000; THETA::Array{Float64} = []
    
    for s = 1:S
        theta_star = randn() * sqrt(delta2) + theta
        log_r = ( sum(logdnorm(y, theta_star, sqrt(s2))) + logdnorm(theta_star, mu, sqrt(t2)) ) - 
                     ( sum(logdnorm(y, theta, sqrt(s2))) + logdnorm(theta, mu, sqrt(t2)) )
    
        if log(rand()) < log_r
            theta = theta_star
        end
        THETA = append!(THETA, theta)
    end    
    return(THETA)
end

#y = [9.37, 10.18, 9.16, 11.60, 10.33]
y = randn(100) * sqrt(10)
res = toy_mh(y, Int64(1e10))
# using PyPlot
# plt[:hist](toy_mh(5000))
# show()