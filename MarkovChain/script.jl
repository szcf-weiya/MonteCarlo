function sumbinom(K::Int)
    res = 0
    for k = 0:K
        res += binomial(K, k)^2
    end
    return res, binomial(2K, K)
end