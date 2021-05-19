function exp_rand()
    # cdf of Possion+(log(2)), since P(X=0) = 1/2
    # for positive Possion, q[k] = sum(log(2)^k / k!)
    q = [
    0.6931471805599453,
    0.9333736875190459,
    0.9888777961838675,
    0.9984959252914960,
    0.9998292811061389,
    0.9999833164100727,
    0.9999985691438767,
    0.9999998906925558,
    0.9999999924734159,
    0.9999999995283275,
    0.9999999999728814,
    0.9999999999985598,
    0.9999999999999289,
    0.9999999999999968,
    0.9999999999999999,
    1.0000000000000000
    ]
    a = 0.0
    u = rand()
    while (u <= 0) | (u >= 1)
        u = rand()
    end
    while true
        u += u;
        if (u > 1.0)
            break
        end
        a += q[1]  # μ = q[1], and the times of loop is M, so a = μM
        # TODO: but why the times of loop follows a geometric (1/2) ??
    end
    u -= 1.0 # an independent uniform variate
    if u <= q[1]
        return a + u
    end
    i = 1 # one-based
    ustar = rand()
    umin = ustar
    while true
        ustar = rand()
        if (umin > ustar)
            umin = ustar
        end
        i += 1
        if u <= q[i]
            return a + umin * q[1]
        end
    end
end
