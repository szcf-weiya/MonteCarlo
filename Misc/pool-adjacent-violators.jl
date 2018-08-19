## Julia program for pool-adjacent-violators algorithm
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-19


function pooladv(f::Array, w::Array)
   n = size(f, 1)
   lag = diff(f)
   if all(i -> i >= 0, lag)
        return(f)
   end
   while true
        for i = 1:(n-1)
            global idx
            idx = i+1
            lag[i] < 0 && break
        end
        newf = (w[idx]*f[idx] + w[idx-1]*f[idx-1])/(w[idx]+w[idx-1])
        f[idx] = newf
        f[idx-1] = newf
        lag = diff(f)
        if all(i -> i>=0, lag)
            return(f)
        end
   end
end

## example
g = pooladv([23, 27, 25, 28], ones(4))
for i in eachindex(g)
    println(g[i])
end