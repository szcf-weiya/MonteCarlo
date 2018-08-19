## Julia program for tree ordering algorithm
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-19

function treeordering(f::Array, w::Array)
    n = size(f, 1)
    lag = diff(f)
    # if f is isotonic
    if all(i -> i >= 0, lag)
        return(f)
    end
    for j = 1:n
        Aj = sum(f[1:j].*w[1:j])/sum(w[1:j])
        if Aj < f[j+1]
            return([Aj*ones(j); f[(j+1):end]])
        end
    end
end

## example
treeordering([18,17,12,21,16], [1,1,3,3,1])
# ans: 14.2 14.2 14.2 21.0 16.0