## Julia program for Metropolization of the Gibbs Sampler
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-25

## example 7.3.8

## data
a = [0.06, 0.14, 0.11, 0.09];
b = [0.17, 0.24, 0.19, 0.20];
x = [9, 15, 12, 7, 8];

function ex738(T)
    z = ones(T+1, 4)
    for t = 1:T
        for i = 1:2
            atrial = 0
            for k = 1:x[i]
                u = rand()
                if u <= a[i]*mu # ? sample mu first
            end
        end
    end
    println(1)
end


