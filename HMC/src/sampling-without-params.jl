using CmdStan
src = "
data {
    real<lower=0, upper=1> theta;
    int<lower=0> K;
    int<lower=0> N;
}
model {

}
generated quantities {
    int<lower=0, upper=K> y[N];
    for (n in 1:N)
        y[n] = binomial_rng(K, theta);
}
"
model = Stanmodel(method = Sample(algorithm=fixed_param), model = src)
data = Dict("theta" => 0.5, "K" => 100, "N" => 50)
stan(model, data)