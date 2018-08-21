## Julia program for Calculation of Gamma function
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-21
##
## refer to "Coefficients for the Lanczos Approximation to the Gamma Function"
##   https://mrob.com/pub/ries/lanczos-gamma.html

const LG_g = 5.0;
const LG_N = 6;
const lct = [
    1.000000000190015,
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5];
const ln_sqrt_2_pi = 0.91893853320467274178;
const g_pi = 3.14159265358979323846;

function lanczos_ln_gamma(z)
    if z < 0.5
        #Use Euler's reflection formula:
        #Gamma(z) = Pi / [Sin[Pi*z] * Gamma[1-z]];
        return log(g_pi / sin(g_pi * z)) - lanczos_ln_gamma(1.0 - z);
    end
    z = z - 1.0;
    base = z + LG_g + 0.5;  # Base of the Lanczos exponential
    sum = 0;
    ## We start with the terms that have the smallest coefficients and largest denominator.
    for i = 2:LG_N+1
      sum += lct[i] / (z + i - 1);
    end
    sum += lct[1];
    ## This printf is just for debugging
#    @printf("ls2p %7g  l(b^e) %7g   -b %7g  l(s) %7g\n", ln_sqrt_2_pi,
#              log(base)*(z+0.5), -base, log(sum));
    ## Gamma[z] = Sqrt(2*Pi) * sum * base^[z + 0.5] / E^base
    return ((ln_sqrt_2_pi + log(sum)) - base) + log(base)*(z+0.5);
end

function lanczos_gamma(z)
    return(exp(lanczos_ln_gamma(z)));
end

## example
lanczos_gamma(1)