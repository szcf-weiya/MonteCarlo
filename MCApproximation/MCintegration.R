## estimate pi
##
## I = \int H(x, y)dxdy
## where H(x, y) = 1 when x^2+y^2 <= 1;
## otherwise H(x, y) = 0

## volume
V = 4

n = 100000
x = runif(n, -1, 1)
y = runif(n, -1, 1)
H = x^2+y^2
H[H<=1] = 1
H[H>1] = 0
I = V* mean(H)
cat("I = ", I, "\n")

## n = 100, I = 2.96
## n = 1000, I = 3.22
## n = 10000, I = 3.1536
## n = 100000, I = 3.14504
