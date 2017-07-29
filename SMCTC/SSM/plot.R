## ########################################
## pre-processing
##
## sed '1d' data.csv > data2.csv
## if want to edit directly in data.csv
## just
##    sed -i '1d' data.csv
## ########################################

data <- read.table('data2.csv')

png('observation.png')
plot(data, xlab = 'x', ylab = 'y')


## results

res <- read.csv('res.csv', header = FALSE)
xm <- res$V1
ym <- res$V2
xv <- res$V3
yv <- res$V4

png('res.png')
plot(xm, ym, type = 'o')
## 95% intervals
points(xm, ym-1.96*sqrt(yv), col = 'red', type = 'o')
points(xm, ym+1.96*sqrt(yv), col = 'blue', type = 'o')
