library(RcppSMC)
res <- pfLineartBS(plot=TRUE)
if (interactive()) ## if not running R CMD check etc
     res <- pfLineartBS(onlinePlot=pfLineartBSOnlinePlot)
