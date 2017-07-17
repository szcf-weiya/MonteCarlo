N = 150
POS = matrix(c(0, 0,
               0, 1), 2, 2)

growPolymer <- function(N)
{
  for (iter in 1:N)
  {
    n = nrow(POS)
    xt = POS[n,]
    
    xtt = rbind(c(xt[1]+1, xt[2]),
                c(xt[1]-1, xt[2]),
                c(xt[1], xt[2]-1),
                c(xt[1], xt[2]+1))
    
    occupied = c()
    
    for (i in c(1:4))
    {
      if (n == 2)
      {
        flag = sum(xtt[i, ] == POS[-n, ]) 
        flagflag = sum(flag == 2)
      }
      else
      {
        flag = colSums(apply(POS[-n,], 1, function(x) x == xtt[i, ]))
        flagflag = sum(flag == 2)
      }
      if (flagflag != 0)
      {
        occupied = c(occupied, i)
      }
    }
    avail = 4 - length(occupied)
    if (avail == 0)
    {
      cat(paste0("Stop at iter = ", iter)) 
      break
    }
    r = runif(1)
    
    unoccupied = xtt[-occupied, ]
    
    if (avail == 1)
      POS = rbind(POS, unoccupied)
    else if (avail == 2)
    {
      if (r < 0.5)
        POS = rbind(POS, unoccupied[1,])
      else
        POS = rbind(POS, unoccupied[2,])
    }
    else if (avail == 3)
    {
      if (r < 1/3)
        POS = rbind(POS, unoccupied[1,])
      else if (r < 2/3)
        POS = rbind(POS, unoccupied[2,])
      else 
        POS = rbind(POS, unoccupied[3,])
    }
  }  
  plot(POS[,1], POS[,2], xlab = 'x', ylab = 'y', type = 'o', main = paste0("Stop at iter = ", iter))
  #return(POS)
}
#set.seed(1)
par(mfrow = c(2,2))
growPolymer(N)

