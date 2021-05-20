createXDummyFull <- function(n,K,xFactor)
{
  
  xDummyFull = matrix(0,nrow = n, ncol = (2^(K-1)-1))
  xDummyBase = matrix(0,nrow = n, ncol = K-1)
  
  backTwos = matrix(0,nrow = (2^(K-1)-1), ncol = K-1)
  backTwos[1,K-1] = 1
  for(i in 2:(2^(K-1)-1))
  {
    breakJ = 1
    for(j in 1:(K-1))
    {
      rj = (K-1-j+1)
      if(rj < breakJ)
        backTwos[i,rj] = backTwos[i-1,rj]
      else
      {
        if(backTwos[i-1,rj] == 0)
        {
          backTwos[i,rj] = 1;
          breakJ = rj
        }
      }
    }
    
  }
  
  for (i in 1:n)
  {
    if(xFactor[i] > 1)
      xDummyBase[i,(xFactor[i]-1)] = 1
  }
  
  for (i in 1:(2^(K-1)-1))
  {
    for(j in 1:(K-1))
    {
      if(backTwos[i,j] == 1)
      {
        xDummyFull[,i] = xDummyFull[,i] + xDummyBase[,j]
      }
    }
  }
  
  return(xDummyFull)
  
}