genData <- function(tmu,tsig2,K,n)
{
  
  prob = rep(0,K)
  for(k in 1:K)
  {
    prob[k] = 1/K;  
  }
  
  xFactor = rep(0,n)
  y = rep(0,n)
  
  for (i in 1:n)
  {
    xFactor[i] = loadedDie(prob) 
    y[i] = tmu[xFactor[i]] + sqrt(tsig2)*rnorm(1)
  }
  
  return(list(y,xFactor))

}