OneWayData <- function(K)
{
  # browser()
  
  # Number of Observations
  n = 5000
  # Number of Levels
  if(K != 6)
  {
    print(paste("------------------------------------"))
    print(paste("Error Expecting Number of Level to be 6 - for data generation purposes"))
    break
  }
  # Means for each level
  tmu = rep(0,K)
  tmu[1] = 1
  tmu[2] = 1.1
  tmu[3] = 4
  tmu[4] = 4
  tmu[5] = -2;
  tmu[6] = 2;
  # Variance (same for all levels)
  tsig2 = 0.5

  # Generate Synthetic Data (each level is equally likely)
  tdData <- genData(tmu,tsig2,K,n)
  y = tdData[[1]]
  xFactor = tdData[[2]]
  
  mleEst = mleOneWay(y,xFactor,K,n)
  mleMu = mleEst[[1]]
  mleSig = mleEst[[2]]
  mleOver = mleEst[[3]]
  mleSSR = mleEst[[4]]

  nLevel = matrix(0, nrow = K, ncol = 1)
  tMuOver = 0
  
  for(i in 1:K)
  {
    nLevel[i] = sum(xFactor == i)
    tMuOver = tMuOver + tmu[i]*nLevel[i]
  }
  
  tMuOver = tMuOver/n
  
  outList = list(n,nLevel,K,y,xFactor,mleMu,mleSig,mleOver,mleSSR,tmu,tsig2,tMuOver)
  
  return(outList)
  
}