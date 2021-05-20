OneWayDataFitModel <- function(tdData,K)
{
  y = tdData[[1]]
  xFactor = tdData[[2]]

  # Number of Observations
  n = length(y)

  mleEst = mleOneWay(y,xFactor,K,n)
  mleMu = mleEst[[1]]
  mleSig = mleEst[[2]]
  mleOver = mleEst[[3]]

  nLevel = matrix(0, nrow = K, ncol = 1)

  outList = list(n,nLevel,K,y,xFactor,mleMu,mleSig,mleOver)
  return(outList)
  
}