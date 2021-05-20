OneWayDataGenData <- function(K)
{
  
  # Number of Observations
  n = 5000

  
  # Means for each level
  tmu = rep(0,K)
  for(i in 1:K)
  {
    tmu[i] = rnorm(1,K,.1)
  }
  
  
  # Variance (same for all levels)
  tsig2 = 0.5

  # Generate Synthetic Data (each level is equally likely)
  tdData <- genData(tmu,tsig2,K,n)
  y = tdData[[1]]
  xFactor = tdData[[2]]
  return(tdData)

}