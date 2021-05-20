mleOneWay <- function(y,xFactor,K,n)
{
  
  mleMu = matrix(0,nrow = K, ncol = 1)
  Denom = matrix(0,nrow = K, ncol = 1)  
  mleSig = 0
  
  for(t in 1:n)
  {
    mleMu[xFactor[t]] = mleMu[xFactor[t]] + y[t]
    Denom[xFactor[t]] = Denom[xFactor[t]] + 1
  }
  
  for(l in 1:K)
    mleMu[l] = mleMu[l]/Denom[l]
  
  for(t in 1:n)
  {
    mleSig = mleSig + (y[t] - mleMu[xFactor[t]])^2
  }
  
  mleSig = sqrt(mleSig/n)
  mleOver = mean(y)
 
  mleEstOut = list(mleMu,mleSig,mleOver)
  
  return(mleEstOut)
  
}