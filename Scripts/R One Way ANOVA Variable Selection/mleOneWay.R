mleOneWay <- function(y,xFactor,K,n)
{
  
  mleMu = matrix(0,nrow = K, ncol = 1)
  Denom = matrix(0,nrow = K, ncol = 1)  
  mleOver = mean(y)
  mleSig = 0
  mleSSR = 0
  
  for(t in 1:n)
  {
    mleMu[xFactor[t]] = mleMu[xFactor[t]] + y[t]
    Denom[xFactor[t]] = Denom[xFactor[t]] + 1
  }
  
  for(l in 1:K)
  {
    mleMu[l] = mleMu[l]/Denom[l]
    mleSSR = mleSSR + Denom[l]*((mleMu[l] - mleOver)^2)
  }
  
  for(t in 1:n)
  {
    mleSig = mleSig + (y[t] - mleMu[xFactor[t]])^2
  }
  
  mleSig = sqrt(mleSig/n)
 
  mleEstOut = list(mleMu,mleSig,mleOver,mleSSR)
  
  return(mleEstOut)
  
}