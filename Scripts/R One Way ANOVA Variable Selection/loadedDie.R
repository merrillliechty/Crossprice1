loadedDie <- function(prob)
{
  
  u = runif(1)
  cumProb = 0
  for(j in 1:length(prob))
  {
    if(u > cumProb)
    {
      cumProb = cumProb + prob[j]
      dieValue = j
    }
    else
    {
      break
    }
  }
 
  return(dieValue)
  
}