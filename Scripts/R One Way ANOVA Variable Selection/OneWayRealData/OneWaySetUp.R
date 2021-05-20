OneWaySetUp <- function(K)
{
  
  # Basics
  # Number of MCMC Samples
  nSample = 100000
  
  # Cut off Post Prob - ignores combos below this cum Prob
  cumPostProb = 0.99
  
  ## Create Matrix Representation of all possible combinations ##
  MatFinal <- genMatFinal(K)
  
  # Priors
  tau = 1000
  sh = 1
  sc = 1
  prob = matrix(1, nrow = length(MatFinal), ncol = 1) * (1/(length(MatFinal)))
  
  priors = list(tau,sh,sc,prob)
  outList <- list(nSample,cumPostProb,MatFinal,priors)
  return(outList)
  
}