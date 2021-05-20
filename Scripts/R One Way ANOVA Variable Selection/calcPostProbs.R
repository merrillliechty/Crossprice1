calcPostProbs <- function(MatFinal,n,priors,y,xDummyFull,groupComboMat,cumPostProb)
{
  pTau = priors[[1]]
  pSh = priors[[2]]
  pSc = priors[[3]]
  pProb = priors[[4]]

  p = length(MatFinal)
  # p = length(MatFinal)
  
  twoPi = 2*pi
  lnTwoPi = log(twoPi)
  
  postProb = matrix(0, nrow = p, ncol = 1)
  lnPostProb = postProb
  lnPriorProb = log(pProb)
  
  for (i in 1:p)
  {
    
    groupIndex = which(groupComboMat[i,] == 1)
    pg = length(groupIndex) + 1
    X = xDummyFull[,groupIndex]
    X = cbind(matrix(1,nrow = n,ncol =1),X)
    B = (t(X) %*% X + (1/pTau)*diag(pg))
    b = (t(X) %*% y)
    betaMean = solve(B,b)
    btB1b = t(b) %*% betaMean
    yty_m_btB1b = t(y) %*% y - btB1b
    
    lnPostProb[i] = -0.5*log(det(B)) - 0.5*(n + p - pg)*lnTwoPi - lnPriorProb[i] + lgamma(0.5*(n-pg) + pSh) - (0.5*(n-pg) + pSh)*log(pSc + 0.5*yty_m_btB1b)
    
  }
  
#  pg = 0
#  B = ((1/pTau)*diag(pg))
#  yty_m_btB1b = t(y) %*% y
  
#  lnPostProb[p] = -0.5*log(det(B)) - 0.5*(n + p - pg)*lnTwoPi - lnPriorProb[i] + lgamma(0.5*(n-pg) + pSh) - (0.5*(n-pg) + pSh)*log(pSc + 0.5*yty_m_btB1b)
  
  postProb = lnPostProb - max(lnPostProb)
  postProb = exp(postProb)
  postProb = postProb/sum(postProb)
  
  orderedPostProbIndex = order(postProb)
  topPostProbOrderedIndex = which(cumsum(postProb[orderedPostProbIndex]) > (1-cumPostProb))
  topPostProbIndex = orderedPostProbIndex[topPostProbOrderedIndex] 
  
  postProbSum = list(postProb,orderedPostProbIndex,topPostProbIndex)
  
  return(postProbSum)
  
}