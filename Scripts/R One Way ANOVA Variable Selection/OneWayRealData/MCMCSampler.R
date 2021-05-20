MCMCSampler <- function(SSParams,postProbSummary,groupComboMat,MatFinal,nSample,n,nLevel,K,y,xDummyFull,priors)
{
  
  pTau = priors[[1]]
  pSh = priors[[2]]
  pSc = priors[[3]]

  twos = matrix(1)
  for(i in 2:(K))
    twos = rbind(2^(i-1),twos)
  
  postProb = postProbSummary[[1]]
  rderedPostProbIndex = postProbSummary[[2]]
  topPostProbIndex = postProbSummary[[3]]
  
  SSBM = SSParams[[1]]
  SSBV = SSParams[[2]]
  SSSigM = SSParams[[3]]
  SSSigV = SSParams[[4]]
  SSBallM = SSParams[[5]]
  SSBallV = SSParams[[6]]
  SSSigAllM = SSParams[[7]]
  SSSigAllV = SSParams[[8]]
  SSBOverM = SSParams[[9]]
  SSBOverV = SSParams[[10]]
  
  nAll = 0
  
  probDiffMat = matrix(1,nrow = K, ncol = K)
  probSameMat = matrix(0,nrow = K, ncol = K)
  
  for (j in 1:length(topPostProbIndex))
  {
    
    i = topPostProbIndex[j]
    
    groupIndex = which(groupComboMat[i,] == 1)
    groupOrderedIndex = groupIndex
    pg = dim(MatFinal[[i]])[2]
    
    for(l in 1:pg)
    {
      s = sum(MatFinal[[i]][,l]*twos)
      if(s < (2^(K-1)))
      {
        if(l == 1)
          X = xDummyFull[,s]
        else
          X = cbind(X,xDummyFull[,s])
      }
      else
      {
        interceptIndex = l
        if(l == 1)
          X = matrix(1,nrow = n,ncol =1)
        else
          X = cbind(X,matrix(1,nrow = n,ncol =1))
      }
      groupOrderedIndex[l] = s
    }

    B = (t(X) %*% X + (1/pTau)*diag(pg))
    b = (t(X) %*% y)
    betaMean = solve(B,b)
    tCholBInv = t(chol(solve(B)))
    btB1b = t(b) %*% betaMean
    yty_m_btB1b = t(y) %*% y - btB1b
    
    for(t in 1:nSample)
    {
      tsig2 = 1/rgamma(1,pSh+0.5*(n-pg),pSc+0.5*yty_m_btB1b)
      tsig = sqrt(tsig2)
      tb = betaMean + tsig*tCholBInv %*% t(t(rnorm(pg)))
      for(l in 1:pg)
      {
        if(l == interceptIndex)
        {
           SSBM[[j]][l] = SSBM[[j]][l] + tb[interceptIndex]
           SSBV[[j]][l] = SSBV[[j]][l] + tb[interceptIndex]*tb[interceptIndex]
        }
        else
        {
          SSBM[[j]][l] = SSBM[[j]][l] + tb[interceptIndex] + tb[l]
          SSBV[[j]][l] = SSBV[[j]][l] + (tb[interceptIndex] + tb[l])*(tb[interceptIndex] + tb[l])
        }
      }

      SSSigM[j] = SSSigM[j] + tsig
      SSSigV[j] = SSSigV[j] + tsig*tsig

      if( runif(1) < postProb[i])
      {
        tBOver = 0
        for(l in 1:pg)
        {
          k = which(MatFinal[[i]][,l] == 1)
          for(m in 1:length(k))
          {
            if(MatFinal[[i]][1,l] == 1)
            {
              SSBallM[k[m]] = SSBallM[k[m]] + tb[interceptIndex]
              SSBallV[k[m]] = SSBallV[k[m]] + tb[interceptIndex]*tb[interceptIndex]
              tBOver = tBOver + tb[interceptIndex]*nLevel[l]
            }
            else
            {
              SSBallM[k[m]] = SSBallM[k[m]] + tb[interceptIndex] + tb[l]
              SSBallV[k[m]] = SSBallV[k[m]] + (tb[interceptIndex] + tb[l])*(tb[interceptIndex] + tb[l])
              tBOver = tBOver + (tb[interceptIndex] + tb[l])*nLevel[l]
            }
          }
        }
        SSSigAllM = SSSigAllM + tsig
        SSSigAllV = SSSigAllV + tsig*tsig
        tBOver = tBOver/n
        SSBOverM = SSBOverM + tBOver
        SSBOverV = SSBOverV + tBOver*tBOver
        nAll = nAll + 1
      }
    }
    
    SSBM[[j]] = SSBM[[j]]*(1/nSample)
    SSBV[[j]] = sqrt(SSBV[[j]]*(1/nSample) - SSBM[[j]]*SSBM[[j]])
    SSSigM[j] = SSSigM[j]*(1/nSample)
    SSSigV[j] = sqrt(SSSigV[j]*(1/nSample) - SSSigM[j]*SSSigM[j])
    
    for(l in 1:pg)
    {
      for(m1 in 1:(K-1))
      {
        for(m2 in (m1+1):K)
        {
          if((MatFinal[[i]][m1,l] == 1) && (MatFinal[[i]][m2,l] == 1))
          {
            probDiffMat[m1,m2] = probDiffMat[m1,m2] - postProb[i]
            probDiffMat[m2,m1] = probDiffMat[m2,m1] - postProb[i]
            probSameMat[m1,m2] = probSameMat[m1,m2] + postProb[i]
            probSameMat[m2,m1] = probSameMat[m2,m1] + postProb[i]        
          }           
        }
      }
    }
    
  }
  
  SSBallM = SSBallM*(1/nAll)
  SSBallV = sqrt(SSBallV*(1/nAll) - SSBallM*SSBallM)
  SSSigAllM = SSSigAllM*(1/nAll)
  SSSigAllV = sqrt(SSSigAllV*(1/nAll) - SSSigAllM*SSSigAllM)
  SSBOverM = SSBOverM*(1/nAll)
  SSBOverV = sqrt(SSBOverV*(1/nAll) - SSBOverM*SSBOverM)
  
  MCMCOut <- list(SSBM,SSBV,SSSigM,SSSigV,SSBallM,SSBallV,SSSigAllM,SSSigAllV,SSBOverM,SSBOverV,probDiffMat,probSameMat)
  
  return(MCMCOut)
  
}