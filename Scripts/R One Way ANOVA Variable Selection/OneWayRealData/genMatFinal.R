genMatFinal <- function(K)
{
  
  m1 = matrix(c(1, 0, 0, 1),nrow=2,ncol=2)
  m2 = matrix(c(1, 1),nrow=2,ncol=1)
  MatBase = list(m1,m2)
  if(K > 2)
  {
    for(i in 2:(K-1))
    {
      nBaseMatrix = length(MatBase)
      nFinalMatrix = 0;
      MatFinal = list()
      for(j in 1:nBaseMatrix)
      {
        nColCurrBase = dim(MatBase[[nBaseMatrix - j + 1]])[2]
        nRowCurrBase = dim(MatBase[[nBaseMatrix - j + 1]])[1]
        nFinalMatrix = nFinalMatrix + nColCurrBase + 1;
        for(l in 1:(nColCurrBase+1))
        {
          tM = rbind(rep(0,nColCurrBase),MatBase[[nBaseMatrix - j + 1]])
          if(l == (nColCurrBase+1))
          {
            tM = cbind(t(t(rep(0,nRowCurrBase+1))),tM)
            tM[1,1] = 1
          }
          else
          {
            tM[1,(nColCurrBase - l + 1)] = 1
          }
          if(length(MatFinal) == 0)
            MatFinal = tM
          else
            if((j == 1) && (l == 2))
              MatFinal = append(list(tM),list(MatFinal))
            else
              MatFinal = append(list(tM),MatFinal)
        }
      }
      
      if(i < (K-1))
        MatBase = MatFinal
      
    }
  }
  if(K == 2)
  {
    MatFinal = MatBase
  }
  
  return(MatFinal)
  
}