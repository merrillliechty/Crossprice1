genMatFinalReducedSpace <- function(K)
{

  m1 = matrix(c(1, 1),nrow=1,ncol=2)
  m2 = matrix(c(1, 0, 0, 1),nrow=2,ncol=2)
  MatFinal = list(m1,m2)
  if(K > 2)
  {
    MatFinal = list()
    for(i in 0:(2^(K-1)-1))
    {
    
      nRowCurrBase = K
      nColCurrBase = sum(as.numeric(intToBits(i))) + 1
      tM = matrix(0,nrow = nRowCurrBase, ncol = nColCurrBase)
      
        nColCount = 1
        tM[1,nColCount] = 1
        for(j in 1:K)
        {
          nColCount = nColCount + as.numeric(intToBits(i))[K-j+1]
          tM[j,nColCount] = 1
        }
      
      MatFinal = append(MatFinal,list((tM)))
    }
    
  }

    if(K == 2)
  {
    MatFinal = MatBase
  }
  
  return(MatFinal)
  
}

