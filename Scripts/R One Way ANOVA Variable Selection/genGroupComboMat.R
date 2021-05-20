genGroupComboMat <-function(K,MatFinal)
{
  groupComboMat = matrix(0,nrow = length(MatFinal), ncol = (2^(K-1)-1))
  
  twos = matrix(1)
  for(i in 2:(K))
    twos = rbind(2^(i-1),twos)
  
  for(i in 1:length(MatFinal))
  {
    for(j in 1:dim(MatFinal[[i]])[2])
    {
      s = sum(MatFinal[[i]][,j]*twos)
      if(s < (2^(K-1)))
        groupComboMat[i,s] = 1;
    }
  }
  
  return(groupComboMat)
  
}
  