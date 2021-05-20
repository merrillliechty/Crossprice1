initSSParams <- function(MatFinal,topPostProbIndex,K)
{

  SSBM = list()
  SSBV = list()
  SSSigM = matrix(0,nrow = length(topPostProbIndex),ncol = 1)
  SSSigV = matrix(0,nrow = length(topPostProbIndex),ncol = 1)
  
  for (i in 1:length(topPostProbIndex))
  {
    tM = matrix(0,nrow = dim(MatFinal[[topPostProbIndex[i]]])[2], ncol = 1) 
    if(i == 1)
    {
      SSBM = list(tM)
      SSBV = list(tM)
    }
    else
    {
      SSBM = append(SSBM,list(tM))
      SSBV = append(SSBV,list(tM))
    }
  }
  
  SSBallM = matrix(0,nrow = K, ncol = 1)
  SSBallV = matrix(0,nrow = K, ncol = 1)
  SSSigAllM = 0
  SSSigAllV = 0
  SSBOverM = 0
  SSBOverV = 0
  SSSSEM = 0
  SSSSEV = 0
  SSMSESSEM = 0
  SSMSESSEV = 0
  
  SSParams = list(SSBM,SSBV,SSSigM,SSSigV,SSBallM,SSBallV,SSSigAllM,SSSigAllV,SSBOverM,SSBOverV,SSSSEM,SSSSEV,SSMSESSEM,SSMSESSEV)
  
  return(SSParams)
  
}