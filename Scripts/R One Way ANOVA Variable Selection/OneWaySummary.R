OneWaySummary <- function(startList,dataList,analysisList)
{

  graphics.off()
 
  MatFinal = startList[[3]]
  
  n = dataList[[1]]
  mleMu = dataList[[6]]
  mleSig = dataList[[7]]
  mleOver = dataList[[8]]
  mleSSR = dataList[[9]]
  tmu = dataList[[10]]
  tsig2 = dataList[[11]]
  tMuOver = dataList[[12]]
  
  postProb = analysisList[[1]]
  orderedPostProbIndex = analysisList[[2]]
  topPostProbIndex = analysisList[[3]]
  SSBM = analysisList[[4]]
  SSBV = analysisList[[5]]
  SSSigM = analysisList[[6]]
  SSSigV = analysisList[[7]]
  SSSigAllM = analysisList[[8]]
  SSSigAllV = analysisList[[9]]
  SSBallM = analysisList[[10]]
  SSBallV = analysisList[[11]]
  SSBOverM = analysisList[[12]]
  SSBOverV = analysisList[[13]]
  SSSSRM = analysisList[[14]]
  SSSSRV = analysisList[[15]]
  twoTailSingleLevelBayesPvalue = analysisList[[16]]
  factorBayesPValue = analysisList[[17]]
  
  ## Summaries ##
  barplot(t(postProb),names.arg=1:length(MatFinal))
  barplot(t(postProb[topPostProbIndex]),names.arg=topPostProbIndex)
  
  for (i in 1:length(topPostProbIndex))
  {
    j = topPostProbIndex[i]
    print(paste("------------------------------------"))
    print(paste("Group:",j))
    print(paste("Top grouping: Post Prob:",postProb[j]))
    if (j == (length(MatFinal) + 1))
      j = length(MatFinal)
    for (l in 1:dim(MatFinal[[j]])[2])
    {
      print(paste("New level:",l))
      currentGroup = which(MatFinal[[j]][,l] == 1)
      print(paste("Contains original levels:",currentGroup))
      print(paste("Posterior Mean:", SSBM[[i]][l]))
      print(paste("Posterior Std:",SSBV[[i]][l]))
    }
    print(paste("------------------------------------"))
    print(paste("------------------------------------"))
    print(paste("Posterior Sig Mean:",SSSigM[i]))
    print(paste("Posterior Sig Std:",SSSigV[i]))  
  }
  
  print(paste("------------------------------------"))
  print(paste("------------------------------------"))
  print(paste("Prob at least one mean different from rest:", 1-(postProb[length(MatFinal)])))
  print(paste("Prob that all the means are the same:", (postProb[length(MatFinal)])))
  print(paste("Bayes P-value at least one mean different from overall mean:", factorBayesPValue))
  
  for(i in 1:K)
  {
    print(paste("------------------------------------"))
    print(paste("Mu:",i))
    print(paste("True:",tmu[i]))
    print(paste("MLE:",mleMu[i]))
    print(paste("Posterior Mean:",SSBallM[i]))
    print(paste("Posterior Std:",SSBallV[i]))  
    print(paste("Bayes two-tail P-value Ho Different from overall mean:",twoTailSingleLevelBayesPvalue[i]))
  }
  
  print(paste("------------------------------------"))
  print(paste("------------------------------------"))
  print(paste("True Sig:",sqrt(tsig2)))
  print(paste("MLE Sig:",mleSig))
  print(paste("Posterior Sig Mean:",SSSigAllM))
  print(paste("Posterior Sig Std:",SSSigAllV))
  
  print(paste("------------------------------------"))
  print(paste("------------------------------------"))
  print(paste("True Mu Overall:",tMuOver))
  print(paste("MLE Mu Overall:",mleOver))
  print(paste("Posterior Mu Overall Mean:",SSBOverM))
  print(paste("Posterior Mu Overall Std:",SSBOverV))
  
  print(paste("------------------------------------"))
  print(paste("------------------------------------"))
  print(paste("Posterior SSR ",SSSSRM))
  print(paste("Posterior SSR Std:",SSSSRV))
  print(paste("MLE SSR ",mleSSR))
  print(paste("MLE sig2*n:",n*mleSig*mleSig))
  
}