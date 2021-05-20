OneWayAnalysis <- function(startList,dataList)
{
  nSample = startList[[1]]
  cumPostProb = startList[[2]]
  MatFinal = startList[[3]]
  priors = startList[[4]]

  n = dataList[[1]]
  nLevel = dataList[[2]]
  K = dataList[[3]]
  y = dataList[[4]]
  xFactor = dataList[[5]]
  
  ## Create groups combo map (May not be strictly neccesary, but follows notes) ##
  groupComboMat <- genGroupComboMat(K,MatFinal)
  ####
  ## Create Dummy Variables for 2,3,..,K and all combos (excluding 1) Have the numbering match the numbering from groupComboMat ##
  xDummyFull <- createXDummyFull(n,K,xFactor)

  ## Generate the Posterior Probabilities for each group ##
  postProbSummary <- calcPostProbs(MatFinal,n,priors,y,xDummyFull,groupComboMat,cumPostProb)
  postProb = postProbSummary[[1]]
  orderedPostProbIndex = postProbSummary[[2]]
  topPostProbIndex = postProbSummary[[3]]
  
  ## Set up summary structures ##
  SSParams <- initSSParams(MatFinal,topPostProbIndex,K)
  
  ## Calculate MCMC based summaries ##
  MCMCOut <- MCMCSampler(SSParams,postProbSummary,groupComboMat,MatFinal,nSample,n,nLevel,K,y,xDummyFull,priors)
  
  SSBM = MCMCOut[[1]]
  SSBV = MCMCOut[[2]]
  SSSigM = MCMCOut[[3]]
  SSSigV = MCMCOut[[4]]
  SSBallM = MCMCOut[[5]]
  SSBallV = MCMCOut[[6]]
  SSSigAllM = MCMCOut[[7]]
  SSSigAllV = MCMCOut[[8]]
  SSBGrandM = MCMCOut[[9]]
  SSBGrandV = MCMCOut[[10]]
  SSSSRM = MCMCOut[[11]]
  SSSSRV = MCMCOut[[12]]
  probDiffMat = MCMCOut[[13]]
  probSameMat = MCMCOut[[14]]
  twoTailSingleLevelBayesPvalue = MCMCOut[[15]]
  factorBayesPValue = MCMCOut[[16]]
  
  outList = list(postProb,orderedPostProbIndex,topPostProbIndex,SSBM,SSBV,SSSigM,SSSigV,SSSigAllM,SSSigAllV,SSBallM,SSBallV,SSBGrandM,SSBGrandV,SSSSRM,SSSSRV,twoTailSingleLevelBayesPvalue,factorBayesPValue,probDiffMat,probSameMat)
  
  return(outList)
    
}