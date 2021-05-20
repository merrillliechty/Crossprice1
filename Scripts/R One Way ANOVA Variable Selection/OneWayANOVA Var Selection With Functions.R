# One Way ANOVA - Variable Dimension
# Clear Work Space
rm(list=ls(all=TRUE))
graphics.off()
# Load package
# Load Functions
source("genData.R")
#debug(genData)
source("loadedDie.R")
#debug(loadedDie)
source("mleOneWay.R")
#debug(mleOneWay)
source("genMatFinal.R")
#debug(genMatFinal)
source("genGroupComboMat.R")
#debug(genGroupComboMat)
source("createXDummyFull.R")
#debug(createXDummyFull)
source("calcPostProbs.R")
#debug(calcPostProbs)
source("initSSParams.R")
#debug(initSSParams)
source("MCMCSampler.R")
#debug(MCMCSampler)

# Basics
# Number of MCMC Samples
nSample = 1000
# Number of Observations
n = 500
# Number of Levels
K = 6
# Means for each level
tmu = rep(0,K)
tmu[1] = 1
tmu[2] = 1.1
tmu[3] = 4
tmu[4] = 4
tmu[5] = -2;
tmu[6] = 2;
# Variance (same for all levels)
tsig2 = 0.5
# Cut off Post Prob - ignores combos below this cum Prob
cumPostProb = 0.99
# Priors
tau = 1000
sh = 1
sc = 1

# Generate Synthetic Data (each level is equally likely)
tdData <- genData(tmu,tsig2,K,n)
y = tdData[[1]]
xFactor = tdData[[2]]
rm(tdData)

mleEst = mleOneWay(y,xFactor,K,n)
mleMu = mleEst[[1]]
mleSig = mleEst[[2]]
mleOver = mleEst[[3]]

nLevel = matrix(0, nrow = K, ncol = 1)
tMuOver = 0

for(i in 1:K) 
{
  nLevel[i] = sum(xFactor == i)
  tMuOver = tMuOver + tmu[i]*nLevel[i]
}

tMuOver = tMuOver/n


## Create Matrix Representation of all possible combinations ##
MatFinal <- genMatFinal(K)

## Create groups combo map (May not be strictly neccesary, but follows notes) ##
groupComboMat <- genGroupComboMat(K,MatFinal)

## Create Dummy Variables for 2,3,..,K and all combos (excluding 1) Have the numbering match the numbering from groupComboMat ##
xDummyFull <- createXDummyFull(n,K,xFactor)

## Generate the Posterior Probabilities for each group ##
prob = matrix(1, nrow = length(MatFinal), ncol = 1) * (1/(length(MatFinal)))
priors = list(tau,sh,sc,prob)
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
SSBOverM = MCMCOut[[9]]
SSBOverV = MCMCOut[[10]]
SSANOVASigM = MCMCOut[[11]]
SSANOVASigV = MCMCOut[[12]]
probDiffMat = MCMCOut[[13]]
probSameMat = MCMCOut[[14]]

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

for(i in 1:K)
{
  print(paste("------------------------------------"))
  print(paste("Mu:",i))
  print(paste("True:",tmu[i]))
  print(paste("MLE:",mleMu[i]))
  print(paste("Posterior Mean:",SSBallM[i]))
  print(paste("Posterior Std:",SSBallV[i]))  
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
print(paste("Posterior ANOVA Sig Mean:",SSANOVASigM))
print(paste("Posterior ANOVA Sig Std:",SSANOVASigV))
