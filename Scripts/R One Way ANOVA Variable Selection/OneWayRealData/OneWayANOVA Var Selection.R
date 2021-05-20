# One Way ANOVA - Variable Dimension
# Clear Work Space
rm(list=ls(all=TRUE))
graphics.off()
# Load package
# library(nlme)
# Load Functions
source("genData.R")
debug(genData)
source("loadedDie.R")
debug(loadedDie)
# Basics
# Number of Observations
trueParam.n = 5000
# Number of Levels
trueParam.K = 6
# Means for each level
trueParam.mu = rep(0,trueParam.K)
trueParam.mu[1] = 1
trueParam.mu[2] = 1.1
trueParam.mu[3] = 4
trueParam.mu[4] = 4
trueParam.mu[5] = -2;
trueParam.mu[6] = 2;
# Variance (same for all levels)
trueParam.sig2 = 0.5
# Generate Synthetic Data (each level is euqlly likely)

prob = rep(0,trueParam.K)
for(k in 1:trueParam.K)
{
  prob[k] = 1/trueParam.K;  
}
## Generate Synthetic Data ##
#tdData <- genData(trueParam,prob)
dData.xBase = rep(0,trueParam.n)
dData.y = rep(0,trueParam.n)

for(i in 1:trueParam.n)
{
  u = runif(1)
  cumProb = 0
  for(j in 1:length(prob))
  {
    if(u > cumProb)
    {
      cumProb = cumProb + prob[j]
      dieValue = j
    }
    else
    {
      break
    }
  }
  dData.xBase[i] = dieValue 
  dData.y[i] = trueParam.mu[dData.xBase[i]] + sqrt(trueParam.sig2)*rnorm(1)
}

## Create Matrix Representation of all possible combinations ##
m1 = matrix(c(1, 0, 0, 1),nrow=2,ncol=2)
m2 = matrix(c(1, 1),nrow=2,ncol=1)
MatBase = list(m1,m2)
if(trueParam.K > 2)
{
  for(i in 2:(trueParam.K-1))
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
    
    if(i < (trueParam.K-1))
      MatBase = MatFinal
    
  }
}
if(trueParam.K == 2)
{
  MatFinal = MatBase
}

## Create groups combo map (May not be strictly neccesary, but follows notes) ##
groupComboMat = matrix(0,nrow = length(MatFinal), ncol = (2^(trueParam.K-1)-1))

twos = matrix(1)
for(i in 2:(trueParam.K))
  twos = rbind(2^(i-1),twos)

for(i in 1:length(MatFinal))
{
  for(j in 1:dim(MatFinal[[i]])[2])
  {
    s = sum(MatFinal[[i]][,j]*twos)
    if(s < (2^(trueParam.K-1)))
      groupComboMat[i,s] = 1;
  }
}

## Create Dummy Variables for 2,3,..,K and all combos (excluding 1) Have the numbering match the numbering from groupComboMat ##
dData.xDummyFull = matrix(0,nrow = trueParam.n, ncol = (2^(trueParam.K-1)-1))
dData.xDummyBase = matrix(0,nrow = trueParam.n, ncol = trueParam.K-1)

backTwos = matrix(0,nrow = (2^(trueParam.K-1)-1), ncol = trueParam.K-1)
backTwos[1,trueParam.K-1] = 1
for(i in 2:(2^(trueParam.K-1)-1))
{
  breakJ = 1
  for(j in 1:(trueParam.K-1))
  {
    rj = (trueParam.K-1-j+1)
    if(rj < breakJ)
      backTwos[i,rj] = backTwos[i-1,rj]
    else
    {
      if(backTwos[i-1,rj] == 0)
      {
        backTwos[i,rj] = 1;
        breakJ = rj
      }
    }
  }
    
}

for (i in 1:trueParam.n)
{
  if(dData.xBase[i] > 1)
    dData.xDummyBase[i,(dData.xBase[i]-1)] = 1
}

for (i in 1:(2^(trueParam.K-1)-1))
{
  for(j in 1:(trueParam.K-1))
  {
    if(backTwos[i,j] == 1)
    {
      dData.xDummyFull[,i] = dData.xDummyFull[,i] + dData.xDummyBase[,j]
    }
  }
}

## Generate the Posterior Probabilities for each group ##
cumPostProb = 0.99
p = length(MatFinal) + 1
n = trueParam.n
tau = 1000
sh = 1
sc = 1
y = dData.y
twoPi = 2*pi
lnTwoPi = log(twoPi)

postProb = matrix(0, nrow = p, ncol = 1)
lnPostProb = postProb
priorProb = matrix(1, nrow = p, ncol = 1) * (1/(p))
lnPriorProb = log(priorProb)

for (i in 1:length(MatFinal))
{

  groupIndex = which(groupComboMat[i,] == 1)
  pg = length(groupIndex) + 1
  X = dData.xDummyFull[,groupIndex]
  X = cbind(matrix(1,nrow = n,ncol =1),X)
  B = (t(X) %*% X + (1/tau)*diag(pg))
  b = (t(X) %*% y)
  betaMean = solve(B,b)
  btB1b = t(b) %*% betaMean
  yty_m_btB1b = t(y) %*% y - btB1b
  
  lnPostProb[i] = -0.5*log(det(B)) - 0.5*(n + p - pg)*lnTwoPi - lnPriorProb[i] + lgamma(0.5*(n-pg) + sh) - (0.5*(n-pg) + sh)*log(sc + 0.5*yty_m_btB1b)
   
}

pg = 0
B = ((1/tau)*diag(pg))
yty_m_btB1b = t(y) %*% y

lnPostProb[p] = -0.5*log(det(B)) - 0.5*(n + p - pg)*lnTwoPi - lnPriorProb[i] + lgamma(0.5*(n-pg) + sh) - (0.5*(n-pg) + sh)*log(sc + 0.5*yty_m_btB1b)

postProb = lnPostProb - max(lnPostProb)
postProb = exp(postProb)
postProb = postProb/sum(postProb)

orderedPostProbIndex = order(postProb)
topPostProbOrderedIndex = which(cumsum(postProb[orderedPostProbIndex]) > (1-cumPostProb))
topPostProbIndex = orderedPostProbIndex[topPostProbOrderedIndex] 

## Calculate Posterior Prob each Level Equal to 0 ##
postProbLevelZero = matrix(0,nrow = trueParam.K, ncol = 1)
tLnPostProbLevelNotZero = matrix(0,nrow = length(MatFinal), ncol = 1)
tLnPostProbLevelZero = matrix(0,nrow = length(MatFinal), ncol = 1)
## Don't think this is worth doing, but to do it run with intercept = 1 and intercept = 0
## then rotate through the base case.  Does not tell you much really.

## Calculate MCMC based summaries ##
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

SSBAllM = matrix(0,nrow = trueParam.K, ncol = 1)
SSBAllV = matrix(0,nrow = trueParam.K, ncol = 1)
SSSigAllM = 0
SSSigAllV = 0
nAll = 0
nSample = 1000
probDiffMat = matrix(1,nrow = trueParam.K, ncol = trueParam.K)
probSameMat = matrix(0,nrow = trueParam.K, ncol = trueParam.K)

for (j in 1:length(topPostProbIndex))
{
  i = topPostProbIndex[j]
  groupIndex = which(groupComboMat[i,] == 1)
  pg = length(groupIndex) + 1
  X = dData.xDummyFull[,groupIndex]
  X = cbind(matrix(1,nrow = n,ncol =1),X)
  B = (t(X) %*% X + (1/tau)*diag(pg))
  b = (t(X) %*% y)
  betaMean = solve(B,b)
  cholBInv = chol(solve(B))
  btB1b = t(b) %*% betaMean
  yty_m_btB1b = t(y) %*% y - btB1b
 
  for(t in 1:nSample)
  {
    tsig2 = 1/rgamma(1,sh+0.5*(n-pg),sc+0.5*yty_m_btB1b)
    tsig = sqrt(tsig2)
    tb = betaMean + tsig*cholBInv %*% t(t(rnorm(pg)))
    ttb = matrix(0,nrow = pg, ncol = 1)
    ttb[1] = tb[1]
    for(l in 2:pg)
    {
      tb[l] = tb[l] + tb[1]
      ttb[pg-l+2] = tb[l]
    }
    SSBM[[j]] = SSBM[[j]] + ttb
    SSBV[[j]] = SSBV[[j]] + ttb*ttb
    SSSigM[j] = SSSigM[j] + tsig
    SSSigV[j] = SSSigV[j] + tsig*tsig
    
    if( runif(1) < postProb[i])
    {
      for(l in 1:trueParam.K)
      {
        k = which(MatFinal[[i]][l,] == 1)
        SSBAllM[l] = SSBAllM[l] + ttb[k]
        SSBAllV[l] = SSBAllV[l] + ttb[k]*ttb[k]
      }
      SSSigAllM = SSSigAllM + tsig
      SSSigAllV = SSSigAllV + tsig*tsig
      nAll = nAll + 1
    }
  }
  
  SSBM[[j]] = SSBM[[j]]*(1/nSample)
  SSBV[[j]] = sqrt(SSBV[[j]]*(1/nSample) - SSBM[[j]]*SSBM[[j]])
  SSSigM[j] = SSSigM[j]*(1/nSample)
  SSSigV[j] = sqrt(SSSigV[j]*(1/nSample) - SSSigM[j]*SSSigM[j])

  for(l in 1:pg)
  {
    for(m1 in 1:(trueParam.K-1))
    {
      for(m2 in (m1+1):trueParam.K)
      {
        if((MatFinal[[i]][m1,l] == 1) && (MatFinal[[i]][m2,l] == 1))
        {
             probDiffMat[m1,m2] = probDiffMat[m1,m2] - postProb[i]
             probDiffMat[m2,m1] = probDiffMat[m2,m1] - postProb[i]
             probSameMat[m1,m2] = probSameMat[m1,m2] + postProb[i]
             probSameMat[m2,m1] = probSameMat[m2,m1] + postProb[i]        }           
      }
    }
  }
  
}

SSBAllM = SSBAllM*(1/nAll)
SSBAllV = sqrt(SSBAllV*(1/nAll) - SSBAllM*SSBAllM)
SSSigAllM = SSSigAllM*(1/nAll)
SSSigAllV = sqrt(SSSigAllV*(1/nAll) - SSSigAllM*SSSigAllM)

## Summaries ##
barplot(t(postProb),names.arg=topPostProbIndex)
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
  print(paste("Posterior Sig Mean:",SSSigM[i]))
  print(paste("Posterior Sig Std:",SSSigV[i]))  
}

print(paste("------------------------------------"))
print(paste("------------------------------------"))
print(paste("Prob all means equal to zero:", postProb[length(MatFinal)+1]))
print(paste("Prob at least one mean different from rest:", 1-(postProb[length(MatFinal)]+postProb[length(MatFinal)+1])))

for(i in 1:trueParam.K)
{
  print(paste("------------------------------------"))
  print(paste("Mu:",i))
  print(paste("True:",trueParam.mu[i]))
  print(paste("Posterior Mean:",SSBAllM[i]))
  print(paste("Posterior Std:",SSBAllV[i]))  
}

print(paste("True Sig:",sqrt(trueParam.sig2)))
print(paste("Posterior Sig Mean:",SSSigAllM))
print(paste("Posterior Sig Std:",SSSigAllV))
