# One Way ANOVA - Variable Dimension
# Clear Work Space
rm(list=ls(all=TRUE))
# Load package
# Load Functions
source("loadOneWayFunctionsRealData.R")
loadOneWayFunctionsRealData()

K = 6
startList = OneWaySetUp(K)
tdData = OneWayDataGenData(K)
dataList = OneWayDataFitModel(tdData,K)

analysisList = OneWayAnalysis(startList,dataList)

OneWaySummaryRealData(startList,dataList,analysisList)

