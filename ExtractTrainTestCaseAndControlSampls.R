#Author : Tina Khajeh
#this is a preprocessing code in order to seperate training, test, case, and control samples for futher procesing in our
#proposed learning algorithm and for DEG verifivations.
 
 
setwd("/Users/Tina/Desktop/Thesis Code")
library(limma)
library(GEOquery)
library(seqinr)
library(Biostrings)
library(rat2302probe)
#Mapping : 
scottData = read.table("Treated to control mapping for SEQC samples-with Run_s.txt",header=TRUE,fill = TRUE,sep = "\t")


#-----------------------------------------------------------------------------|
#Trainig, Test, Case and Control data:                                        |
#-----------------------------------------------------------------------------|
trainingDataIndex = as.character(scottData$Original.SET.Annotation)=="TRAINING SET"#62 
testDataIndex = as.character(scottData$Original.SET.Annotation)=="TEST SET"#42 
trainingData = scottData[trainingDataIndex,]
testData = scottData[testDataIndex,]
#Case :
trainCaseIndex = trainingData$Control.Group.Mapping!="X"
trainCase = trainingData[trainCaseIndex,]

testCaseIndex = testData$Control.Group.Mapping != "X"
testCase = testData[testCaseIndex,]

#Control :
trainControlIndex = trainingData$Control.Group!="X"
trainControl = trainingData[trainControlIndex,]

testControlIndex = testData$Control.Group != "X"
testControl = testData[testControlIndex,]

controlData = rbind(trainControl, testControl)
caseData = rbind(trainCase, testCase)
samples = getGEO(GEO = "GSE47875", destdir = "/Users/Tina/Desktop/Thesis Code",GSEMatrix = FALSE)
sampleList = GSMList(samples)
namesList = c()
longNames = c()
for( i in 1:length(sampleList))
{
  metaSampleList =Meta(sampleList[[i]])
  longNames = c(longNames,metaSampleList$title)
  name = names(sampleList)[i]
  namesList = c(namesList,name)
}

namesAndGSMMicroarray = data.frame(longNames,namesList)

chemical = "NIT"#change it manually for each chemical
indexes = which(as.character(caseData$Identifier.for.Chemical.model) == chemical)
group = unique(caseData$Control.Group.Mapping[indexes])
if(length(group)!=1){
  print("Oh!")
}
tmpmicCase = as.character(caseData$Microarray.ID[indexes])
micIndex = which(namesAndGSMMicroarray$longNames %in% tmpmicCase)
micCase =  as.character(namesAndGSMMicroarray$namesList[micIndex])
RNACase = as.character(caseData$RNA.seq.ID[indexes])

#Control:
controlIndex = which(controlData$Control.Group==group)
tmpmicControl = as.character(controlData$Microarray.ID[controlIndex])
ind = which(namesAndGSMMicroarray$longNames %in% tmpmicControl)
micControl = as.character(namesAndGSMMicroarray$namesList[ind])
RNAControl = as.character(controlData$RNA.seq.ID[controlIndex])

load("3DataScaledAndRemoveMean.RData")
microarrayData = newRMAScaleremoveMean
RNAData = newRNAFileScale

micTestIndex = c(micControl , micCase)
RNATestIndex = c(RNAControl ,RNACase)


micTestControl = microarrayData[micControl,]
micTestCase = microarrayData[micCase,]
RNATestControl = RNAData[RNAControl,]
RNATestCase = RNAData[RNACase,]

indexTodeleteMic = rownames(microarrayData) %in% micTestIndex
indexTodeleteRNA = rownames(RNAData) %in% RNATestIndex

micTrain = microarrayData[!indexTodeleteMic,]
RNATrain = RNAData[!indexTodeleteMic,]

MicData = rbind(micTrain , micTestControl,micTestCase)
RNAData = rbind(RNATrain , RNATestControl, RNATestCase)

shuffledMicData<-MicData
shuffledRNAData<-RNAData


MicCaseAndControls = rbind(micTestControl,micTestCase)
RNACaseAndControls = rbind(RNATestControl, RNATestCase)

#save data:
write.csv(shuffledMicData, file = "shuffledMicData1_22.csv",row.names=FALSE)
write.csv(shuffledRNAData, file = "shuffledRNAData1_22.csv",row.names=FALSE)
write.csv(MicCaseAndControls, file = "MicCaseAndControls.csv",row.names=FALSE)
write.csv(RNACaseAndControls, file = "RNACaseAndControls.csv",row.names=FALSE)

