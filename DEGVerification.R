#Author: Tina Khajeh
#Using Limma package the differentially expressed genes(DEGs) will be calculated.
#DEG is one of the criteria for checking the performance of proposed approach
 
setwd("/Users/Tina/Desktop/Thesis Code")
library(limma)
library(GEOquery)
library(seqinr)
library(Biostrings)
library(rat2302probe)

#"test.txt" is the normalized gene expression values predicted by the proposed method(the proposed method is in python)
predictedData = read.table(file = "test.txt",fill = TRUE,sep = "\t")
predictedData = as.matrix(predictedData)

MicCaseAndControls = read.table(file = "MicCaseAndControls.csv",fill = TRUE,sep = "\t")
MicCaseAndControls = as.matrix(MicCaseAndControls)

RNACaseAndControls = read.table(file = "RNACaseAndControls.csv",fill = TRUE,sep = "\t")
RNACaseAndControls = as.matrix(RNACaseAndControls)

#find genes with max variance:
diffArray = c()
for(i in 1:dim(predictedData)[2]){
  geneCol = predictedData[,i]
  tmpDiff = max(geneCol) - min(geneCol)
  diffArray = c(diffArray , tmpDiff)
}
sortedDiffArray = sort(diffArray, decreasing = TRUE)
top = sortedDiffArray
topIndexes = match(top , diffArray )

colnames(RNACaseAndControls) <- colnames(MicCaseAndControls)
colnames(predictedData) <- colnames(MicCaseAndControls)
colnames(predictedData2) <- colnames(MicCaseAndControls)

#Find DEGs
design <- model.matrix(~ 0+factor(c(1,1,1,1,1,1,2,2,2)))#6 control and 3 case for each chemical create a matrix for control samples column 1 is 1 and column 2 0 and for case samples col 1 is 0 and com 2 is 1
colnames(design) <- c("group1", "group2")#name columns of created matrix

#4 different sets to calculate DEG for each:
ExpressionSetMic = ExpressionSet(assayData = t(MicCaseAndControls))
ExpressionSetRNA = ExpressionSet(assayData = t(RNACaseAndControls))
ExpressionSetPredict = ExpressionSet(assayData = t(predictedData))
ExpressionSetPredict2 = ExpressionSet(assayData = t(predictedData2))

contrast.matrix <- makeContrasts(group1-group2, levels=design)#create a contrast matrix which is needed for DEG verification process
#assigned method = lmfit
fitMic <- lmFit(ExpressionSetMic, design)
fitRNA = lmFit(ExpressionSetRNA, design)
fitPredict = lmFit(ExpressionSetPredict, design)
fitPredict2 = lmFit(ExpressionSetPredict2, design)

fit2Mic <- contrasts.fit(fitMic, contrast.matrix)#calculate contrast matrix for data
fit2RNA<-contrasts.fit(fitRNA, contrast.matrix)
fit2Predict<-contrasts.fit(fitPredict, contrast.matrix)
fit2Predict2<-contrasts.fit(fitPredict2, contrast.matrix)

fit2Mic <- eBayes(fit2Mic)
fit2RNA <-eBayes(fit2RNA)
fit2Predict <-eBayes(fit2Predict)
fit2Predict2 <-eBayes(fit2Predict2)

resMic = topTable(fit2Mic,lfc = 1.5,p.value = 0.05, sort.by = "logFC", n=7000,adjust.method="BH")
resRNA = topTable(fit2RNA,lfc = 1.5,p.value = 0.05, sort.by = "logFC",n=7000,adjust.method="BH")
resPredict = topTable(fit2Predict, lfc = 1.5,p.value = 0.05, sort.by = "logFC",n=7000,adjust.method="BH")
resPredict2 = topTable(fit2Predict2, lfc = 1.5,p.value = 0.05, sort.by = "logFC",n=7000,adjust.method="BH")

MicName = rownames(resMic)
RNAName = rownames(resRNA)
PredictName = rownames(resPredict)
PredictName2 = rownames(resPredict2)

intersectMicAndRNA = intersect(RNAName, MicName)
intersectPredictAndRNA = intersect(RNAName, PredictName2)
moreThanMicroarray = sum(!(intersectPredictAndRNA %in% intersectMicAndRNA))
microarrayMore = sum(!(intersectMicAndRNA %in%intersectPredictAndRNA))
TP = length(intersectPredictAndRNA)
FP = length(PredictName2) - TP
#compare DEGs from different approaches(gene expressions from a microarray normalization method(RMA or MAS5), RNA-Seq, and proposed approach)
intersect1 = intersect(RNAName, MicName)
intersect2 = intersect(RNAName, PredictName)
intersect2_2 = intersect(RNAName, PredictName2)
intersect3 = intersect(PredictName , MicName)
intersect3_2 = intersect(PredictName2 , MicName)

names1 = intersect( MicName, PredictName2)
names2 = intersect( names1, RNAName)
mic2 = resMic[names2,]
RNA2 = resRNA[names2,]
predict2 = resPredict[names2,]
micVec = resMic[names2,]$logFC
micRNA = resRNA[names2,]$logFC
micPredict = resPredict[names2,]$logFC
micPredict2 = resPredict2[names2,]$logFC
cor(micVec,micRNA)
cor(micPredict,micRNA)
cor(micPredict2,micRNA)
sum(abs(micRNA -micVec))
sum(abs(micRNA -micPredict2))

micVecP = resMic[names2,]$adj.P.Val
micRNAP = resRNA[names2,]$adj.P.Val
micPredictP = resPredict2[names2,]$adj.P.Val
cor(micVecP,micRNAP)
cor(micPredictP,micRNAP)
sum(abs(micRNAP -micVecP))
sum(abs(micRNAP -micPredictP))

micVecPVal = resMic[names2,]$P.Value
micRNAPVal = resRNA[names2,]$P.Value
micPredictPVal = resPredict[names2,]$P.Value
cor(micVecPVal,micRNAPVal)
cor(micPredictPVal,micRNAPVal)
sum(abs(micRNAPVal -micVecPVal))
sum(abs(micRNAPVal -micPredictPVal))
