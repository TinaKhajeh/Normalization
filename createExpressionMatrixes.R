#Author: Tina Khajeh
#this code will create matrixes for microarray and RNA-Seq expressions
 
setwd("/Users/Tina/Desktop/Thesis Code")
library(GEOquery)
library(seqinr)
library(affy)
library(Biostrings)
library(rat2302probe)

CommonGenes = read.table(file = "FinalCommonGenesWithoutNa.txt",header=TRUE,fill = TRUE,sep = "\t")
#MAA5:
tmp=getGEO(GEO = "GSE47875", destdir = "/Users/Tina/Desktop/Thesis Code",AnnotGPL = TRUE)
microExpressionMas5 = tmp[1]$GSE47875_series_matrix.txt.gz# MAS5 normalization
expressionMatrixMas5 = exprs(microExpressionMas5)
##RMA:
data=ReadAffy(celfile.path="/Users/Tina/Desktop/Thesis Code/data/") # Import expression raw data and stores them as AffyBatch object.
eset_rma <- rma(data) # Normalizes the data with 'rma' function 
expressionMatrixRMA = exprs(eset_rma)
colNames = colnames(expressionMatrixRMA)
tmpcolNames = substr(colNames, start=0, stop=10)
colnames(expressionMatrixRMA) = tmpcolNames

save(expressionMatrixMas5, expressionMatrixRMA, file = "FinalRMA&Mas5MicroarrayExpression.RData")
#load("FinalRMA&Mas5MicroarrayExpression.RData")
#rows are not in the same order:
rowIndexes = match(rownames(expressionMatrixRMA), rownames(expressionMatrixMas5))
expressionMatrixRMA2 = expressionMatrixRMA[rowIndexes,]

samples = getGEO(GEO = "GSE47875", destdir = "/Users/Tina/Desktop/Thesis Code",GSEMatrix = FALSE)
sampleList = GSMList(samples)
namesList = c()#GSMList
longNames = c()
for( i in 1:length(sampleList))
{
  metaSampleList =Meta(sampleList[[i]])
  longNames = c(longNames,metaSampleList$title)
  name = names(sampleList)[i]
  namesList = c(namesList,name)
}

namesAndGSMMicroarray = data.frame(longNames,namesList)#contains long and short name
#probMapGeneSymbol = read.table(file = "Filtered1ProbeSymbolsToManyGenes.txt",header=TRUE,fill = TRUE,sep = "\t")
probMapGeneSymbolOneToOne = read.table(file = "FinalIntersectGenes_at11NaFilter2OneToOne.txt",header=TRUE,fill = TRUE,sep = "\t")
probMapGeneSymbolOneToMany = read.table(file = "FinalIntersectGenes_at11NaFilter2OneToMany.txt",header=TRUE,fill = TRUE,sep = "\t")

RNAFile = read.table(file = "GSE55347_TGxSEQC_GeneExpressionIndex_Magic_20120831_116samples.txt",header=TRUE,fill = TRUE,sep = "\t")
SamplesTable = read.table(file = "Treated to control mapping for SEQC samples-with Run_s.txt",header=TRUE,fill = TRUE,sep = "\t")
#eliminate to training and test samples:
indexes = SamplesTable$Original.SET.Annotation %in% c("TEST SET", "TRAINING SET")
newSampleTable = SamplesTable[indexes,]#train and test

#Corresponding Samples
micSample = c()
RNASample = c()
for (sampleCount in 1:nrow(newSampleTable)) {#sample is in file for mapping treated and control
  sample = newSampleTable[sampleCount,]
  sampleMic = as.character(sample$Microarray.ID)
  sampleRNA = as.character(sample$RNA.seq.ID)
  tmpIndex = (namesAndGSMMicroarray$longNames ==as.character(sampleMic))
  sampleMicShort = as.character(namesAndGSMMicroarray$namesList[tmpIndex])
  micSample = c(micSample,sampleMicShort)
  RNASample = c(RNASample,sampleRNA)
  #set name for columns and rows
}

#firts: oneToone Mapping
oneToOneMappingGenes = as.character(unique(probMapGeneSymbolOneToOne$SYMBOL))
oneToManyMappingGenes = as.character(unique(probMapGeneSymbolOneToMany$SYMBOL))
micCol= c()
RNAColIndex = c()
RNAColGeneName = c()
for(counterIndex in 1:(length(oneToOneMappingGenes))){
  gene = as.character(probMapGeneSymbolOneToOne$SYMBOL[counterIndex])
  probeId = as.character(probMapGeneSymbolOneToMany$PROBEID[counterIndex])
  micCol = c(micCol,probeId)
  #micExpression = expressionMatrix[probeId,sampleMicShort]
  RNAIndex = which(RNAFile$Gene_Name == as.character(gene))
  if(length(RNAIndex)!=1){
    print(length(RNAIndex))
    print(gene)
    print(counterIndex)
  }
  RNAgeneName = as.character(RNAFile$Gene_Name[RNAIndex])
  tmpIndex = rep(RNAIndex,length(probeId))
  tmpRNAgeneName = rep(RNAgeneName, length(probeId))
  RNAColIndex = c(RNAColIndex,tmpIndex)
  RNAColGeneName = c(RNAColGeneName, tmpRNAgeneName)
}

oneToOneMas5Matrix = expressionMatrixMas5[micCol,micSample]
oneToOneMas5Matrix = t(oneToOneMas5Matrix)
rownames(oneToOneMas5Matrix) = micSample
colnames(oneToOneMas5Matrix) = micCol

oneToOneRMAMatrix = expressionMatrixRMA[micCol,micSample]
oneToOneRMAMatrix = t(oneToOneRMAMatrix)
rownames(oneToOneRMAMatrix) = micSample
colnames(oneToOneRMAMatrix) = micCol

RNAOneToOneMatrix = RNAFile[RNAColIndex, RNASample]
RNAOneToOneMatrix = t(RNAOneToOneMatrix)
rownames(RNAOneToOneMatrix) = RNASample
colnames(RNAOneToOneMatrix) = RNAColGeneName


micMatrixOneToOne = expressionMatrix[micCol,micSample]
micMatrix = t(micMatrix)
rownames(micMatrix) = micSample
colnames(micMatrix) = micCol
RNAMatrix = RNAFile[RNAColIndex, RNASample]
RNAMatrix = t(RNAMatrix)
rownames(RNAMatrix) = RNASample
colnames(RNAMatrix) = RNAColGeneName
micMatrix1_2 = micMatrix
RNAMatrix1_2 = RNAMatrix
#save(micMatrix1_2, RNAMatrix1_2, file = "2MatrixFilteredGenes1_2.RData")
#write.table(micMatrix1_2, "MicExpression1_2.txt",sep = "\t")
#write.table(RNAMatrix1_2, "RNAExpression1_2.txt",sep = "\t")
cols = colnames(micMatrix1_2)
rows = rownames(micMatrix1_2)
RMAMicMatrix = t(expressionMatrixRMA[cols,rows])
plot(micMatrix1_2[1,],RNAMatrix1_2[1,] )
plot(RMAMicMatrix[1,], RNAMatrix1_2[1,])

#find columns with Na in RNA:
indexes = is.na(RNAMatrix1_2)
indexes = t(indexes)
indexesRowSum = rowSums(indexes)
toDelete = (indexesRowSum>0)
keepGenes = !(toDelete)
newRNAFile = RNAMatrix1_2[,keepGenes]
newMas5 = micMatrix1_2[,keepGenes]
newRMA = RMAMicMatrix[,keepGenes]
cor(newMas5[1,],newRNAFile[1,])
cor(newMas5[1,],newRNAFile[1,])
#save(newRNAFile, newMas5, newRMA, file = "3DataAfterPreprocessing.RData")
load("3DataAfterPreprocessing.RData")

#Scale Data:
newMas5Scale = scaleFunction(newMas5,100)
newRMAScale = scaleFunction(newRMA,100)
newRNAFileScale = scaleFunction(newRNAFile,100)

#remove mean:
meanDiffRMA = mean(newRMAScale) - mean(newRNAFileScale)
meanDiffMas5 = mean(newMas5Scale) - mean(newRNAFileScale)
newMas5ScaleremoveMean = newMas5Scale - meanDiffMas5
newRMAScaleremoveMean = newRMAScale - meanDiffRMA
save(newRNAFileScale, newMas5ScaleremoveMean, newRMAScaleremoveMean, file = "3DataScaledAndRemoveMean.RData")

#recognize repeated probes:
probeIds = colnames(newMas5ScaleremoveMean)
geneNames = colnames(newRNAFileScale)
mapping = unique(cbind(probeIds, geneNames))
geneColumn = mapping[,2]
geneColumnTable = sort(table(geneColumn), decreasing = TRUE)
geneColumnNames = names(geneColumnTable)
geneColumnFrequency = as.vector(geneColumnTable)
indexesroup2 = (geneColumnFrequency>1)
geneGroup2 = geneColumnNames[indexesroup2]
geneGroup1 = geneColumnNames[!indexesroup2]
group2MatrixColumn = geneNames %in% geneGroup2
group1MatrixColumn = geneNames %in% geneGroup1
newMas5Group1 = newMas5ScaleremoveMean[,group1MatrixColumn]
newMas5Group2 = newMas5ScaleremoveMean[,group2MatrixColumn]
newRMAGroup1 = newRMAScaleremoveMean[,group1MatrixColumn]
newRMAGroup2 = newRMAScaleremoveMean[,group2MatrixColumn]
newRNAFileGroup1 = newRNAFileScale[,group1MatrixColumn]
newRNAFileGroup2 = newRNAFileScale[,group2MatrixColumn]
plot(newMas5[1,], newRNAFile[1,])
plot(newMas5Group1[1,], newRNAFileGroup1[1,])
plot(newMas5Group2[1,], newRNAFileGroup2[1,])

plot(newRMA[1,], newRNAFile[1,])
plot(newRMAGroup1[1,], newRNAFileGroup1[1,])
plot(newRMAGroup2[1,], newRNAFileGroup2[1,])
save(newMas5Group1, newMas5Group2,newRMAGroup1,newRMAGroup2, newRNAFileGroup1,newRNAFileGroup2, file = "6GroupData.RData")
#load("6GroupData.RData")
save(newMas5Group1, newMas5Group2,newRMAGroup1,newRMAGroup2, newRNAFileGroup1,newRNAFileGroup2, file = "6GroupData_2.RData")
#load("6GroupData_2.RData")
