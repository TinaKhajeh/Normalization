#Author: Tina Khajeh
#use cell files to acquire microarray data
 
setwd("/Users/Tina/Desktop/Thesis Code")
library(GEOquery)
library(affy)

#Download the CEL file package for this dataset
getGEOSuppFiles("GSE47875")

#Unpack the CEL files:
setwd("/Users/Tina/Desktop/Thesis Code/GSE47875")
untar("GSE47875_RAW.tar", exdir="data")

celData=ReadAffy(celfile.path = "/Users/Tina/Desktop/Thesis Code/data")
normalizedData = rma(celData)
x1= exprs(normalizedData)
normalized2Data = mas5(celData)
logNormal2 = log2(exprs(normalized2Data))
x2 = logNormal2
x[1,1:10]
expressionMatrix[1,1:10]
x2[1,1:10]
#save 2 normalized data RMA and MAS5:
rmaNormal = normalizedData
MAS5Normal = normalized2Data
saveData = c(rmaNormal , MAS5Normal)
save(saveData , file= "NormalizedRmaAndMas5.RData")

pmProbesIntensities = probes(celData,which="pm")
rowNames = row.names(pmProbesIntensities)
save (pmProbesIntensities, file= "Intensities.RData")
