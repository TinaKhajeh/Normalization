#Author: Tina Khajeh
#Detect genes with outlier values in microarray and RNA-Seq:
 
setwd("/Users/Tina/Desktop/Thesis Code")
CommonGenes = read.table(file = "CommonGenes.txt",header=TRUE,fill = TRUE,sep = "\t")

diffMatrix = abs(micScale - RNAScale)
numberOfOutliers = 200
outlierGenes = c()
for(i in 1:nrow(diffMatrix)) {
  selectedGenes = c()
  tmpVector = diffMatrix[i,]
  selectedGenes = order(tmpVector, decreasing=TRUE)[1:numberOfOutliers]#this will return index
  outlierGenes = rb
  
  
  ind(outlierGenes,selectedGenes ) 
}

intersectArray = intersect(outlierGenes[1,],outlierGenes[2,])
for (i in 3:nrow(diffMatrix)){
  intersectArray=intersect(outlierGenes[i,],intersectArray)
}
#find genes from obtained indexes
commonOutliers = length(intersectArray)
outlierGeneNames = CommonGenes$Genes[intersectArray] 
probMapGeneSymbol = read.table(file = "FilteredMappingProbeSymbols.txt",header=TRUE,fill = TRUE,sep = "\t")
indexOutlierGenes = which(probMapGeneSymbol$SYMBOL %in% outlierGeneNames)
outlierProbeNames200 = probMapGeneSymbol$PROBEID[indexOutlierGenes]


#save data to a file
save(outlierProbeNames200, file = "outlierProbeName200.RData")


#-----------------------------------------------------------------------------|
#plot microarray and RNA for two samples with outlier genes:                  |
#-----------------------------------------------------------------------------|
#
#*********Sample 1:************
mean(micScale)
mean(RNAScale)
plot(micScale[1,], RNAScale[1,])
#all 500:
tmpselectedGenes500 = order(diffMatrix[1,], decreasing=TRUE)[1:500]
points(micScale[1,tmpselectedGenes500],RNAScale[1,tmpselectedGenes500],col = 'blue')
#outliers:
points(micScale[1,intersectArray],RNAScale[1,intersectArray],col = 'red')

#*********Sample50**********
mean(micScale)
mean(RNAScale)
plot(micScale[50,], RNAScale[50,])
#all 500:
tmpselectedGenes500 = order(diffMatrix[50,], decreasing=TRUE)[1:500]
points(micScale[50,tmpselectedGenes500],RNAScale[50,tmpselectedGenes500],col = 'blue')
#outliers:
points(micScale[50,intersectArray],RNAScale[50,intersectArray],col = 'red')


