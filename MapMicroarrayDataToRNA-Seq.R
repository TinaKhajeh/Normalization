#Author : Tina Khajeh
#This code will map probes to genes and microarray to RNA-Seq. 
#eliminate probes to prebsets with _at postfix
#considering only 1 to 1 and 1 to many probe mapping
setwd("/Users/Tina/Desktop/Thesis Code")
library(rat2302.db)
 

#|-------------------------------------------------------------|
#|first select _at probsets with 11 probes for each set:  |
#|-------------------------------------------------------------|
library(rat2302probe)
probIds = rat2302probe$Probe.Set.Name
geneNameAndEntrezIdAndSymbol0 = select(rat2302.db, keys = probIds, keytype = "PROBEID" , columns = c("ENTREZID", "GENENAME","SYMBOL","REFSEQ","UNIGENE") )

allpostfixes = substr(probIds, start=8, stop=nchar(probIds))
uniquePostfixes = unique(allpostfixes)
#number of all probs:
length(probIds)
#number of all probsets:
length(unique(probIds) )
#number of all genes:
length(unique(geneNameAndEntrezIdAndSymbol0$ENTREZID))

#filter probes and get _at, _a_at, _s_at, _x_at probes:
indexes1 = (allpostfixes=="_at")
indexes2 = (allpostfixes=="_a_at")
indexes3 = (allpostfixes=="_s_at")
indexes4 = (allpostfixes=="_x_at")
atProbeSets = probIds[indexes1]
a_atProbeSets = probIds[indexes2]
s_atProbeSets = probIds[indexes3]
x_atProbeSets = probIds[indexes4]


selectedProbes = probIds[indexes1]
#number of probes with _at:
length(selectedProbes)
#number of probsets with _at:
length(unique(selectedProbes))
#not divisible by 11:
frequencyArray = table(selectedProbes)
TableNames = names(frequencyArray)
TableAsVector = as.vector(frequencyArray)
#recognize data which is not 11
indexOfNot11 = (TableAsVector!=11)
#Detect name of elements which are not 11
toDeleteNames = TableNames[indexOfNot11]
#count number of elements to delete which is equal to 19
numOfProbes = sum(TableAsVector[indexOfNot11])
#delete those from probesets:
newIndexes = ! selectedProbes %in% toDeleteNames
selectedProbesId2 = selectedProbes[newIndexes]
#number of probes with _at and 11 probes
length(selectedProbesId2)
#number of probe sets with _at and 11 probes
length(unique(selectedProbesId2))

#-------------------------------------------------------------
#|Delete Probesets which their gene symbol is Na             |
#|                                                          |
#-------------------------------------------------------------
geneSymbolAndProbID_at11 = select(rat2302.db , keys = selectedProbesId2 , keytype = "PROBEID", columns = c( "SYMBOL"))
probIdSymbol_at11 = cbind(geneSymbolAndProbID_at11$PROBEID ,geneSymbolAndProbID_at11$SYMBOL )
notNaGeneIndexes = !is.na(probIdSymbol_at11[,2])
data_at11Na = probIdSymbol_at11[notNaGeneIndexes,]
length(unique(data_at11Na[,1]))#num of probe sets



#-------------------------------------------------------------------|
#|delete probesets measure more than 1 gene                         |
#-------------------------------------------------------------------|
#find unique tuples:to delete probes of one probe set
uniqueTuples = unique(data_at11Na)
dim(uniqueTuples)
#recognize probIds which is repeated for many genes and delete it:
uniqueTuplesProbIds = uniqueTuples[,1]
tableProbeIds = table(uniqueTuplesProbIds)
sortedtableProbeIds = sort(tableProbeIds,decreasing = TRUE)
sum(sortedtableProbeIds>1)
#Conclusion : number of probes which check expression of more than 1 gene :813


idTemp = names(sortedtableProbeIds)
idVector = as.vector(sortedtableProbeIds)
indexesMoreThan1 = idVector>1
probSetsToManyGenes = idTemp[indexesMoreThan1]
probes_at11NA = data_at11Na[,1]
Index_at11NAFilter1 = ! probes_at11NA %in% probSetsToManyGenes
data_at11NaFilter1 = data_at11Na[Index_at11NAFilter1,]
PROBEID = data_at11NaFilter1[,1]
SYMBOL = data_at11NaFilter1[,2]
dataFrame_at11NAFilter1 = data.frame(PROBEID,SYMBOL)
write.table(dataFrame_at11NAFilter1,"Final_at11NaFilter1.txt",sep="\t",row.names=FALSE)

#-----------------------------------------------------------------------------|
#|Seperate data in to two data sets:                                          |
#|First : many to One mapping(one gene is measure with more than one probeset)|
#|Second : one to one mapping                                                 |
#|                                                                            |
#-----------------------------------------------------------------------------|
#data_at11NaFilter1
uniqueTuples_at11NAFilter2 = unique(data_at11NaFilter1)
geneColumn = uniqueTuples_at11NAFilter2[,2]
tableGeneColumn = table(geneColumn)
sortedTableGeneColumn = sort(tableGeneColumn,decreasing = TRUE)
tmpGeneNames = names(sortedTableGeneColumn)
geneRepeatsVector = as.vector(sortedTableGeneColumn)
geneIndexesMoreThan1 = geneRepeatsVector>1
sum(geneIndexesMoreThan1)
genesToManyProbeSet = tmpGeneNames[geneIndexesMoreThan1]
#one to one mapping
ProbeIdIndexGroup1 = !(uniqueTuples_at11NAFilter2[,2]  %in% genesToManyProbeSet)
data_at11NaFilter1Group1 = uniqueTuples_at11NAFilter2[ProbeIdIndexGroup1,]
PROBEID = data_at11NaFilter1Group1[,1]
SYMBOL = data_at11NaFilter1Group1[,2]
dataFrame_at11NAFilter2G1 = data.frame(PROBEID,SYMBOL)
write.table(dataFrame_at11NAFilter2G1,"Final_at11NaFilter2OneToOne.txt",sep="\t",row.names=FALSE)
#one to many mapping
ProbeIdIndexGroup2 = uniqueTuples_at11NAFilter2[,2]  %in% genesToManyProbeSet
data_at11NaFilter1Group2 = uniqueTuples_at11NAFilter2[ProbeIdIndexGroup2,]
PROBEID = data_at11NaFilter1Group2[,1]
SYMBOL = data_at11NaFilter1Group2[,2]
dataFrame_at11NAFilter2G2 = data.frame(PROBEID,SYMBOL)
write.table(dataFrame_at11NAFilter2G2,"Final_at11NaFilter2OneToMany.txt",sep="\t",row.names=FALSE)


#-----------------------------------------------------------------------------|
#Intersect genes between RNA-Seq and Microarray                               |
#-----------------------------------------------------------------------------|
RNAFile = read.table(file = "GSE55347_TGxSEQC_GeneExpressionIndex_Magic_20120831_116samples.txt",header=TRUE,fill = TRUE,sep = "\t")
RNAFile1 = cbind(RNAFile[,1:2], RNAFile[,10:125])
#delete genes that measured with 2 different anootation methods:

dataGeneAnnotation = cbind(as.character(RNAFile$Gene_Name),as.character(RNAFile$Annotation_Method))
uniqueDataGeneAnnotation = unique(dataGeneAnnotation)
Genefreq = sort(table(uniqueDataGeneAnnotation[,1]), decreasing = TRUE)
deleteNames = names(Genefreq[1:53])
indexesToKeep = !(RNAFile1$Gene_Name %in% deleteNames)
RNAFile2 = RNAFile1[indexesToKeep,]# deleted genes with two annotation method
RNAGenes = RNAFile2$Gene_Name
MicGroup1 = read.table(file = "Final_at11NaFilter2OneToOne.txt",header=TRUE,fill = TRUE,sep = "\t")
MicGroup2 = read.table(file = "Final_at11NaFilter2OneToMany.txt",header=TRUE,fill = TRUE,sep = "\t")
MicroarrayGenes = unique(c(as.character(MicGroup1$SYMBOL) , as.character(MicGroup2$SYMBOL)))
commonIndexes =RNAGenes %in% MicroarrayGenes 

RNAFile3 = RNAFile2[commonIndexes,]

#-----------------------------------------------------------------------------|
#Delete NA from RNA-Seq File                                                  |
#-----------------------------------------------------------------------------|
NAIndexes = is.na(RNAFile3)
rowSum = as.vector(rowSums(NAIndexes))
NAIndexes = rowSum>0
#Delete these genes :
RNAFile4 = RNAFile3[!NAIndexes,]
genes = RNAFile4[,1]
commonGENES = as.data.frame(genes)
write.table(commonGENES,"FinalCommonGenesWithoutNa.txt",sep="\t",row.names=FALSE)

#-----------------------------------------------------------------------------|
#Save Mapping between RNA-Seq and Microarray                                  |
#-----------------------------------------------------------------------------|
MicGroup1Symbol  = MicGroup1$SYMBOL
MicGroup1Index = MicGroup1Symbol %in% commonGENES$genes
newMicGroup1 = MicGroup1[MicGroup1Index,]

MicGroup2Symbol  = MicGroup2$SYMBOL
MicGroup2Index = MicGroup2Symbol %in% commonGENES$genes
newMicGroup2 = MicGroup2[MicGroup2Index,]

write.table(newMicGroup1,"FinalIntersectGenes_at11NaFilter2OneToOne.txt",sep="\t",row.names=FALSE)
write.table(newMicGroup2,"FinalIntersectGenes_at11NaFilter2OneToMany.txt",sep="\t",row.names=FALSE)

