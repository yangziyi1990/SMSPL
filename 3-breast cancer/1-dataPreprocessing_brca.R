#################################################
# Code version V1.0
# 1-dataPreprocessing_brca.R
# Date: December 26, 2019
# updated: December 26, 2019
# Author: Zi-Yi Yang
#
#################################################
WhereAmI <- "D:/path"

library(tidyverse)

###----------------------------
#
# Import datasets
#
###----------------------------
load(paste0(WhereAmI, "/datatrainTestDatasets.RDATA"))

trainSubj <- rownames(pam50Train)[pam50Train$Call != "Normal"]
testSubj <- rownames(pam50Test)[pam50Test$Call != "Normal"]
###################################################################
#
# 1) Normalize mRNA dataset : 1.1 filter 1.2 select common genes (train & test)
# 1.3 Normalize 1.4 remove PAM50 genes (options)
#
###################################################################
mrnaTrain2 <- apply(mrnaTrain, 2, as.numeric)  
rownames(mrnaTrain2) <- rownames(mrnaTrain)
mrnaTest2 <- apply(mrnaTest, 2, as.numeric)   
rownames(mrnaTest2) <- rownames(mrnaTest)

## remove genes with more than 70% zeros
mrnaTrain2 <- mrnaTrain2[rowSums(mrnaTrain2 != 0)/ncol(mrnaTrain2) >= 0.70, ]  # 16913   387
mrnaTest2 <- mrnaTest2[rowSums(mrnaTest2 != 0)/ncol(mrnaTest2) >= 0.70, ]  # 17045   638

## filter both datasets to have the common set of genes
comGenes <- intersect(rownames(mrnaTrain2), rownames(mrnaTest2))  
mrnaTrain2 <- mrnaTrain2[comGenes, ]
mrnaTest2 <- mrnaTest2[comGenes, ]
dim(mrnaTrain2); # 16901   387
dim(mrnaTest2); # 16901   638

## Normalize entire dataset
lib.size <- colSums(mrnaTrain2)
mrnaTrain3 <- log2(t(mrnaTrain2+0.5)/(lib.size+1) * 1e+06)
lib.size <- colSums(mrnaTest2)
mrnaTest3 <- log2(t(mrnaTest2+0.5)/(lib.size+1) * 1e+06)

all(colnames(mrnaTrain3) == colnames(mrnaTest3))
dim(mrnaTrain3); #387 16901
dim(mrnaTest3); #638 16901

## remove PAM50 genes
pam50genes <- read.delim(paste0(WhereAmI, "/pam50_annotation.txt"))
as.character(pam50genes$GeneName)

mrnaTrain4 <- mrnaTrain3[, setdiff(colnames(mrnaTrain3), as.character(pam50genes$GeneName))]
mrnaTest4 <- mrnaTest3[, setdiff(colnames(mrnaTest3), as.character(pam50genes$GeneName))]
dim(mrnaTrain4);  # 387 16851
dim(mrnaTest4);  # 638 16851

## train and test datasets
mrnaTrain5 <- mrnaTrain4[trainSubj, ]   # 379 16851
mrnaTest5 <- mrnaTest4[testSubj, colnames(mrnaTrain5)]  # 610 16851


## missing data?
sum(colSums(is.na(mrnaTrain5)))
sum(colSums(is.na(mrnaTest5)))
dim(mrnaTrain5); # 379 16851
dim(mrnaTest5); # 610 16851

###################################################################
#
# 2) Normalize miRNA dataset : 2.1 filter 2.2 select common genes (train & test)
# 2.3 Normalize
#
###################################################################

## remove miRNA with more than 70% zeros
mirnaTrain <- mirnaTrain[rowSums(mirnaTrain != 0)/ncol(mirnaTrain) >= 0.70, ]  # 350 387
mirnaTest <- mirnaTest[rowSums(mirnaTest != 0)/ncol(mirnaTest) >= 0.70, ] # 383 638

## filter both datasets to have the common set of miRNAs
comGenes <- intersect(rownames(mirnaTrain), rownames(mirnaTest))
mirnaTrain <- mirnaTrain[comGenes, ]
mirnaTest <- mirnaTest[comGenes, ]
dim(mirnaTrain);  # 349 387
dim(mirnaTest); #  349 638

## Normalize entire dataset
lib.size <- colSums(mirnaTrain)
mirnaTrain2 <- log2(t(mirnaTrain+0.5)/(lib.size+1) * 1e+06)
lib.size <- colSums(mirnaTest)
mirnaTest2 <- log2(t(mirnaTest+0.5)/(lib.size+1) * 1e+06)

## select trainSubj samples
mirnaTrain3 <- mirnaTrain2[trainSubj, ]
mirnaTest3 <- mirnaTest2[testSubj, colnames(mirnaTrain3)]
all(colnames(mirnaTrain3) == colnames(mirnaTest3))

## missing data?
sum(colSums(is.na(mirnaTrain3)))
sum(colSums(is.na(mirnaTest3)))
dim(mirnaTrain3); # 379 349
dim(mirnaTest3); #  610 349

###################################################################
#
# 3) Normalize methylation dataset
#
###################################################################
## summarise to gene-symbols
all(rownames(methTrain) == rownames(methAnnotation))
## Train
methTrain2 <- methTrain %>% 
  as.data.frame %>% 
  mutate(genSym = methAnnotation$Gene_Symbol) %>% 
  group_by(genSym) %>% 
  summarise_all(funs(mean))
methTrain3 <- as.matrix(methTrain2[, -1])
rownames(methTrain3) <- methTrain2$genSym  # 14270   387
## Test
methTest2 <- methTest %>% 
  as.data.frame %>% 
  mutate(genSym = methAnnotation$Gene_Symbol) %>% 
  group_by(genSym) %>% 
  summarise_all(funs(mean))
methTest3 <- as.matrix(methTest2[, -1])
rownames(methTest3) <- methTest2$genSym  # 14270   638

## remove cps with more than 70% missing data
# methTrain4 <- methTrain3[rowSums(!is.na(methTrain3))/ncol(methTrain3) >= 0.70, ]
# methTest4 <- methTest3[rowSums(!is.na(methTest3))/ncol(methTest3) >= 0.70, ]

## for now remove variables with any missing data
methTrain4 <- na.omit(methTrain3)
methTest4 <- na.omit(methTest3)

commonCpGs <- intersect(rownames(methTrain4), rownames(methTest4))
methTrain5 <- t(methTrain4[commonCpGs, trainSubj])
methTest5 <- t(methTest4[commonCpGs, testSubj])

###################################################################
#
# 4) Proteomics 4.1 check annoation 4.2 compare protein names with gene ID
#
###################################################################
## are all proteins in the annotation file?
missingProt <- setdiff(rownames(protTrain), rownames(protAnnotation))

## add missing proteins to annotation file   ## according to the need to modify
missingProteins <- data.frame(Composite.Element.REF = missingProt,
                              Gene.Name = c("EIF4EBP1","ANLN","ANXA1","MAPK8","KRAS","LBK1","PEA15","RBM3","SCD1"))
rownames(missingProteins) <- missingProteins$Composite.Element.REF

protAnnotation <- rbind(protAnnotation, missingProteins)
## are all proteins in the annotation file?
setdiff(rownames(protTrain), rownames(protAnnotation))  # yes  character(0)

## average proteins with the same gene ID
protTrain2 <- protTrain[, trainSubj] %>%
  as.data.frame %>% 
  mutate(genSym = as.character(protAnnotation[rownames(protTrain), "Gene.Name"])) %>% 
  group_by(genSym) %>% 
  summarise_all(funs(mean))
protTrain3 <- t(as.matrix(protTrain2[, -1]))
colnames(protTrain3) <- protTrain2$genSym

## any missing data in the proteomics dataset?
sum(colSums(is.na(protTrain)))
dim(protTrain3)  # 379 115

###################################################################
#
# Finalize datasets
#
###################################################################
clinVar <- rownames(na.omit(t(rbind(clinTrain[trainSubj,], clinTest[testSubj,]))))
###----------- Training dataset
clinTrain0 <- clinTrain[trainSubj, clinVar]
mrnaTrain0 <- mrnaTrain5
mirnaTrain0 <- mirnaTrain3
methTrain0 <- methTrain5
protTrain0 <- protTrain3
pam50Train0 <- pam50Train[trainSubj, ]
dim(clinTrain0); #  379 217
dim(mrnaTrain0); # 379 16851
dim(mirnaTrain0); # 379 349
dim(methTrain0);  # 379 9482
dim(protTrain0); # 379 115

###----------- Training dataset
clinTest0 <- clinTest[testSubj, clinVar]
mrnaTest0 <- mrnaTest5
mirnaTest0 <- mirnaTest3
methTest0 <- methTest5
pam50Test0 <- pam50Test[testSubj, ]
dim(clinTest0); # 610 217
dim(mrnaTest0); # 610 16851
dim(mirnaTest0); # 610 349
dim(methTest0); # 610 9482


###################################################################
#
# Data distributions
#
###################################################################
pdf(paste0(WhereAmI, "/figure/OverlapExpressionData_brca.pdf"))
par(mfrow = c(2, 2))
hist(mrnaTest0, col=rgb(1,0,0,0.5), main="mRNA", xlab="log2 cpm")
hist(mrnaTrain0, col=rgb(0,0,1,0.5), add=T)
box()
hist(mirnaTest0, col=rgb(1,0,0,0.5), main="miRNA", xlab="log2 cpm")
hist(mirnaTrain0, col=rgb(0,0,1,0.5), add=T)
box()
hist(as.matrix(methTest0), col=rgb(1,0,0,0.5), main="CpG", xlab="beta values")
hist(as.matrix(methTrain0), col=rgb(0,0,1,0.5), add=T)
box()
hist(as.matrix(protTrain0), col=rgb(0,0,1,0.5), main="Proteomics", xlab="centered data")
legend("topright", c("Train", "Test"),
       col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), pch = 19,
       bty = "n")
box()
dev.off()

###################################################################
#
# Save to file
#
###################################################################
save(clinTrain0=clinTrain0, mrnaTrain0=mrnaTrain0, mirnaTrain0=mirnaTrain0, protTrain0=protTrain0, 
     methTrain0=methTrain0, pam50Train0=pam50Train0,
     clinTest0=clinTest0, mrnaTest0=mrnaTest0, mirnaTest0=mirnaTest0, methTest0=methTest0, pam50Test0=pam50Test0,
     file = paste0(WhereAmI, "trainTestDatasetsNormalized.RDATA"))


