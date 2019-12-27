##----------------------------------
##  Code Version 1.0
##  This breast cancer code for sPLSDA embeded in the concatenation-based framework.
##  sPLSDA
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

library(glmnet)
library(caret)
library(spls)

##--------------
# 1. Load data
##--------------
source("Function_performance.R")
load("/1-data/1-datatrainTestDatasetsNormalized.RDATA")

## Preparing data
train_data <- list("mrna" = mrnaTest0,
                   "mirna" = mirnaTest0,
                   "meth" = methTest0)
colnames(train_data$mrna) <- paste("mrna", colnames(train_data$mrna), sep = "_")
colnames(train_data$mirna) <- paste("mirna", colnames(train_data$mirna), sep = "_")
colnames(train_data$meth) <- paste("meth", colnames(train_data$meth), sep = "_")
train_group <- pam50Test0$Call


test_data <- list("mrna" = mrnaTrain0,
                  "mirna" = mirnaTrain0,
                  "meth" = methTrain0)
colnames(test_data$mrna) <- paste("mrna", colnames(test_data$mrna), sep = "_")
colnames(test_data$mirna) <- paste("mirna", colnames(test_data$mirna), sep = "_")
colnames(test_data$meth) <- paste("meth", colnames(test_data$meth), sep = "_")
test_group <- pam50Train0$Call

train_label <-as.numeric(train_group)
test_label <-as.numeric(test_group)

##--------------
# 2.Ensemble_EN
##--------------
J=3
y_train <- as.factor(train_label)
y_test <- as.factor(test_label)

## Initialize 
coefmatrix <- list()
valpredmatrix <- matrix(0, nrow = length(y_train), ncol = J)
evlpredmatrix <- matrix(0, nrow = length(y_test), ncol = J)

## Ensemble_sPLSDA
ensemblePanels <- lapply(train_data, function(i){
  cvfit_ensem <- splsda(i,y_train, K=3, eta=0.85, scale.x=FALSE)
})

ensembleValiPredction <- mapply(function(cvfit,x){
  valprediction <- predict(cvfit, newx = x)
}, cvfit = ensemblePanels, x = train_data)

ensembleTestPredction <- mapply(function(cvfit,x){
  tsprediction <- predict(cvfit, newx = x)
}, cvfit = ensemblePanels, x = test_data)

ensembleCoef <- lapply(ensemblePanels, function(i){
  coefprediction <- coef(i)[,1]
})

ensembleValiPredction <- as.data.frame(ensembleValiPredction)
ensembleTestPredction <- as.data.frame(ensembleTestPredction)

valprediction <- evaluate.ensemble.sPLSDA(ensembleValiPredction, train_label)
tsprediction <- evaluate.ensemble.sPLSDA(ensembleTestPredction, test_label)


##----------------------
## feature
##----------------------
coef.idx.mRNA <- which(ensembleCoef$mrna!=0)
coef.name.mRNA <- colnames(train_data[[1]])[coef.idx.mRNA]
coef.value.mRNA <- ensembleCoef$mrna[coef.idx.mRNA]
coef.mRNA <- cbind(coef.name.mRNA,coef.value.mRNA) 

coef.idx.miRNA <- which(ensembleCoef$mirna!=0)
coef.name.miRNA <- colnames(train_data[[2]])[coef.idx.miRNA]
coef.value.miRNA <- ensembleCoef$mirna[coef.idx.miRNA]
coef.miRNA <- cbind(coef.name.miRNA,coef.value.miRNA) 

coef.idx.cpg <- which(ensembleCoef$meth!=0)
coef.name.cpg <- colnames(train_data[[3]])[coef.idx.cpg]
coef.value.cpg <- ensembleCoef$meth[coef.idx.cpg]
coef.cpg <- cbind(coef.name.cpg,coef.value.cpg) 

coef.name <- c(coef.name.mRNA, coef.name.miRNA, coef.name.cpg)

## calculate features ##
num_feature <- matrix(0, nrow = 1, ncol = 3)
num_feature[,1] <- length(grep("mrna", coef.name))
num_feature[,2] <- length(grep("mirna", coef.name))
num_feature[,3] <- length(grep("meth", coef.name))


##----------------------------
# 3. Evaluate the performance
##----------------------------
conMatrix.train <- confusionMatrix(as.factor(train_label),as.factor(valprediction))
conMatrix.test <- confusionMatrix(as.factor(test_label),as.factor(tsprediction),mode = "prec_recall")

perf.Ensemble.sPLSDA <- list("Feature" = coef.name,
                             "Perf.Train" = conMatrix.train,
                             "Perf.Test" = conMatrix.test)

