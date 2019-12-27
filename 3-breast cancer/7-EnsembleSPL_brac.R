##----------------------------------
##  Code Version 1.0
##  This breast cancer code for self-paced learning with L1 penalty embeded in the Ensemble-based framework.
##  Ensemble_SPL
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

library(glmnet)
library(caret)

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
# 2.Ensemble_SPL
##--------------
ensemblePanels <- lapply(train_data, function(i){
  cvfit_ensem <- cv.glmnet(i, train_label, alpha=0.9,family="multinomial",type.multinomial="grouped", nfolds = 4, lambda = seq(0.045,0.05,by=0.005))
})


## Initializing parameters
sample.select <- list()
sample.add <- list()
Times <- 25
Iter_num = 100
for(i in 1: length(unique(train_label))){
  sample.select[[i]] <- 4
  sample.add[[i]] <- ceiling(length(which(train_label==i))/Times)
}

cat("Starting training.\n", sep = "")

ensemblePanels <- lapply(train_data, function(i){
  cvfit_ensem <- selfpaced.EN.multiclass(X_train = i,
                                         Y_train = train_label,
                                         lambda = sample.select,
                                         uplambda = sample.add,
                                         Iter_num = Iter_num)
})

ensembleValiPredction <- mapply(function(cvfit,x){
  valprediction <- as.numeric(predict(cvfit,x,type="class",s="lambda.min"))
}, cvfit = ensemblePanels, x = train_data)

ensembleTestPredction <- mapply(function(cvfit,x){
  tsprediction <- as.numeric(predict(cvfit,x,type="class",s="lambda.min"))
}, cvfit = ensemblePanels, x = test_data)

ensembleCoef <- lapply(ensemblePanels, function(i){
  coefprediction <- predict(i,type="coefficients",s="lambda.min")
})

ensembleCoef.idx <- lapply(ensembleCoef, function(i){
  coef.idx <- which(i$`1`[-1]!=0)
})

ensembleCoef.name <- lapply(ensembleCoef, function(i){
  coef.idx <- which(i$`1`!=0)[-1]
  coef.name <- rownames(i$`1`)[coef.idx]
})

valprediction <- evaluate.ensemble(ensembleValiPredction, train_label)
tsprediction <- evaluate.ensemble(ensembleTestPredction, test_label)

##----------------------------
# 3. Evaluate the performance
##----------------------------
conMatrix.train <- confusionMatrix(as.factor(train_label),as.factor(valprediction))
conMatrix.test <- confusionMatrix(as.factor(test_label),as.factor(tsprediction),mode = "prec_recall")


##----------------------
## feature
##----------------------
coef.name.mRNA <- ensembleCoef.name$mrna
coef.value.mRNA <- ensembleCoef$mrna$`1`[which(ensembleCoef$mrna$`1`!=0)][-1]
coef.mRNA <- cbind(coef.name.mRNA,coef.value.mRNA) 

coef.name.miRNA <- ensembleCoef.name$mirna
coef.value.miRNA <- ensembleCoef$mirna$`1`[which(ensembleCoef$mirna$`1`!=0)][-1]
coef.miRNA <- cbind(coef.name.miRNA,coef.value.miRNA) 

coef.name.cpg <- ensembleCoef.name$meth
coef.value.cpg <- ensembleCoef$meth$`1`[which(ensembleCoef$meth$`1`!=0)][-1]
coef.cpg <- cbind(coef.name.cpg,coef.value.cpg) 

coef.name <- c(coef.name.mRNA, coef.name.miRNA, coef.name.cpg)

## calculate features ##
num_feature <- matrix(0, nrow = 1, ncol = 3)
num_feature[,1] <- length(grep("mrna", coef.name))
num_feature[,2] <- length(grep("mirna", coef.name))
num_feature[,3] <- length(grep("meth", coef.name))



perf.Ensemble_SPL <- list("coef.mRNA" = coef.mRNA,
                          "coef.miRNA" = coef.miRNA,
                          "coef.cpg" = coef.cpg,
                          "feature.num" = num_feature,
                          "Perf.Train" = conMatrix.train,
                          "Perf.Test" = conMatrix.test)

