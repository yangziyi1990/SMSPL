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
load("1-data/1-datatrainTestDatasetsNormalized.RDATA")

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

##------------------------
## 2. concatenation_sPLSDA
##------------------------
combined_train_data <- do.call(cbind, train_data)
combined_test_data <- do.call(cbind, test_data)

y_train <- as.factor(train_label)

fit.splsda <- splsda(combined_train_data,y_train, K=3, eta=0.75, scale.x=FALSE)
valprediction <- predict(fit.splsda, newx = combined_train_data)
tsprediction <- predict(fit.splsda, newx = combined_test_data)
coefmatrix <- coef(fit.splsda)

valprediction <- as.numeric(as.character(valprediction))
tsprediction <- as.numeric(as.character(tsprediction))

## calculate features ##
coef.idx <- which(coefmatrix!=0)
coef.value <- coefmatrix[which(coefmatrix!=0)]
coef.name <- colnames(combined_train_data)[coef.idx]

coef <- cbind(coef.name, coef.value)
coef.mRNA <- coef[grep("mrna", coef),]
coef.miRNA <- coef[grep("mirna", coef),]
coef.cpg <- coef[grep("meth", coef),]

num_feature <- matrix(0, nrow = 1, ncol = 3)
num_feature[,1] <- length(grep("mrna", coef.name))
num_feature[,2] <- length(grep("mirna", coef.name))
num_feature[,3] <- length(grep("meth", coef.name))

##----------------------------
# 3. Evaluate the performance
##----------------------------

conMatrix.train <- confusionMatrix(as.factor(train_label),as.factor(valprediction))
conMatrix.test <- confusionMatrix(as.factor(test_label),as.factor(tsprediction),mode = "prec_recall")

perf.concatenation_sPLSDA <- list("Feature" = coef.name,
                                  "Feature.idx" = coef.idx,
                                  "Perf.Train" = conMatrix.train,
                                  "Perf.Test" = conMatrix.test)


