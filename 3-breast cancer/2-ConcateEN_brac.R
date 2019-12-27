##----------------------------------
##  Code Version 1.0
##  This breast cancer code for logistic regression with Elastic net penalty embeded in the concatenation-based framework.
##  Concate_EN
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
library(glmnet)
library(caret)

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

##--------------
# 2. Concatenation
##--------------
combined_train_data <- do.call(cbind, train_data)
combined_test_data <- do.call(cbind, test_data)

cvfit <- cv.glmnet(combined_train_data,train_label,alpha=0.9, family="multinomial", type.multinomial="grouped") 
valprediction <- predict(cvfit,combined_train_data,type="class",s="lambda.min")  # train
tsprediction <- predict(cvfit,combined_test_data,type="class",s="lambda.min") # test

valprediction <- as.numeric(valprediction)
tsprediction <- as.numeric(tsprediction)

coefprediction <- predict(cvfit,type="coefficients",s="lambda.min")  # coef
coef.idx <- which(coefprediction$`1`!=0)[-1]
coef.name <- rownames(coefprediction$`1`[-1])[coef.idx]

##----------------------------
# 3. Evaluate the performance
##----------------------------

conMatrix.train <- confusionMatrix(as.factor(train_label),as.factor(valprediction))
conMatrix.test <- confusionMatrix(as.factor(test_label),as.factor(tsprediction),mode = "prec_recall")

perf.ConcatenationEN <- list("Feature" = coef.name,
                             "Feature.idx" = coef.idx,
                             "Perf.Train" = conMatrix.train,
                             "Perf.Test" = conMatrix.test)


