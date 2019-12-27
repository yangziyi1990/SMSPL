##----------------------------------
##  Code Version 1.0
##  This breast cancer code for self-paced learning with L1 penalty embeded in the concatenation-based framework.
##  Concate_SPL
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

##--------------------
# 2. Concatenation_SPL
##--------------------
combined_train_data <- do.call(cbind, train_data)
combined_test_data <- do.call(cbind, test_data)

## Initializing parameters
sample.select <- list()
sample.add <- list()
Times <- 25
Iter_num = 100
for(i in 1: length(unique(train_label))){
  sample.select[[i]] <- 4
  sample.add[[i]] <- ceiling(length(which(train_label==i))/Times)
}

results <- selfpaced.muliticlass(X_train = combined_train_data, 
                                 Y_train = train_label,
                                 X_test = combined_test_data, 
                                 Y_test = test_label,
                                 lambda = sample.select,
                                 uplambda = sample.add,
                                 Iter_num = Iter_num)

## best results ##
best.iter <- which(results$valmaps == min(results$valmaps[1:length(results$Coef)]))
best_valperf <- results$valpredmatrix[[best.iter]]
best_evlperf <- results$evlpredmatrix[[best.iter]]
best_coef <- results$Coef[[best.iter]]
best_numcoef <- results$NumbernonzeroCoef[[best.iter]]
best_coef_name <- results$Coef.name[[best.iter]]

## Record best result
results.breast<- list("best.valperf" = best_valperf, 
                      "best.evlperf" = best_evlperf,
                      "best.coef" = best_coef, 
                      "best.numcoef" = best_numcoef,
                      "best.coef.name" = best_coef_name)

# 3. Evaluate the performance
##----------------------------

conMatrix.train <- confusionMatrix(as.factor(train_label),as.factor(results.breast$best.valperf))
conMatrix.test <- confusionMatrix(as.factor(test_label),as.factor(results.breast$best.evlperf),mode = "prec_recall")


##----------------------
## feature
##----------------------
## calculate features ##
coef.idx <- which(results.breast$best.coef$`1`[-1]!=0)
coef.value <- results.breast$best.coef$`1`[which(results.breast$best.coef$`1`!=0)][-1]
coef.name <- colnames(combined_train_data)[coef.idx]

coef <- cbind(coef.name, coef.value)
coef.mRNA <- coef[grep("mrna", coef),]
coef.miRNA <- coef[grep("mirna", coef),]
coef.cpg <- coef[grep("meth", coef),]

num_feature <- matrix(0, nrow = 1, ncol = 3)
num_feature[,1] <- length(grep("mrna", coef.name))
num_feature[,2] <- length(grep("mirna", coef.name))
num_feature[,3] <- length(grep("meth", coef.name))

perf.ConcateSPL <- list("coef.mRNA" = coef.mRNA,
                        "coef.miRNA" = coef.miRNA,
                        "coef.cpg" = coef.cpg,
                        "feature.num" = num_feature,
                        "Perf.Train" = conMatrix.train,
                        "Perf.Test" = conMatrix.test)
