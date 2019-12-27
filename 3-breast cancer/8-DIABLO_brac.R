##----------------------------------
##  Code Version 1.0
##  This breast cancer code for DIABLO.
##  DIABLO
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

## load libraries
library(mixOmics)
library(ROCR)
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

##--------------------------
## 2.1 Parameter choice
##--------------------------
design = matrix(0.1, ncol = length(train_data), nrow = length(train_data), 
                dimnames = list(names(train_data), names(train_data)))
diag(design) = 0

# 2.1.1 Tuning the number of components
sgccda.res = block.splsda(X = train_data, Y = train_label, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# 2.1.2 Tuning keepX
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune.TCGA = tune.block.splsda(X = train_data, Y = train_label, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX


##------------------------------
# 3.1 Final model
##------------------------------
brca.res <- block.splsda(X = train_data, Y = train_label, ncomp = ncomp, 
                         keepX = list.keepX, design = design)

list_selectVar_mrna = selectVar(brca.res, comp = ncomp, block = "mrna")
list_selectVar_mirna <- selectVar(brca.res, comp = ncomp, block = "mirna")
list_selectVar_meth <- selectVar(brca.res, comp = ncomp, block = "meth")

Var.DIABLO <- list("Var_mrna" = list_selectVar_mrna$mrna$name,
                   "Var_mrna" = list_selectVar_mirna$mirna$name,
                   "Vat_meth" = list_selectVar_meth$meth$name)

##-------------------------------
# 4. Prediction
##-------------------------------
# 4.1 Performance of the model
perf.brca = perf(brca.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = 'centroids.dist')

# 4.2 Prediction on dataset
predict.brca_train = predict(brca.res, newdata = train_data)
predict.brca_test = predict(brca.res, newdata = test_data)

conMatrix.train <- confusionMatrix(as.factor(train_label),
                                   as.factor(predict.brca_train$WeightedVote$centroids.dist[,2]),
                                   mode = "prec_recall")
conMatrix.test <- confusionMatrix(as.factor(test_label),
                                  as.factor(predict.brca_test$WeightedVote$centroids.dist[,2]),
                                  mode = "prec_recall")



perf.DIABLO <- list("Feature" = Var.DIABLO,
                    "Perf.Train" = conMatrix.train,
                    "Perf.Test" = conMatrix.test)


