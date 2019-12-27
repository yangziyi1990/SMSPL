##----------------------------------
##  This real data code for DIABLO.
##  DIABLO
##  Code Version 1.0
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------


## load libraries
library(mixOmics)
library(ROCR)

setwd("D:/path")
source("Function_performance.R")

##--------------
# 1. Load data
##--------------
load("1-data/SNFdatasets.RDATA")


##-----------------
# colon dataset
##-----------------
##-----------------
## 1.1 Import data
##-----------------
# 1.1.1 Preparing data
data_colon <- snf_data$colon
label_colon <- as.factor(snf_group$colon)

# 1.1.2 Random select samples and Setting training and testing
randidx_colon = sample(c(1:length(label_colon)),size=length(label_colon))
splits_colon = vector(mode = "numeric", length = length(label_colon))
splits_colon_trainidx = randidx_colon[1:(0.7*length(randidx_colon))]
splits_colon_testidx = randidx_colon[(0.7*length(randidx_colon)+1):length(randidx_colon)]

splits_colon[splits_colon_trainidx] = 0
splits_colon[splits_colon_testidx] = 1
splits_colon = as.matrix(splits_colon)

trainidx_colon = which(splits_colon[,1]==0)
testidx_colon = which(splits_colon[,1]==1)

train_colon_mrna<-data_colon$mrna[trainidx_colon,]
train_colon_mirna<-data_colon$mirna[trainidx_colon,]
train_colon_cpg<-data_colon$cpg[trainidx_colon,]

test_colon_mrna<-data_colon$mrna[testidx_colon,]
test_colon_mirna<-data_colon$mirna[testidx_colon,]
test_colon_cpg<-data_colon$cpg[testidx_colon,]

train_colon_label <- label_colon[trainidx_colon]
test_colon_label <- label_colon[testidx_colon]

## 1.1.3 Generate Colon data
colon_train_data = list(mrna=train_colon_mrna, mirna=train_colon_mirna, cpg=train_colon_cpg)
colon_train_lable = train_colon_label
colon_test_data = list(mrna=test_colon_mrna, mirna=test_colon_mirna, cpg=test_colon_cpg)
colon_test_lable = test_colon_label

##--------------------------
## 2.1 Parameter choice
##--------------------------
design = matrix(0.1, ncol = length(colon_train_data), nrow = length(colon_train_data), 
                dimnames = list(names(colon_train_data), names(colon_train_data)))
diag(design) = 0


# 2.1.1 Tuning the number of components
sgccda.res = block.splsda(X = colon_train_data, Y = colon_train_lable, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# 2.1.2 Tuning keepX
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune.TCGA = tune.block.splsda(X = colon_train_data, Y = colon_train_lable, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX

##------------------------------
# 3.1 Final model
##------------------------------
colon.res <- block.splsda(X = colon_train_data, Y = colon_train_lable, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

# colon.res$design

# mrna variables selected on component 1
list_selectVar_mrna = selectVar(colon.res, comp = 2, block = "mrna")
list_selectVar_mirna <- selectVar(colon.res, comp = 2, block = "mirna")
list_selectVar_cpg <- selectVar(colon.res, comp = 2, block = "cpg")

# list_selectVar_mrna$mrna$name
# list_selectVar_mrna$mrna$value

##-------------------------------
# 4. Prediction
##-------------------------------
# 4.1 Performance of the model
perf.colon = perf(colon.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = 'centroids.dist')

perf.colon$MajorityVote.error.rate
auc.splsda = auroc(colon.res, roc.block = "miRNA", roc.comp = 2)

# 4.2 Prediction on an external test set
predict.colon_train = predict(colon.res, newdata = colon_train_data)
confusion.mat.train = get.confusion_matrix(truth = colon_train_lable,
                                           predicted = predict.colon_train$WeightedVote$centroids.dist[,2])

predict.colon_test = predict(colon.res, newdata = colon_test_data)
confusion.mat.test = get.confusion_matrix(truth = colon_test_lable, 
                                          predicted = predict.colon_test$WeightedVote$centroids.dist[,2])

# 4.3 Calculate accuracy, sensitivity, specificity and AUC
perf.colon.train <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.train, true_lable = colon_train_lable, 
                                                predict_label = predict.colon_train$WeightedVote$centroids.dist[,2])
perf.colon.test <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.test, true_lable = colon_test_lable,
                                               predict_label = predict.colon_test$WeightedVote$centroids.dist[,2])

results.colon <- list("valiperformance" = perf.colon.train, "evlperformance" = perf.colon.test,
                      "Var_mrna" = list_selectVar_mrna, "Var_mirna" = list_selectVar_mirna,
                      "Var_cpg" = list_selectVar_cpg)





##-----------------
# GBM dataset
##-----------------
##-----------------
## 1.1 Import data
##-----------------
# 1.1.1 Preparing data
data_gbm <- snf_data$gbm
label_gbm <- as.factor(snf_group$gbm)

# 1.1.2 Random select samples and Setting training and testing
randidx_gbm = sample(c(1:length(label_gbm)),size=length(label_gbm))
splits_gbm = vector(mode = "numeric", length = length(label_gbm))
splits_gbm_trainidx = randidx_gbm[1:(0.7*length(randidx_gbm))]
splits_gbm_testidx = randidx_gbm[(0.7*length(randidx_gbm)+1):length(randidx_gbm)]

splits_gbm[splits_gbm_trainidx] = 0
splits_gbm[splits_gbm_testidx] = 1
splits_gbm = as.matrix(splits_gbm)

trainidx_gbm = which(splits_gbm[,1]==0)
testidx_gbm = which(splits_gbm[,1]==1)

train_gbm_mrna<-data_gbm$mrna[trainidx_gbm,]
train_gbm_mirna<-data_gbm$mirna[trainidx_gbm,]
train_gbm_cpg<-data_gbm$cpg[trainidx_gbm,]

test_gbm_mrna<-data_gbm$mrna[testidx_gbm,]
test_gbm_mirna<-data_gbm$mirna[testidx_gbm,]
test_gbm_cpg<-data_gbm$cpg[testidx_gbm,]

train_gbm_label <- label_gbm[trainidx_gbm]
test_gbm_label <- label_gbm[testidx_gbm]

## 1.1.3 Generate gbm data
gbm_train_data = list(mrna=train_gbm_mrna, mirna=train_gbm_mirna, cpg=train_gbm_cpg)
gbm_train_lable = train_gbm_label
gbm_test_data = list(mrna=test_gbm_mrna, mirna=test_gbm_mirna, cpg=test_gbm_cpg)
gbm_test_lable = test_gbm_label

##--------------------------
## 2.1 Parameter choice
##--------------------------
design = matrix(0.1, ncol = length(gbm_train_data), nrow = length(gbm_train_data), 
                dimnames = list(names(gbm_train_data), names(gbm_train_data)))
diag(design) = 0


# 2.1.1 Tuning the number of components
sgccda.res = block.splsda(X = gbm_train_data, Y = gbm_train_lable, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# 2.1.2 Tuning keepX
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune.TCGA = tune.block.splsda(X = gbm_train_data, Y = gbm_train_lable, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX

##------------------------------
# 3.1 Final model
##------------------------------
gbm.res <- block.splsda(X = gbm_train_data, Y = gbm_train_lable, ncomp = ncomp, 
                        keepX = list.keepX, design = design)

# gbm.res$design

# mrna variables selected on component 1
list_selectVar_mrna = selectVar(gbm.res, comp = 2, block = "mrna")
list_selectVar_mirna <- selectVar(gbm.res, comp = 2, block = "mirna")
list_selectVar_cpg <- selectVar(gbm.res, comp = 2, block = "cpg")

# list_selectVar_mrna$mrna$name
# list_selectVar_mrna$mrna$value

##-------------------------------
# 4. Prediction
##-------------------------------
# 4.1 Performance of the model
perf.gbm = perf(gbm.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = 'centroids.dist')

perf.gbm$MajorityVote.error.rate
auc.splsda = auroc(gbm.res, roc.block = "miRNA", roc.comp = 2)

# 4.2 Prediction on an external test set
predict.gbm_train = predict(gbm.res, newdata = gbm_train_data)
confusion.mat.train = get.confusion_matrix(truth = gbm_train_lable,
                                           predicted = predict.gbm_train$WeightedVote$centroids.dist[,2])

predict.gbm_test = predict(gbm.res, newdata = gbm_test_data)
confusion.mat.test = get.confusion_matrix(truth = gbm_test_lable, 
                                          predicted = predict.gbm_test$WeightedVote$centroids.dist[,2])

# 4.3 Calculate accuracy, sensitivity, specificity and AUC
perf.gbm.train <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.train, true_lable = gbm_train_lable, 
                                              predict_label = predict.gbm_train$WeightedVote$centroids.dist[,2])
perf.gbm.test <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.test, true_lable = gbm_test_lable,
                                             predict_label = predict.gbm_test$WeightedVote$centroids.dist[,2])

results.gbm <- list("valiperformance" = perf.gbm.train, "evlperformance" = perf.gbm.test,
                    "Var_mrna" = list_selectVar_mrna, "Var_mirna" = list_selectVar_mirna,
                    "Var_cpg" = list_selectVar_cpg)





##-----------------
# Kidney dataset
##-----------------
##-----------------
## 1.1 Import data
##-----------------
# 1.1.1 Preparing data
data_kidney <- snf_data$kidney
label_kidney <- as.factor(snf_group$kidney)

# 1.1.2 Random select samples and Setting training and testing
randidx_kidney = sample(c(1:length(label_kidney)),size=length(label_kidney))
splits_kidney = vector(mode = "numeric", length = length(label_kidney))
splits_kidney_trainidx = randidx_kidney[1:(0.7*length(randidx_kidney))]
splits_kidney_testidx = randidx_kidney[(0.7*length(randidx_kidney)+1):length(randidx_kidney)]

splits_kidney[splits_kidney_trainidx] = 0
splits_kidney[splits_kidney_testidx] = 1
splits_kidney = as.matrix(splits_kidney)

trainidx_kidney = which(splits_kidney[,1]==0)
testidx_kidney = which(splits_kidney[,1]==1)

train_kidney_mrna<-data_kidney$mrna[trainidx_kidney,]
train_kidney_mirna<-data_kidney$mirna[trainidx_kidney,]
train_kidney_cpg<-data_kidney$cpg[trainidx_kidney,]

test_kidney_mrna<-data_kidney$mrna[testidx_kidney,]
test_kidney_mirna<-data_kidney$mirna[testidx_kidney,]
test_kidney_cpg<-data_kidney$cpg[testidx_kidney,]

train_kidney_label <- label_kidney[trainidx_kidney]
test_kidney_label <- label_kidney[testidx_kidney]

## 1.1.3 Generate kidney data
kidney_train_data = list(mrna=train_kidney_mrna, mirna=train_kidney_mirna, cpg=train_kidney_cpg)
kidney_train_lable = train_kidney_label
kidney_test_data = list(mrna=test_kidney_mrna, mirna=test_kidney_mirna, cpg=test_kidney_cpg)
kidney_test_lable = test_kidney_label

##--------------------------
## 2.1 Parameter choice
##--------------------------
design = matrix(0.1, ncol = length(kidney_train_data), nrow = length(kidney_train_data), 
                dimnames = list(names(kidney_train_data), names(kidney_train_data)))
diag(design) = 0


# 2.1.1 Tuning the number of components
sgccda.res = block.splsda(X = kidney_train_data, Y = kidney_train_lable, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# 2.1.2 Tuning keepX
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune.TCGA = tune.block.splsda(X = kidney_train_data, Y = kidney_train_lable, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX

##------------------------------
# 3.1 Final model
##------------------------------
kidney.res <- block.splsda(X = kidney_train_data, Y = kidney_train_lable, ncomp = ncomp, 
                           keepX = list.keepX, design = design)

# kidney.res$design

# mrna variables selected on component 1
list_selectVar_mrna = selectVar(kidney.res, comp = 1, block = "mrna")
list_selectVar_mirna <- selectVar(kidney.res, comp = 1, block = "mirna")
list_selectVar_cpg <- selectVar(kidney.res, comp = 1, block = "cpg")

# list_selectVar_mrna$mrna$name
# list_selectVar_mrna$mrna$value

##-------------------------------
# 4. Prediction
##-------------------------------
# 4.1 Performance of the model
perf.kidney = perf(kidney.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = 'centroids.dist')

perf.kidney$MajorityVote.error.rate
auc.splsda = auroc(kidney.res, roc.block = "miRNA", roc.comp = 1)

# 4.2 Prediction on an external test set
predict.kidney_train = predict(kidney.res, newdata = kidney_train_data)
confusion.mat.train = get.confusion_matrix(truth = kidney_train_lable,
                                           predicted = predict.kidney_train$WeightedVote$centroids.dist[,1])

predict.kidney_test = predict(kidney.res, newdata = kidney_test_data)
confusion.mat.test = get.confusion_matrix(truth = kidney_test_lable, 
                                          predicted = predict.kidney_test$WeightedVote$centroids.dist[,1])

# 4.3 Calculate accuracy, sensitivity, specificity and AUC
perf.kidney.train <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.train, true_lable = kidney_train_lable, 
                                                 predict_label = predict.kidney_train$WeightedVote$centroids.dist[,1])
perf.kidney.test <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.test, true_lable = kidney_test_lable,
                                                predict_label = predict.kidney_test$WeightedVote$centroids.dist[,1])

results.kidney <- list("valiperformance" = perf.kidney.train, "evlperformance" = perf.kidney.test,
                       "Var_mrna" = list_selectVar_mrna, "Var_mirna" = list_selectVar_mirna,
                       "Var_cpg" = list_selectVar_cpg)







##-----------------
# Lung dataset
##-----------------
##-----------------
## 1.1 Import data
##-----------------
# 1.1.1 Preparing data
data_lung <- snf_data$lung
label_lung <- as.factor(snf_group$lung)

# 1.1.2 Random select samples and Setting training and testing
randidx_lung = sample(c(1:length(label_lung)),size=length(label_lung))
splits_lung = vector(mode = "numeric", length = length(label_lung))
splits_lung_trainidx = randidx_lung[1:(0.7*length(randidx_lung))]
splits_lung_testidx = randidx_lung[(0.7*length(randidx_lung)+1):length(randidx_lung)]

splits_lung[splits_lung_trainidx] = 0
splits_lung[splits_lung_testidx] = 1
splits_lung = as.matrix(splits_lung)

trainidx_lung = which(splits_lung[,1]==0)
testidx_lung = which(splits_lung[,1]==1)

train_lung_mrna<-data_lung$mrna[trainidx_lung,]
train_lung_mirna<-data_lung$mirna[trainidx_lung,]
train_lung_cpg<-data_lung$cpg[trainidx_lung,]

test_lung_mrna<-data_lung$mrna[testidx_lung,]
test_lung_mirna<-data_lung$mirna[testidx_lung,]
test_lung_cpg<-data_lung$cpg[testidx_lung,]

train_lung_label <- label_lung[trainidx_lung]
test_lung_label <- label_lung[testidx_lung]

## 1.1.3 Generate lung data
lung_train_data = list(mrna=train_lung_mrna, mirna=train_lung_mirna, cpg=train_lung_cpg)
lung_train_lable = train_lung_label
lung_test_data = list(mrna=test_lung_mrna, mirna=test_lung_mirna, cpg=test_lung_cpg)
lung_test_lable = test_lung_label

##--------------------------
## 2.1 Parameter choice
##--------------------------
design = matrix(0.1, ncol = length(lung_train_data), nrow = length(lung_train_data), 
                dimnames = list(names(lung_train_data), names(lung_train_data)))
diag(design) = 0


# 2.1.1 Tuning the number of components
sgccda.res = block.splsda(X = lung_train_data, Y = lung_train_lable, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# 2.1.2 Tuning keepX
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune.TCGA = tune.block.splsda(X = lung_train_data, Y = lung_train_lable, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX

##------------------------------
# 3.1 Final model
##------------------------------
lung.res <- block.splsda(X = lung_train_data, Y = lung_train_lable, ncomp = ncomp, 
                         keepX = list.keepX, design = design)

# lung.res$design

# mrna variables selected on component 1
list_selectVar_mrna = selectVar(lung.res, comp = 1, block = "mrna")
list_selectVar_mirna <- selectVar(lung.res, comp = 1, block = "mirna")
list_selectVar_cpg <- selectVar(lung.res, comp = 1, block = "cpg")

# list_selectVar_mrna$mrna$name
# list_selectVar_mrna$mrna$value

##-------------------------------
# 4. Prediction
##-------------------------------
# 4.1 Performance of the model
perf.lung = perf(lung.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = 'centroids.dist')

perf.lung$MajorityVote.error.rate
auc.splsda = auroc(lung.res, roc.block = "miRNA", roc.comp = 1)

# 4.2 Prediction on an external test set
predict.lung_train = predict(lung.res, newdata = lung_train_data)
confusion.mat.train = get.confusion_matrix(truth = lung_train_lable,
                                           predicted = predict.lung_train$WeightedVote$centroids.dist[,1])

predict.lung_test = predict(lung.res, newdata = lung_test_data)
confusion.mat.test = get.confusion_matrix(truth = lung_test_lable, 
                                          predicted = predict.lung_test$WeightedVote$centroids.dist[,1])

# 4.3 Calculate accuracy, sensitivity, specificity and AUC
perf.lung.train <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.train, true_lable = lung_train_lable, 
                                               predict_label = predict.lung_train$WeightedVote$centroids.dist[,1])
perf.lung.test <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.test, true_lable = lung_test_lable,
                                              predict_label = predict.lung_test$WeightedVote$centroids.dist[,1])

results.lung <- list("valiperformance" = perf.lung.train, "evlperformance" = perf.lung.test,
                     "Var_mrna" = list_selectVar_mrna, "Var_mirna" = list_selectVar_mirna,
                     "Var_cpg" = list_selectVar_cpg)



###-----------------------
# 5. Benchmark Performance
###-----------------------
performance.benchmark <- list(perf.colon = results.colon,
                              perf.gbm = results.gbm,
                              perf.kidney = results.kidney,
                              perf.lung = results.lung)
