##----------------------------------
##  This simulation code for DIABLO.
##  DIABLO
##  Code Version 1.0
##  Created by Zi-Yi Yang
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

##---------------------
## load libraries
##---------------------
library(mixOmics)
library(ROCR)

setwd("D:/path")
source("Function_simu_performance.R")

##--------------------------------
## 1. Generate simulation data
##--------------------------------
source("Generate_Simdata.R")
#source("Generate_Simdata_confidence.R")

##--------------------------
## 2.1 Parameter choice
##--------------------------
data.train <- list(mrna = data.train[[1]],mirna = data.train[[2]],cpg = data.train[[3]])
data.test <- list(mrna = data.test[[1]],mirna = data.test[[2]],cpg = data.test[[3]])

design = matrix(0.1, ncol = length(data.train), nrow = length(data.train), 
                dimnames = list(names(data.train), names(data.train)))
diag(design) = 0


# 2.1.1 Tuning the number of components
sgccda.res = block.splsda(X = data.train, Y = y_train, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# 2.1.2 Tuning keepX
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
tune.TCGA = tune.block.splsda(X = data.train, Y = y_train, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              cpus = 2, dist = "centroids.dist")

list.keepX = tune.TCGA$choice.keepX

##------------------------------
# 3. Final model
##------------------------------
res <- block.splsda(X = data.train, Y = y_train, ncomp = ncomp,
                    keepX = list.keepX, design = design)

##-------------------------------
# 4. Prediction
##-------------------------------

predict_train = predict(res, newdata = data.train)
confusion.mat.train = get.confusion_matrix(truth = y_train,
                                           predicted = predict_train$WeightedVote$centroids.dist[,1])

predict_test = predict(res, newdata = data.test)
confusion.mat.test = get.confusion_matrix(truth = y_test, 
                                          predicted = predict_test$WeightedVote$centroids.dist[,1])

list_selectVar_mrna = selectVar(res, comp = 1, block = "mrna")
list_selectVar_mirna <- selectVar(res, comp = 1, block = "mirna")
list_selectVar_cpg <- selectVar(res, comp = 1, block = "cpg")

Var.mrna <- gsub("X","",list_selectVar_mrna$mrna$name)
Var.mirna <- gsub("X","",list_selectVar_mirna$mirna$name)
Var.cpg <- gsub("X","",list_selectVar_cpg$cpg$name)

gsub("X","",Var.mrna)
gsub("X","",Var.mirna)
gsub("X","",Var.cpg)

Coef.idx <- list("coef.mrna" = Var.mrna,
                 "coef.mirna" = Var.mirna,
                 "coef.cpg" = Var.cpg)

##-------------------------------
# 4. Performance
##-------------------------------
perf.train <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.train, 
                                          true_lable = y_train, 
                                          predict_label = predict_train$WeightedVote$centroids.dist[,1])

perf.test <- evaluate.DIABLO.performance(confusion.mat = confusion.mat.test, 
                                         true_lable = y_test, 
                                         predict_label = predict_test$WeightedVote$centroids.dist[,1])

perf.beta <- evaluate.beta.DIABLO.performance(coef = Coef.idx, beta = beta.matrix)


DIABLO.performance <- list("perf.train" = perf.train, "perf.test" = perf.test, "perf.beta" = perf.beta)
