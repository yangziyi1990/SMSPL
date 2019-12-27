##----------------------------------
##  This simulation code for logistic regression with Elastic net penalty embeded in the ensemble-based framework.
##  Ensemble_EN
##  Code Version 1.0
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

##---------------------
## load libraries
##---------------------
library(amritr)
library(glmnet)
library(Matrix)
library(tseries)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)

setwd("D:/path")
source("Function_simu_performance.R")

##--------------------------------
## 1. Generate simulation data
##--------------------------------
source("Generate_Simdata.R")
#source("Generate_Simdata_confidence.R")


##--------------
# 2.Ensemble Enet
##--------------
ensemblePanels <- lapply(data.train, function(i){
  cvfit_ensem <- cv.glmnet(i, y_train, alpha=0.9,family="binomial",type.measure = "class")
})

ensembleValiPredction <- mapply(function(cvfit,x){
  valprediction <- predict(cvfit,x,type="response",s="lambda.min")
}, cvfit = ensemblePanels, x = data.train)

ensembleTestPredction <- mapply(function(cvfit,x){
  tsprediction <- predict(cvfit,x,type="response",s="lambda.min")
}, cvfit = ensemblePanels, x = data.test)

ensembleCoef <- lapply(ensemblePanels, function(i){
  coefprediction <- as.vector(coef(i,s="lambda.min"))
})

ensembleCoefNum <- lapply(ensembleCoef, function(i){
  numbercoef <- length(which(i[-1]!=0))
})

results <- list("valpredmatrix" = ensembleValiPredction, 
                "evlpredmatrix" = ensembleTestPredction, 
                "Coef" = ensembleCoef, 
                "NumbernonzeroCoef" = ensembleCoefNum)

##--------------
# 3. Performance
##--------------
para.runs = list()
para.runs = evaluate.EnsembleEN.performance(results, y_train, y_test, beta.matrix)

