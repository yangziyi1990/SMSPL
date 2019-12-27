##----------------------------------
##  This simulation code for sPLSDA embeded in the ensemble-based framework.
##  Ensemble_sPLSDA
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
library(spls)

setwd("D:/path")
source("Function_simu_performance.R")


##--------------------------------
## 1. Generate simulation data
##--------------------------------
source("Generate_Simdata.R")
#source("Generate_Simdata_confidence.R")


##-----------------------
# 2. Ensemble_sPLSDA
##-----------------------
y_train <- as.factor(y_train)
y_test <- as.factor(y_test)

## Initialize 
coefmatrix <- list()
valpredmatrix <- matrix(0, nrow = length(y_train), ncol = J)
evlpredmatrix <- matrix(0, nrow = length(y_test), ncol = J)

## Ensemble_sPLSDA
ensemblePanels <- lapply(data.train, function(i){
  cvfit_ensem <- splsda(i,y_train, K=5, eta=0.85, scale.x=FALSE)
})

ensembleValiPredction <- mapply(function(cvfit,x){
  valprediction <- predict(cvfit, newx = x)
}, cvfit = ensemblePanels, x = data.train)

ensembleTestPredction <- mapply(function(cvfit,x){
  tsprediction <- predict(cvfit, newx = x)
}, cvfit = ensemblePanels, x = data.test)

ensembleCoef <- lapply(ensemblePanels, function(i){
  coefprediction <- coef(i)
})

ensembleValiPredction <- as.data.frame(ensembleValiPredction)
ensembleTestPredction <- as.data.frame(ensembleTestPredction)
y_train <- as.numeric(as.character(y_train))
y_test <- as.numeric(as.character(y_test))

results <- list("valpredmatrix" = ensembleValiPredction, 
                "evlpredmatrix" = ensembleTestPredction, 
                "Coef" = ensembleCoef)

##--------------
# 3. Performance
##--------------
para.runs = list()
para.runs = evaluate.Ensemble.splsda.performance(results, y_train, y_test, beta.matrix)

