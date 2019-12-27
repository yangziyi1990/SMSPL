##----------------------------------
##  This simulation code for logistic regression with Elastic net penalty embeded in the concatenation-based framework.
##  Concate_EN
##  Code Version 1.0
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

##---------------------
## load libraries
##---------------------
library(mixOmics)
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

##-----------------
# 2. Concatenation
##-----------------
combined_train_data <- do.call(cbind, data.train)
combined_test_data <- do.call(cbind, data.test)
combined_beta <-c(beta.matrix[[1]],beta.matrix[[2]],beta.matrix[[3]])

cvfit<-cv.glmnet(x = combined_train_data,
                 y = y_train,
                 alpha = 0.9,
                 family = "binomial",
                 type.measure = "class")

valprediction <- predict(cvfit, combined_train_data, type="response", s="lambda.min")
tsprediction <- predict(cvfit, combined_test_data, type="response", s="lambda.min")
coefprediction <- as.vector(coef(cvfit,s="lambda.min"))
numbernonzerocoef <- length(which(coefprediction[-1]!=0))

results <- list("valpredmatrix" = valprediction, 
                "evlpredmatrix" = tsprediction, 
                "Coef" = coefprediction, 
                "NumbernonzeroCoef" = numbernonzerocoef,
                "cvfit" = cvfit)

##--------------
# 3. Performance
##--------------
para.runs = list()
para.runs = evaluate.ConcatenationEN.performance(results, y_train, y_test, combined_beta)

