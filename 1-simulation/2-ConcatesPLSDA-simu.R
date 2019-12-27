##----------------------------------
##  This simulation code for sPLSDA embeded in the concatenation-based framework.
##  Concate_sPLSDA
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
# 2. Concatenation_sPLSDA
##-----------------------
combined_train_data <- do.call(cbind, data.train)
combined_test_data <- do.call(cbind, data.test)
combined_beta <-c(beta.matrix[[1]],beta.matrix[[2]],beta.matrix[[3]])

X <- combined_train_data
Y <- as.factor(y_train)
X_test <- combined_test_data
Y_test <- as.factor(y_test)

## SPLSDA classification model ##
fit.splsda <- splsda(X,Y, K=5, eta=0.85, scale.x=FALSE)
coef <- coef(fit.splsda)
predict_train <- predict(fit.splsda, newx = X)
predict_test <- predict(fit.splsda, newx = X_test)

predict_train <- as.numeric(as.character(predict_train))
predict_test <- as.numeric(as.character(predict_test))

results <- list("valpredmatrix" = predict_train, "evlpredmatrix" = predict_test, 
                "Coef" = coef)

##--------------
# 3. Performance
##--------------
para.runs = list()
para.runs = evaluate.Concatenation.splsda.performance(results, y_train, y_test, combined_beta)
  
