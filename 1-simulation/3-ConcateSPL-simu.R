##----------------------------------
##  This simulation code for self-paced learning with L1 penalty embeded in the concatenation-based framework.
##  Concate_SPL
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

# tuning the parameters
sample.select = 4
sample.add = 2
Iter_num = 100  

results <- selfpaced.sim.single(X_train = combined_train_data, 
                                Y_train = y_train,
                                X_test = combined_test_data, 
                                Y_test = y_test,
                                lambda = sample.select,
                                uplambda = sample.add,
                                Iter_num = Iter_num)

## best results ##
best.iter <- which(results$valmaps == min(results$valmaps[1:length(results$itervidx)]))
best_valperf <- results$valpredmatrix[[best.iter]]
best_evlperf <- results$evlpredmatrix[[best.iter]]
best_coef <- results$Coef[[best.iter]]
best_numcoef <- results$NumbernonzeroCoef[[best.iter]]

## Record best result
results.best <- list("valpredmatrix" = best_valperf, 
                     "evlpredmatrix" = best_evlperf,
                     "Coef" = best_coef, 
                     "NumbernonzeroCoef" = best_numcoef)

##--------------
# 3. Performance
##--------------
para.runs = list()
para.runs = evaluate.ConcatenationSPL.performance(results.best, y_train, y_test, combined_beta)

