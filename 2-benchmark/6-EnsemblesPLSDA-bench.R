##----------------------------------
##  This benchmark code for sPLSDA embeded in the ensemble-based framework.
##  Ensemble_sPLSDA
##  Code Version 1.0
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

## Load libraries
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
source("Function_performance.R")

##--------------
# 1. Load data
##--------------
load("1-data/SNFdatasets.RDATA")
data <- Data.Random.select(snf_data,snf_group)


##----------------
# 2. Ensemble_EN
##----------------
J <- 3  # number of omics
##------------------
## 2.1 Colon dataset
##------------------
results.colon <- Ensemble_sPLSDA(train_data = data$colon_train_data,
                                 y_train = data$colon_train_lable,
                                 test_data = data$colon_test_data,
                                 y_test = data$colon_test_lable,
                                 J = J)

##----------------
## 2.2 GBM dataset
##----------------
results.gbm <- Ensemble_sPLSDA(train_data = data$gbm_train_data,
                           y_train = data$gbm_train_lable,
                           test_data = data$gbm_test_data,
                           y_test = data$gbm_test_lable,
                           J = J)

##-------------------
## 2.3 Kidney dataset
##-------------------
results.kidney <- Ensemble_sPLSDA(train_data = data$kidney_train_data,
                              y_train = data$kidney_train_lable,
                              test_data = data$kidney_test_data,
                              y_test = data$kidney_test_lable,
                              J = J)

##------------------
## 2.4 Lung dataset
##------------------
results.lung <- Ensemble_sPLSDA(train_data = data$lung_train_data,
                            y_train = data$lung_train_lable,
                            test_data = data$lung_test_data,
                            y_test = data$lung_test_lable,
                            J = J)


results.benchmark <- list(results.colon = results.colon,
                          results.gbm = results.gbm,
                          results.kidney = results.kidney,
                          results.lung = results.lung)

##--------------
## 2.5 features
##--------------
## colon
num.mrna <- length(which(results.benchmark$results.colon$Coef[[1]]!=0))
num.mirna <- length(which(results.benchmark$results.colon$Coef[[2]]!=0))
num.cpg <- length(which(results.benchmark$results.colon$Coef[[3]]!=0))

feature.colon.num <- c(num.mrna,num.mirna,num.cpg)

## gbm
num.mrna <- length(which(results.benchmark$results.gbm$Coef[[1]]!=0))
num.mirna <- length(which(results.benchmark$results.gbm$Coef[[2]]!=0))
num.cpg <- length(which(results.benchmark$results.gbm$Coef[[3]]!=0))

feature.gbm.num <- c(num.mrna,num.mirna,num.cpg)

## kidney
num.mrna <- length(which(results.benchmark$results.kidney$Coef[[1]]!=0))
num.mirna <- length(which(results.benchmark$results.kidney$Coef[[2]]!=0))
num.cpg <- length(which(results.benchmark$results.kidney$Coef[[3]]!=0))

feature.kidney.num <- c(num.mrna,num.mirna,num.cpg)

## lung
num.mrna <- length(which(results.benchmark$results.lung$Coef[[1]]!=0))
num.mirna <- length(which(results.benchmark$results.lung$Coef[[2]]!=0))
num.cpg <- length(which(results.benchmark$results.lung$Coef[[3]]!=0))

feature.lung.num <- c(num.mrna,num.mirna,num.cpg)


##--------------
# 3. Performance
##--------------
para.runs.colon = list()
para.runs.colon = evaluate.ensemble.sPLSDA.performance(results.colon, 
                                                       as.numeric(as.character(data$colon_train_lable)),
                                                       as.numeric(as.character(data$colon_test_lable)))

para.runs.gbm = list()
para.runs.gbm = evaluate.ensemble.sPLSDA.performance(results.gbm, as.numeric(as.character(data$gbm_train_lable)), 
                                              as.numeric(as.character(data$gbm_test_lable)))

para.runs.kidney = list()
para.runs.kidney = evaluate.ensemble.sPLSDA.performance(results.kidney, as.numeric(as.character(data$kidney_train_lable)), 
                                                 as.numeric(as.character(data$kidney_test_lable)))

para.runs.lung = list()
para.runs.lung = evaluate.ensemble.sPLSDA.performance(results.lung, as.numeric(as.character(data$lung_train_lable)), 
                                               as.numeric(as.character(data$lung_test_lable)))

performance.benchmark <- list(perf.colon = para.runs.colon,
                              perf.gbm = para.runs.gbm,
                              perf.kidney = para.runs.kidney,
                              perf.lung = para.runs.lung)


