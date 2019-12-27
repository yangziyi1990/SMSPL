##----------------------------------
##  This benchmark code for sPLSDA embeded in the concatenation-based framework.
##  Concate_sPLSDA
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
library(spls)

setwd("D:/path")
source("Function_performance.R")


##--------------
# 1. Load data
##--------------
load("1-data/SNFdatasets.RDATA")
data <- Data.Random.select(snf_data,snf_group)

##--------------
# 2. Concatenation_sPLSDA
##--------------
##----------------
## 2.1 Colon dataset
#----------------
combined_train_colon <- do.call(cbind, data$colon_train_data)
combined_test_colon <- do.call(cbind, data$colon_test_data)

results.colon <- Concatenation_sPLSDA(combined_train_data = combined_train_colon,
                                      y_train = data$colon_train_lable,
                                      combined_test_data = combined_test_colon)

##----------------
## 2.2 GBM dataset
#----------------
combined_train_gbm <- do.call(cbind, data$gbm_train_data)
combined_test_gbm <- do.call(cbind, data$gbm_test_data)

results.gbm <- Concatenation_sPLSDA(combined_train_data = combined_train_gbm,
                                    y_train = data$gbm_train_lable,
                                    combined_test_data = combined_test_gbm)

##----------------
## 2.3 Kidney dataset
#----------------
combined_train_kidney <- do.call(cbind, data$kidney_train_data)
combined_test_kidney <- do.call(cbind, data$kidney_test_data)

results.kidney <- Concatenation_sPLSDA(combined_train_data = combined_train_kidney,
                                       y_train = data$kidney_train_lable,
                                       combined_test_data = combined_test_kidney)

#----------------
## 2.4 Lung dataset
#----------------
combined_train_lung <- do.call(cbind, data$lung_train_data)
combined_test_lung <- do.call(cbind, data$lung_test_data)

results.lung <- Concatenation_sPLSDA(combined_train_data = combined_train_lung,
                                     y_train = data$lung_train_lable,
                                     combined_test_data = combined_test_lung)


results.benchmark <- list(results.colon = results.colon,
                          results.gbm = results.gbm,
                          results.kidney = results.kidney,
                          results.lung = results.lung)

##--------------
## 2.5 features
##--------------
## colon
combined.colon <- do.call(cbind, data$colon_train_data)
coef.colon <- results.benchmark$results.colon$Coef
feature.colon.name <- colnames(combined.colon)[which(coef.colon!=0)]
num.mrna <- length(grep("mrna*", feature.colon.name))
num.mirna <- length(grep("mirna*", feature.colon.name))
num.cpg <- length(grep("cpg*", feature.colon.name))
feature.colon.num <- c(num.mrna,num.mirna,num.cpg)

## gbm
combined.gbm <- do.call(cbind, data$gbm_train_data)
coef.gbm <- results.benchmark$results.gbm$Coef
feature.gbm.name <- colnames(combined.gbm)[which(coef.gbm!=0)]
num.mrna <- length(grep("mrna*", feature.gbm.name))
num.mirna <- length(grep("mirna*", feature.gbm.name))
num.cpg <- length(grep("cpg*", feature.gbm.name))
feature.gbm.num <- c(num.mrna,num.mirna,num.cpg)

## kidney
combined.kidney <- do.call(cbind, data$kidney_train_data)
coef.kidney <- results.benchmark$results.kidney$Coef
feature.kidney.name <- colnames(combined.kidney)[which(coef.kidney!=0)]
num.mrna <- length(grep("mrna*", feature.kidney.name))
num.mirna <- length(grep("mirna*", feature.kidney.name))
num.cpg <- length(grep("cpg*", feature.kidney.name))
feature.kidney.num <- c(num.mrna,num.mirna,num.cpg)

## lung
combined.lung <- do.call(cbind, data$lung_train_data)
coef.lung <- results.benchmark$results.lung$Coef
feature.lung.name <- colnames(combined.lung)[which(coef.lung!=0)]
num.mrna <- length(grep("mrna*", feature.lung.name))
num.mirna <- length(grep("mirna*", feature.lung.name))
num.cpg <- length(grep("cpg*", feature.lung.name))
feature.lung.num <- c(num.mrna,num.mirna,num.cpg)

##--------------
# 3. Performance
##--------------
para.runs.colon = list()
para.runs.colon = evaluate.performance(results.colon, as.numeric(as.character(data$colon_train_lable)), 
                                       as.numeric(as.character(data$colon_test_lable)))

para.runs.gbm = list()
para.runs.gbm = evaluate.performance(results.gbm, as.numeric(as.character(data$gbm_train_lable)), 
                                     as.numeric(as.character(data$gbm_test_lable)))

para.runs.kidney = list()
para.runs.kidney = evaluate.performance(results.kidney, as.numeric(as.character(data$kidney_train_lable)), 
                                        as.numeric(as.character(data$kidney_test_lable)))

para.runs.lung = list()
para.runs.lung = evaluate.performance(results.lung, as.numeric(as.character(data$lung_train_lable)), 
                                      as.numeric(as.character(data$lung_test_lable)))

performance.benchmark <- list(perf.colon = para.runs.colon,
                              perf.gbm = para.runs.gbm,
                              perf.kidney = para.runs.kidney,
                              perf.lung = para.runs.lung)


