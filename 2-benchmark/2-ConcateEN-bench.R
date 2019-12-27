##----------------------------------
##  This real data code for logistic regression with Elastic net penalty embeded in the concatenation-based framework.
##  Concate_EN
##  Code Version 1.0
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

# Load libraries
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
## Load data
##--------------
load("1-data/SNFdatasets.RDATA")

  
#----------------
## Colon dataset
#----------------
# Preparing data
data_colon <- snf_data$colon
label_colon <- snf_group$colon

# random select samples and Setting training and testing
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

label_colon[which(label_colon=="high")] <- 1
label_colon[which(label_colon=="low")] <- 0

train_colon_label <- label_colon[trainidx_colon]
test_colon_label <- label_colon[testidx_colon]

## Generate Colon data
colon_train_data = list(mrna=train_colon_mrna, mirna=train_colon_mirna, cpg=train_colon_cpg)
colon_train_lable = train_colon_label
colon_test_data = list(mrna=test_colon_mrna, mirna=test_colon_mirna, cpg=test_colon_cpg)
colon_test_lable = test_colon_label

#----------------
## GBM dataset
#----------------
# Preparing data
data_gbm <- snf_data$gbm
label_gbm <- snf_group$gbm

# random select samples and Setting training and testing
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

label_gbm[which(label_gbm=="high")] <- 1
label_gbm[which(label_gbm=="low")] <- 0

train_gbm_label <- label_gbm[trainidx_gbm]
test_gbm_label <- label_gbm[testidx_gbm]

## Generate gbm data
gbm_train_data = list(mrna=train_gbm_mrna, mirna=train_gbm_mirna, cpg=train_gbm_cpg)
gbm_train_lable = train_gbm_label
gbm_test_data = list(mrna=test_gbm_mrna, mirna=test_gbm_mirna, cpg=test_gbm_cpg)
gbm_test_lable = test_gbm_label

#----------------
## Kidney dataset
#----------------
# Preparing data
data_kidney <- snf_data$kidney
label_kidney <- snf_group$kidney

# random select samples and Setting training and testing
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

label_kidney[which(label_kidney=="high")] <- 1
label_kidney[which(label_kidney=="low")] <- 0

train_kidney_label <- label_kidney[trainidx_kidney]
test_kidney_label <- label_kidney[testidx_kidney]

## Generate kidney data
kidney_train_data = list(mrna=train_kidney_mrna, mirna=train_kidney_mirna, cpg=train_kidney_cpg)
kidney_train_lable = train_kidney_label
kidney_test_data = list(mrna=test_kidney_mrna, mirna=test_kidney_mirna, cpg=test_kidney_cpg)
kidney_test_lable = test_kidney_label


#----------------
## Lung dataset
#----------------
# Preparing data
data_lung <- snf_data$lung
label_lung <- snf_group$lung

# random select samples and Setting training and testing
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

label_lung[which(label_lung=="high")] <- 1
label_lung[which(label_lung=="low")] <- 0

train_lung_label <- label_lung[trainidx_lung]
test_lung_label <- label_lung[testidx_lung]

## Generate lung data
lung_train_data = list(mrna=train_lung_mrna, mirna=train_lung_mirna, cpg=train_lung_cpg)
lung_train_lable = train_lung_label
lung_test_data = list(mrna=test_lung_mrna, mirna=test_lung_mirna, cpg=test_lung_cpg)
lung_test_lable = test_lung_label

##--------------
# 2. Concatenation
##--------------

##----------------
## 2.1 Colon dataset
#----------------
combined_train_colon <- do.call(cbind, colon_train_data)
combined_test_colon <- do.call(cbind, colon_test_data)

cvfit<-cv.glmnet(combined_train_colon,colon_train_lable,alpha=0.9,family="binomial",type.measure = "class") 
valprediction <- predict(cvfit,combined_train_colon,type="response",s="lambda.min")
tsprediction <- predict(cvfit,combined_test_colon,type="response",s="lambda.min")
coefprediction <- as.vector(coef(cvfit,s="lambda.min"))
numbernonzerocoef <- length(which(coefprediction[-1]!=0))

results.colon <- list("valpredmatrix" = valprediction, "evlpredmatrix" = tsprediction, 
                      "Coef" = coefprediction, "NumbernonzeroCoef" = numbernonzerocoef)

##----------------
## 2.2 GBM dataset
#----------------
combined_train_gbm <- do.call(cbind, gbm_train_data)
combined_test_gbm <- do.call(cbind, gbm_test_data)

cvfit<-cv.glmnet(combined_train_gbm,gbm_train_lable,alpha=0.9,family="binomial",type.measure = "class") 
valprediction <- predict(cvfit,combined_train_gbm,type="response",s="lambda.min")
tsprediction <- predict(cvfit,combined_test_gbm,type="response",s="lambda.min")
coefprediction <- as.vector(coef(cvfit,s="lambda.min"))
numbernonzerocoef <- length(which(coefprediction[-1]!=0))

results.gbm <- list("valpredmatrix" = valprediction, "evlpredmatrix" = tsprediction, 
                    "Coef" = coefprediction, "NumbernonzeroCoef" = numbernonzerocoef)

##----------------
## 2.3 Kidney dataset
#----------------
combined_train_kidney <- do.call(cbind, kidney_train_data)
combined_test_kidney <- do.call(cbind, kidney_test_data)

cvfit<-cv.glmnet(combined_train_kidney,kidney_train_lable,alpha=0.9,family="binomial",type.measure = "class") 
valprediction <- predict(cvfit,combined_train_kidney,type="response",s="lambda.min")
tsprediction <- predict(cvfit,combined_test_kidney,type="response",s="lambda.min")
coefprediction <- as.vector(coef(cvfit,s="lambda.min"))
numbernonzerocoef <- length(which(coefprediction[-1]!=0))

results.kidney <- list("valpredmatrix" = valprediction, "evlpredmatrix" = tsprediction, 
                       "Coef" = coefprediction, "NumbernonzeroCoef" = numbernonzerocoef)

#----------------
## 2.4 Lung dataset
#----------------

combined_train_lung <- do.call(cbind, lung_train_data)
combined_test_lung <- do.call(cbind, lung_test_data)

cvfit<-cv.glmnet(combined_train_lung,lung_train_lable,alpha=0.9,family="binomial",type.measure = "class") 
valprediction <- predict(cvfit,combined_train_lung,type="response",s="lambda.min")
tsprediction <- predict(cvfit,combined_test_lung,type="response",s="lambda.min")
coefprediction <- as.vector(coef(cvfit,s="lambda.min"))
numbernonzerocoef <- length(which(coefprediction[-1]!=0))

results.lung <- list("valpredmatrix" = valprediction, "evlpredmatrix" = tsprediction, 
                     "Coef" = coefprediction, "NumbernonzeroCoef" = numbernonzerocoef)

results.benchmark <- list(results.colon = results.colon,
                          results.glm = results.gbm,
                          results.kidney = results.kidney,
                          results.lung = results.lung)

##--------------
# 3. Performance
##--------------
para.runs.colon = list()
para.runs.colon = evaluate.performance(results.colon, as.numeric(as.character(colon_train_lable)), 
                                       as.numeric(as.character(colon_test_lable)))

para.runs.gbm = list()
para.runs.gbm = evaluate.performance(results.gbm, as.numeric(as.character(gbm_train_lable)), 
                                     as.numeric(as.character(gbm_test_lable)))

para.runs.kidney = list()
para.runs.kidney = evaluate.performance(results.kidney, as.numeric(as.character(kidney_train_lable)), 
                                        as.numeric(as.character(kidney_test_lable)))

para.runs.lung = list()
para.runs.lung = evaluate.performance(results.lung, as.numeric(as.character(lung_train_lable)), 
                                      as.numeric(as.character(lung_test_lable)))

performance.benchmark <- list(perf.colon = para.runs.colon,
                              perf.gbm = para.runs.gbm,
                              perf.kidney = para.runs.kidney,
                              perf.lung = para.runs.lung)



