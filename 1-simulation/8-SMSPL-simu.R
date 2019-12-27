##----------------------------------
##  This simulation code for multimodal self-paced learning with soft weighting scheme.
##  SMSPL
##  Code Version 1.0
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------

setwd("D:/path")
source("Function_simu_performance.R")

##---------------------
## load libraries
##---------------------

library(Matrix)
library(tseries)
library(glmnet)
library(caret)
library(ncvreg)
library(pROC)
library(ROCR)
library(ggplot2)
library(reticulate)

##--------------------------------
## 1. Generate simulation data
##--------------------------------
source("Generate_Simdata.R")
#source("Generate_Simdata_confidence.R")

##-----------------
## Step 1 : Initialization parameters
##-----------------
View_num <- 3
iter_num <- 100
gamma <- 0.66                        # adjusting parameters based on the real problem
lambda <- c(0.76, 0.76, 0.76)        # adjusting parameters
uplambda <- 0.02                     # adjusting parameters
num_sample <- 4                      # adjusting selected samaples in each iteration
num_up <- 2


valpredmatrix <- list()
evlpredmatrix <- list()
coefmatrix <- list()
nonzerocoefmatrix <- list()
coef_idx <- list()
coef_value <- list()
iter.weight <- list()                
selectedidx <- list()
sample.weight <- list()

loss <- matrix(0, nrow = length(y_train), ncol = View_num)
v_iter <- matrix(0, nrow = length(y_train), ncol = View_num)

for(iter in 1:iter_num) {
  valpredmatrix[[iter]] <- matrix(0, nrow = length(y_train), ncol = View_num)
  evlpredmatrix[[iter]] <- matrix(0, nrow = length(y_test), ncol = View_num)
  coefmatrix[[iter]] <-  list()
  nonzerocoefmatrix[[iter]] <- matrix(0, nrow = 1, ncol = View_num)
  iter.weight[[iter]] <- list()
}

val_labels <- matrix(rep(y_train,each =3),ncol = 3, byrow = T)
evl_labels <- matrix(rep(y_test,each =3),ncol = 3, byrow = T)
valmaps <- replicate(iter_num,0)
evlmaps <- replicate(iter_num,0)

##-------------------------------------
## Step 2.1: Initialization classifiers
##-------------------------------------
## using all samples to initialize classifiers
for(i in 1:View_num){
  cvfit <- cv.glmnet(x = data.train[[i]],
                   y = y_train,
                   alpha = 1,
                   family = "binomial",
                   type.measure = "class")
  valpredmatrix[[1]][,i] <- predict(cvfit,data.train[[i]],type="response",s="lambda.min")
  evlpredmatrix[[1]][,i] <- predict(cvfit,data.test[[i]],type="response",s="lambda.min")
}

## select weight with non-zero samples under the initialization
v_iter <- SMSPL.rank(dev_decval = valpredmatrix[[1]],
                    dev_labels = y_train,
                    v_iter = v_iter,
                    lambda = lambda,
                    gamma = gamma,
                    View_num = View_num,
                    num_sample = num_sample)


for(i in 1:View_num){
  selectedidx[[i]] <- which(v_iter[,i]!=0)              # selected samples index
  sample.weight[[i]] <- v_iter[which(v_iter[,i]!=0),i]  # the weight of selected samples
  iter.weight[[1]][[i]] <- sample.weight[[i]]
}
val_loss <- sum((valpredmatrix[[1]] - val_labels)^2)
evl_loss <- sum((evlpredmatrix[[1]] - evl_labels)^2)
valmaps[1] <- val_loss
evlmaps[1] <- evl_loss

##-----------------------
## Step 2.2: Optimization
##-----------------------
for(iter in 1:iter_num){
  cat("Starting the ",iter,"-th iteration.\n", sep = "")
  
  ## Training Ensemble classifier ##
  for(i in 1:View_num){
    cvfit <- cv.glmnet(x = data.train[[i]][selectedidx[[i]],],
                     y = y_train[selectedidx[[i]]],
                     weights = sample.weight[[i]],
                     alpha = 1,                           ## Lasso penalty
                     family = "binomial",
                     type.measure = "class")
    
    valpredmatrix[[iter]][,i] <- predict(cvfit, data.train[[i]], type = "response", s = "lambda.min")
    evlpredmatrix[[iter]][,i] <- predict(cvfit, data.test[[i]], type = "response", s = "lambda.min")
    coefmatrix[[iter]][[i]] <- as.vector(coef(cvfit,s = "lambda.min"))
    nonzerocoefmatrix[[iter]][[i]] <- length(which(coefmatrix[[iter]][[i]]!=0)-1)
  }
  
  ## evaluate the training and test error
  val_loss <- sum((valpredmatrix[[iter]] - val_labels)^2)
  evl_loss <- sum((evlpredmatrix[[iter]] - evl_labels)^2)
  valmaps[iter] <- val_loss
  evlmaps[iter] <- evl_loss
  
  ## if all the sample used to tranining the classifiers, then stop the iteration.
  if(length(unlist(selectedidx)) == (length(y_train)*View_num)){break}
  
  ## update lambda and valpredmatrix for next iteriation
  lambda <-  uplambda + lambda
  num_sample <- num_sample + num_up

  ## select samples
  v_iter <- SMSPL.rank(dev_decval = valpredmatrix[[iter]],
                      dev_labels = y_train,
                      v_iter = v_iter,
                      lambda = lambda,
                      gamma = gamma,
                      View_num = View_num,
                      num_sample = num_sample)

  for(i in 1:View_num){
    selectedidx[[i]] <- which(v_iter[,i]!=0)              # selected samples index
    sample.weight[[i]] <- v_iter[which(v_iter[,i]!=0),i]  # the weight of selected samples
    iter.weight[[iter]][[i]] <- sample.weight[[i]]
  }
}

##----------------------------------------------------
## Step 3: Record performance in each iteration
##----------------------------------------------------
## best results ##
best.iter <- which(valmaps == min(valmaps[1:length(which(valmaps!=0))]))
best_valperf <- valpredmatrix[[best.iter]]
best_evlperf <- evlpredmatrix[[best.iter]]
best_coef <- coefmatrix[[best.iter]]
best_numcoef <- nonzerocoefmatrix[[best.iter]]

MVSPL.result <- list("best.valperf" = best_valperf, 
                     "best.evlperf" = best_evlperf,
                     "best.coef" = best_coef, 
                     "best.numcoef" = best_numcoef)

##-------------------------------------------
## Step 4: Evaluate the best performance
##-------------------------------------------
pref.train <- evaluate.performance(pred_label = MVSPL.result$best.valperf, true_label = y_train, View_num)
pref.test <- evaluate.performance(pred_label = MVSPL.result$best.evlperf, true_label = y_test, View_num)
perf.beta <- evaluate.beta.performance(pred_beta = MVSPL.result$best.coef, true_beta = beta.matrix)


SMSPL.performance <- list("pref.train" = pref.train, "perf.test" = pref.test, "perf.beta" = perf.beta)

