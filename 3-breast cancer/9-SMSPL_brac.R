##----------------------------------
##  Code Version 1.0
##  This benchmark code for multimodal self-paced learning with soft weighting scheme.
##  SMSPL
##  Created by Zi-Yi Yang 
##  Modified by Zi-Yi Yang on December 26, 2019
##  Concact: yangziyi091100@163.com
##----------------------------------
library(caret)
library(glmnet)

##--------------
## 1. Load data
##--------------
source("Function_performance.R")
load("/1-data/1-datatrainTestDatasetsNormalized.RDATA")

## Preparing data
data.train <- list("mrna" = mrnaTest0,
                   "mirna" = mirnaTest0,
                   "meth" = methTest0)
colnames(data.train$mrna) <- paste("mrna", colnames(data.train$mrna), sep = "_")
colnames(data.train$mirna) <- paste("mirna", colnames(data.train$mirna), sep = "_")
colnames(data.train$meth) <- paste("meth", colnames(data.train$meth), sep = "_")
train_group <- pam50Test0$Call


data.test <- list("mrna" = mrnaTrain0,
                  "mirna" = mirnaTrain0,
                  "meth" = methTrain0)
colnames(data.test$mrna) <- paste("mrna", colnames(data.test$mrna), sep = "_")
colnames(data.test$mirna) <- paste("mirna", colnames(data.test$mirna), sep = "_")
colnames(data.test$meth) <- paste("meth", colnames(data.test$meth), sep = "_")
test_group <- pam50Train0$Call

y_train <-as.numeric(train_group)
y_test <-as.numeric(test_group)


##--------------
# 2. SMSPL
##--------------
##-----------------
## Step 1 : Initialization parameters
##-----------------
View_num = 3
iter_num = 50
gamma <- 0.56                        # adjusting parameters
lambda <- c(0.66, 0.76, 0.76)        # adjusting parameters
uplambda <- 0.02                     # adjusting parameters

## setting selected sample number in each iteration
sample.select <- list()
sample.add <- list()
Times <- 50         #approximate iteration times

for(i in 1: length(unique(y_train))){
  sample.select[[i]] <- 4
  sample.add[[i]] <- ceiling(length(which(y_train==i))/Times)
}

valpredmatrix = list()
valprobmatrix = list()
evlpredmatrix = list()
evlprobmatrix = list()
coefmatrix = list()
nonzerocoefmatrix = list()
coef_idx = list()
coef_coef = list()
coef_value = list()
coef_name = list()
iter.weight = list()
selectedidx <- list()
sample.weight <- list()

loss = matrix(0, nrow = length(y_train), ncol = View_num)
v_iter = matrix(0, nrow = length(y_train), ncol = View_num)

for(iter in 1:iter_num) {
  valpredmatrix[[iter]] = matrix(0, nrow = length(y_train), ncol = View_num)
  valprobmatrix[[iter]] = matrix(0, nrow = length(y_train), ncol = View_num)
  evlpredmatrix[[iter]] = matrix(0, nrow = length(y_test), ncol = View_num)
  evlprobmatrix[[iter]] = matrix(0, nrow = length(y_test), ncol = View_num)
  coefmatrix[[iter]] =  list()
  nonzerocoefmatrix[[iter]] = matrix(0, nrow = 1, ncol = View_num)
  iter.weight[[iter]] = list()
}

val_labels <- matrix(1, nrow = length(y_train), ncol = 3)
evl_labels <- matrix(1, nrow = length(y_test), ncol = 3)
valmaps <- replicate(iter_num,0)
evlmaps <- replicate(iter_num,0)


##------------------------------------
## Step 2.1: Initialization classifier
##------------------------------------
## using all samples to initialize classifiers
for(i in 1:View_num){
  cvfit <- cv.glmnet(x = data.train[[i]],
                     y = y_train,
                     alpha = 1,
                     family = "multinomial",
                     type.multinomial="grouped")
  valpredmatrix[[1]][,i] <- as.numeric(predict(cvfit, data.train[[i]], type = "class", s = "lambda.min"))
  valprobmatrix[[1]][,i] <- apply(predict(cvfit,data.train[[i]], type = "response", s = "lambda.min"), 1, max)
  evlpredmatrix[[i]][,i] <- as.numeric(predict(cvfit, data.test[[i]], type = "class", s = "lambda.min"))
  evlprobmatrix[[1]][,i] <- apply(predict(cvfit,data.test[[i]], type = "response", s = "lambda.min"), 1, max)
}

v_iter = SMSPL.multiclass.rank(dev_prob = valprobmatrix[[1]], 
                               dev_labels = y_train, 
                               v_iter = v_iter,
                               lambda = lambda,
                               gamma = gamma,
                               View_num = View_num,
                               num_sample = sample.select)

for(i in 1:View_num){
  selectedidx[[i]] <- which(v_iter[,i]!=0)              # selected samples index
  sample.weight[[i]] <- v_iter[which(v_iter[,i]!=0),i]  # the weight of selected samples
  iter.weight[[1]][[i]] <- sample.weight[[i]]
}
val_loss <- sum((valprobmatrix[[1]] - val_labels)^2)
evl_loss <- sum((evlprobmatrix[[1]] - evl_labels)^2)
valmaps[1] <- val_loss
evlmaps[1] <- evl_loss

##--------------------------
## Step 2.2: Optimization
##--------------------------
for(iter in 1:iter_num){

  cat("Starting the ",iter, "-th iterations.\n", sep = "")
  ## Training Ensemble classifier ##
  for(i in 1:View_num){
    cvfit <- cv.glmnet(x = data.train[[i]][selectedidx[[i]],],
                       y = y_train[selectedidx[[i]]],
                       weights = sample.weight[[i]],
                       alpha = 1,
                       family = "multinomial",
                       type.multinomial = "grouped")
    
    valpredmatrix[[iter]][,i] <- as.numeric(predict(cvfit, data.train[[i]], type = "class", s = "lambda.min"))
    valprobmatrix[[iter]][,i] <- apply(predict(cvfit,data.train[[i]], type = "response", s = "lambda.min"), 1, max)
    evlpredmatrix[[iter]][,i] <- as.numeric(predict(cvfit, data.test[[i]], type = "class", s = "lambda.min"))
    evlprobmatrix[[iter]][,i] <- apply(predict(cvfit,data.test[[i]], type = "response", s = "lambda.min"), 1, max)
    coefmatrix[[iter]][[i]] <- predict(cvfit, type = "coefficients", s = "lambda.min")
    nonzerocoefmatrix[[iter]][[i]] <- length(which(coefmatrix[[iter]][[i]]$`1`!=0)-1)
  }
  
  ## evaluate the training and test error
  val_loss <- sum((valprobmatrix[[iter]] - val_labels)^2)
  evl_loss <- sum((evlprobmatrix[[iter]] - evl_labels)^2)
  valmaps[iter] <- val_loss
  evlmaps[iter] <- evl_loss
  
  ## if all the sample used to tranining the classifiers, then stop the iteration.
  if(length(unlist(selectedidx)) == (length(y_train)*View_num)){break}
  
  ## update lambda and valpredmatrix for next iteriation
  lambda <- uplambda + lambda
  for(i in 1:length(sample.select)){
    sample.select[[i]] <- sample.select[[i]] + sample.add[[i]]
  }
  
  ## select samples
  v_iter = SMSPL.multiclass.rank(dev_prob = valprobmatrix[[iter]], 
                                 dev_labels = y_train, 
                                 v_iter = v_iter,
                                 lambda = lambda,
                                 gamma = gamma,
                                 View_num = View_num,
                                 num_sample = sample.select)
  
  for(i in 1:View_num){
    selectedidx[[i]] <- which(v_iter[,i]!=0)              # selected samples index
    sample.weight[[i]] <- v_iter[which(v_iter[,i]!=0),i]  # the weight of selected samples
    iter.weight[[iter]][[i]] <- sample.weight[[i]]
  }

}

##----------------------------------------------------
# Step 3: Find the run with the best valudation map
##----------------------------------------------------
## best results ##
best.iter <- which(valmaps == min(valmaps[1:length(which(valmaps!=0))]))
best_valperf <- valpredmatrix[[best.iter]]
best_evlperf <- evlpredmatrix[[best.iter]]
best_coef <- coefmatrix[[best.iter]]
best_numcoef <- nonzerocoefmatrix[[best.iter]]

## record label
final_val_label <- calculate.final.label(pred_label = best_valperf, true_label = y_train)
final_evl_label <- calculate.final.label(pred_label = best_evlperf, true_label = y_test)

## record selected features

for(i in 1:View_num){
  coef_idx[[i]] <- which(best_coef[[i]]$`1`!=0)[-1]
  coef_coef[[i]] <- best_coef[[i]]$`1`[coef_idx[[i]]]
  coef_name[[i]] <- rownames(best_coef[[i]]$`1`)[coef_idx[[i]]]
  
}
coef.mRNA <- cbind(coef_name[[1]], coef_coef[[1]])
coef.miRNA <- cbind(coef_name[[2]], coef_coef[[2]])
coef.meth <- cbind(coef_name[[3]], coef_coef[[3]])

## calculate features ##
num_feature <- matrix(0, nrow = 1, ncol = 3)
for(i in 1:View_num){
  num_feature[,i] <- length(coef_name[[i]])
}

##-------------------------------------------
## Step 4: Evaluate the best performance
##-------------------------------------------
conMatrix.train <- confusionMatrix(as.factor(y_train),
                                   as.factor(final_val_label),
                                   mode = "prec_recall")
conMatrix.test <- confusionMatrix(as.factor(y_test),
                                  as.factor(final_evl_label),
                                  mode = "prec_recall")

## record results
perf.SMSPL <- list("coef.mRNA" = coef.mRNA,
                   "coef.miRNA" = coef.miRNA,
                   "coef.cpg" = coef.meth,
                   "feature.num" = num_feature,
                   "Perf.Train" = conMatrix.train,
                   "Perf.Test" = conMatrix.test)
