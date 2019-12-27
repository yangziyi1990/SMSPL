Data.Random.select <- function(snf_data, snf_group)
{
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
  
  Data.RS <- list("colon_train_data" = colon_train_data,
                  "colon_train_lable" = colon_train_lable,
                  "colon_test_data" = colon_test_data,
                  "colon_test_lable" = colon_test_lable,
                  "gbm_train_data" = gbm_train_data,
                  "gbm_train_lable" = gbm_train_lable,
                  "gbm_test_data" = gbm_test_data,
                  "gbm_test_lable" = gbm_test_lable,
                  "kidney_train_data" = kidney_train_data,
                  "kidney_train_lable" = kidney_train_lable,
                  "kidney_test_data" = kidney_test_data,
                  "kidney_test_lable" = kidney_test_lable,
                  "lung_train_data" = lung_train_data,
                  "lung_train_lable" = lung_train_lable,
                  "lung_test_data" = lung_test_data,
                  "lung_test_lable" = lung_test_lable)
  
  return(Data.RS)
}

evaluate.performance <- function(results, train_label, test_label){

  valiperformance = matrix(0, nrow=1, ncol=5)  #maps on the validation set
  testiperformance = matrix(0, nrow=1, ncol=5)   #maps on the test set
  colnames(valiperformance) = c("accuracy","sensitivity","specificity","recall","AUC") 
  colnames(testiperformance) = c("accuracy","sensitivity","specificity","recall","AUC")
  
  valiperformance[1,] = evaluate.map(results$valpredmatrix,train_label)
  testiperformance[1,] = evaluate.map(results$evlpredmatrix,test_label)
  
  return(list("valiperformance" = valiperformance, "testiperformance" = testiperformance))
}

evaluate.map <- function(predictlabels,truelabels){
  pre_label_final <- as.integer(predictlabels>0.5)
  
  # calculate AUC
  pred1 <- prediction(pre_label_final, truelabels)
  perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred1, measure = "auc")
  auc <- auc@y.values[[1]]

  TN = sum((1 - pre_label_final)*(1 - truelabels)) # A:TN
  FP = sum(pre_label_final*(1 - truelabels)) # B:FP
  FN = sum((1 - pre_label_final)*truelabels) # C:FN
  TP = sum(pre_label_final*truelabels) # D:TP
  accuracy = (TP + TN)/(TN + TP + FP + FN)
  sensitivity = TP/(TP + FN)
  specificity = TN/(TN + FP)
  recall = TP/(TP + FP)
  perf <- c(accuracy,sensitivity,specificity,recall,auc)
  return(perf)
}

evaluate.ensemble.performance <- function(results, train_label, test_label){
  
  valiperformance = matrix(0, nrow=1, ncol=5)  #maps on the validation set
  testiperformance = matrix(0, nrow=1, ncol=5)   #maps on the test set
  colnames(valiperformance) = c("accuracy","sensitivity","specificity","recall","AUC") 
  colnames(testiperformance) = c("accuracy","sensitivity","specificity","recall","AUC")
  
  valiperformance[1,] = evaluate.ensemble.map(results$valpredmatrix,train_label)
  testiperformance[1,] = evaluate.ensemble.map(results$evlpredmatrix,test_label)
  
  return(list("valiperformance" = valiperformance, "testiperformance" = testiperformance))
}

evaluate.ensemble.map <- function(predictlabels,truelabels){
  
  ## initial predict matirx
  Prelabel_matrix = matrix(0, nrow = length(truelabels), ncol = dim(predictlabels)[2])
  colnames(Prelabel_matrix) <- colnames(predictlabels)
  
  for(i in 1:length(truelabels)){
    for(j in 1:dim(predictlabels)[2]){
      Prelabel_matrix[i,j] <- as.integer(predictlabels[i,j]>0.5)
    }
  }
  
  ## voting
  prelable_final = matrix(0, nrow = length(truelabels), ncol = 1)
  for(i in 1:length(truelabels)){
    freq_table <- as.data.frame(table(Prelabel_matrix[i,]))
    vote_index <- which.max(freq_table$Freq)
    prelable_final[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
  }
  
  # calculate AUC
  pred1 <- prediction(prelable_final, truelabels)
  perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred1, measure = "auc")
  auc <- auc@y.values[[1]]
  
  TN = sum((1 - prelable_final)*(1 - truelabels)) # A:TN
  FP = sum(prelable_final*(1 - truelabels)) # B:FP
  FN = sum((1 - prelable_final)*truelabels) # C:FN
  TP = sum(prelable_final*truelabels) # D:TP
  accuracy = (TP + TN)/(TN + TP + FP + FN)
  sensitivity = TP/(TP + FN)
  specificity = TN/(TN + FP)
  recall = TP/(TP + FP)
  perf <- c(accuracy,sensitivity,specificity,recall,auc)
  return(perf)
}

evaluate.DIABLO.performance <- function(confusion.mat, true_lable, predict_label){
  tp <- confusion.mat[2,2]
  tn <- confusion.mat[1,1]
  fp <- confusion.mat[2,1]
  fn <- confusion.mat[1,2]
  
  Accuracy <- (tp + tn)/(tp + tn + fp + fn)
  Sensitivity <- tp/(tp + fn)
  Specificity <- tn/(tn + fp)
  Recall = tp/(tp + fp)
  
  predict_label <- as.factor(predict_label)
  true_lable <- as.numeric(true_lable)
  predict_label <- as.numeric(predict_label)
  pred <- prediction(predict_label, true_lable)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred, measure = "auc")
  AUC <- auc@y.values[[1]]
  
  perf <- c(Accuracy, Sensitivity, Specificity, Recall, AUC)
  return(perf)
  
}


evaluate.randforest.performance  <- function(confusion.mat, true_lable, predict_label){
  tp <- confusion.mat[2,2]
  tn <- confusion.mat[1,1]
  fp <- confusion.mat[2,1]
  fn <- confusion.mat[1,2]
  
  Accuracy <- (tp + tn)/(tp + tn + fp + fn)
  Sensitivity <- tp/(tp + fn)
  Specificity <- tn/(tn + fp)
  Recall = tp/(tp + fp)
  
  predict_label <- as.factor(predict_label)
  true_lable <- as.numeric(true_lable)
  predict_label <- as.numeric(predict_label)
  pred <- prediction(predict_label, true_lable)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred, measure = "auc")
  AUC <- auc@y.values[[1]]
  
  perf <- c(Accuracy, Sensitivity, Specificity, Recall, AUC)
  return(perf)
  
}

selfpaced.single <- function(X_train, Y_train, X_test, Y_test, lambda, uplambda, Iter_num) {
  #the list storing the result for each iteration
  valpredmatrix = list()
  evlpredmatrix = list()
  coefmatrix = list()
  nonzerocoefmatrix = list()
  iters.V_idx = list()
  valmaps <- replicate(Iter_num,0)
  evlmaps <- replicate(Iter_num,0)
  #the starting values
  cvfit<-cv.glmnet(X_train,Y_train,alpha=1,family="binomial",type.measure = "class") 
  valpred <- predict(cvfit,X_train,type="response",s="lambda.min")
  
  v.idx = selfpace2.rank(dev_decval = valpred, dev_labels = Y_train, lambda = lambda)
  this.training.vidx = v.idx
  iters.vidx = list()	
  
  for(iter in 1:Iter_num) {
    if(length(this.training.vidx) == length(Y_train)){break}
    cat("Starting the ",iter,"-th iteration.\t", sep = "")
    iters.vidx[[iter]] = this.training.vidx
    
    # glmnet (Lasso, Elastic net & L2)
    cvfit<-cv.glmnet(data.matrix(X_train[this.training.vidx,]),
                     Y_train[this.training.vidx],
                     alpha=1,
                     family="binomial",
                     type.measure = "class")
    valprediction <- predict(cvfit,X_train,type="response",s="lambda.min")
    trprediction = valprediction
    tsprediction <- predict(cvfit,X_test,type="response",s="lambda.min")
    coefprediction = as.vector(coef(cvfit,s="lambda.min")[-1])
    numbernonzerocoef = length(which(coefprediction!=0))
    
    #self-paced learning
    selectedidx = selfpace2.rank(dev_decval = trprediction, dev_labels = Y_train, lambda)
    this.training.vidx = selectedidx
    cat("Select ", length(selectedidx), " samples.\t", sep = "")
    
    #change the parameter accoding to the step size
    lambda = lambda + uplambda
    
    #store the prediction for this iteration
    coefmatrix[[iter]]= coefprediction
    nonzerocoefmatrix[[iter]] = numbernonzerocoef
    valpredmatrix[[iter]] = as.numeric(valprediction)
    evlpredmatrix[[iter]] = tsprediction
    iters.V_idx[[iter]] = selectedidx
    
    #evaluate the training and test error
    val_loss <- sum((valpredmatrix[[iter]] - Y_train)^2)
    evl_loss <- sum((evlpredmatrix[[iter]] - Y_test)^2)
    valmaps[iter] <- val_loss
    evlmaps[iter] <- evl_loss
    
    cat("Finish the ",iter,"-th iteration.\n", sep = "")
  }
  
  results <- list("valpredmatrix" = valpredmatrix, 
                  "evlpredmatrix" = evlpredmatrix, 
                  "valmaps" = valmaps,
                  "evlmaps" = evlmaps,
                  "itervidx" = iters.V_idx, 
                  "Coef"=coefmatrix, 
                  "NumbernonzeroCoef"=nonzerocoefmatrix)
  return(results)
}


selfpace2.rank <- function(dev_decval, dev_labels, lambda) {
  #calculate the loss
  loss = (dev_decval-dev_labels)^2	#squared error
  #loss = 1/(1+e^(-1*loss))			#logistic
  
  posidx = which(dev_labels==1)	#postive id mapping
  negidx = which(dev_labels==0)	#negative id mapping
  
  #calculate pos_lambda1 and neg_lambda2 according to the rank
  pos_lambda = sort(loss[posidx,1])[min(length(posidx), lambda)]
  neg_lambda = sort(loss[negidx,1])[min(length(negidx), lambda)]
  
  #it is like first sorting sampled based on the metric and then select top lambda1_rank
  if(length(unique(loss[posidx]))!=1){
    selectedposidx <- posidx[which(loss[posidx,1] <= pos_lambda)]
  }else{
    selectedposidx <- sample(posidx, size = min(lambda, length(posidx)), replace = FALSE)
  }
  if(length(unique(loss[negidx]))!=1){
    selectednegidx <- negidx[which(loss[negidx,1] <= neg_lambda)]
  }else{
    selectednegidx <- sample(negidx, size = min(lambda, length(negidx)), replace = FALSE)
  }
  
  #selectedposidx <- posidx[which(loss[posidx,1] <= pos_lambda)]
  #selectednegidx <- negidx[which(loss[negidx,1] <= neg_lambda)]

  selecedidx = c(selectedposidx, selectednegidx)
  
  return(selecedidx)
}


evaluate.SPL.performance <- function(results, y_train, y_test){
  
  valiperformance = matrix(0, nrow=1, ncol=5)  #maps on the validation set
  testiperformance = matrix(0, nrow=1, ncol=5)   #maps on the test set
  
  
  colnames(valiperformance) = c("accuracy","sensitivity","specificity","recall","AUC") 
  colnames(testiperformance) = c("accuracy","sensitivity","specificity","recall","AUC")
  
  valiperformance[1,] = evaluate.SPL.map(results$best.valperf,y_train)
  testiperformance[1,] = evaluate.SPL.map(results$best.evlperf,y_test)
  
  return(list("valiperformance" = valiperformance, "testiperformance" = testiperformance))
}

evaluate.SPL.map <- function(predictlabels,truelabels){
  pre_label_final <- as.integer(predictlabels>0.5)
  
  # calculate AUC
  pred1 <- prediction(pre_label_final, truelabels)
  perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred1, measure = "auc")
  auc <- auc@y.values[[1]]
  
  TN = sum((1 - pre_label_final)*(1 - truelabels)) # A:TN
  FP = sum(pre_label_final*(1 - truelabels)) # B:FP
  FN = sum((1 - pre_label_final)*truelabels) # C:FN
  TP = sum(pre_label_final*truelabels) # D:TP
  accuracy = (TP + TN)/(TN + TP + FP + FN)
  sensitivity = TP/(TP + FN)
  specificity = TN/(TN + FP)
  recall = TP/(TP + FP)
  perf <- c(accuracy,sensitivity,specificity,recall,auc)
  return(perf)
}


selfpaced.EN.single <- function(X_train, Y_train, lambda, uplambda, Iter_num) {
  
  #the list storing the result for each iteration
  cvfitmatrix = list()
  valmaps <- replicate(Iter_num,0)

  #the starting values
  cvfit<-cv.glmnet(X_train,Y_train,alpha=1,family="binomial",type.measure = "class") 
  valpred <- predict(cvfit,X_train,type="response",s="lambda.min")
  
  this.training.vidx = selfpace3.rank(dev_decval = valpred, dev_labels = Y_train, lambda = lambda)

  for(iter in 1:Iter_num){
    
    if(length(this.training.vidx) == length(Y_train)){break}
    cat("Starting the ",iter,"-th iteration.\t", sep = "")
    # glmnet (Lasso, Elastic net & L2)
    cvfit<-cv.glmnet(data.matrix(X_train[this.training.vidx,]),
                     Y_train[this.training.vidx],
                     alpha=1,
                     family="binomial",
                     type.measure = "class")
    valprediction <- predict(cvfit,X_train,type="response",s="lambda.min")

    #self-paced learning
    selectedidx = selfpace3.rank(dev_decval = valprediction, dev_labels = Y_train, lambda)
    this.training.vidx = selectedidx
    cat("Select ", length(selectedidx), " samples.\t", sep = "")
    
    #change the parameter accoding to the step size
    lambda = lambda + uplambda
    
    #store the prediction for this iteration
    cvfitmatrix[[iter]] <- cvfit

    #evaluate the training and test error
    val_loss <- sum((valprediction - Y_train)^2)
    valmaps[iter] <- val_loss

    cat("Finish the ",iter,"-th iteration.\n", sep = "")
  }

  ## best results ##
  best.iter <- which(valmaps == min(valmaps[1:length(cvfitmatrix)]))
  best.cvfit <- cvfitmatrix[[best.iter]]

  return(best.cvfit)
}


selfpace3.rank <- function(dev_decval, dev_labels, lambda){
  #calculate the loss
  loss = (dev_decval-dev_labels)^2	#squared error
  #loss = 1/(1+e^(-1*loss))			#logistic
  
  posidx = which(dev_labels==1)	#postive id mapping
  negidx = which(dev_labels==0)	#negative id mapping
  
  #calculate pos_lambda1 and neg_lambda2 according to the rank
  pos_lambda = sort(loss[posidx,1])[min(length(posidx), lambda)]
  neg_lambda = sort(loss[negidx,1])[min(length(negidx), lambda)]
  
  #it is like first sorting sampled based on the metric and then select top lambda1_rank
  if(length(unique(loss[posidx]))!=1){
    selectedposidx <- posidx[which(loss[posidx,1] <= pos_lambda)]
  }else{
    selectedposidx <- sample(posidx, size = min(lambda, length(posidx)), replace = FALSE)
  }
  if(length(unique(loss[negidx]))!=1){
    selectednegidx <- negidx[which(loss[negidx,1] <= neg_lambda)]
  }else{
    selectednegidx <- sample(negidx, size = min(lambda, length(negidx)), replace = FALSE)
  }
  selecedidx = c(selectedposidx, selectednegidx)
  
  return(selecedidx)
}

Concatenation_sPLSDA <- function(combined_train_data, y_train, combined_test_data){
  ## R package
  X <- combined_train_data
  Y <- as.factor(y_train)
  X_test <- combined_test_data
  
  ## SPLSDA classification model ##
  fit.splsda <- splsda(X,Y, K=3, eta=0.92, scale.x=FALSE)
  valprediction <- predict(fit.splsda, newx = X)
  tsprediction <- predict(fit.splsda, newx = X_test)
  coefprediction <- coef(fit.splsda)
  
  valprediction <- as.numeric(as.character(valprediction))
  tsprediction <- as.numeric(as.character(tsprediction))
  
  results <- list("valpredmatrix" = valprediction, "evlpredmatrix" = tsprediction,
                  "Coef" = coefprediction)
  
  return(results)
}

Ensemble_sPLSDA <- function(train_data, y_train, test_data, y_test, J){
  y_train <- as.factor(y_train)
  
  coefmatrix <- list()
  valpredmatrix <- matrix(0, nrow = length(y_train), ncol = J)
  evlpredmatrix <- matrix(0, nrow = length(y_test), ncol = J)
  
  ## Ensemble_sPLSDA
  ensemblePanels <- lapply(train_data, function(i){
    cvfit_ensem <- splsda(i,y_train, K=3, eta=0.93, scale.x=FALSE)
  })
  
  ensembleValiPredction <- mapply(function(cvfit,x){
    valprediction <- predict(cvfit, newx = x)
  }, cvfit = ensemblePanels, x = train_data)
  
  ensembleTestPredction <- mapply(function(cvfit,x){
    tsprediction <- predict(cvfit, newx = x)
  }, cvfit = ensemblePanels, x = test_data)
  
  ensembleCoef <- lapply(ensemblePanels, function(i){
    coefprediction <- coef(i)
  })
  
  ensembleValiPredction <- as.data.frame(ensembleValiPredction)
  ensembleTestPredction <- as.data.frame(ensembleTestPredction)
  
    results <- list("valpredmatrix" = ensembleValiPredction, 
                  "evlpredmatrix" = ensembleTestPredction,
                  "Coef" = ensembleCoef)
  
}

evaluate.ensemble.sPLSDA.performance <- function(results, train_label, test_label){
  
  valiperformance = matrix(0, nrow=1, ncol=5)  #maps on the validation set
  testiperformance = matrix(0, nrow=1, ncol=5)   #maps on the test set
  colnames(valiperformance) = c("accuracy","sensitivity","specificity","recall","AUC") 
  colnames(testiperformance) = c("accuracy","sensitivity","specificity","recall","AUC")
  
  valiperformance[1,] = evaluate.ensemble.sPLSDA.map(results$valpredmatrix,train_label)
  testiperformance[1,] = evaluate.ensemble.sPLSDA.map(results$evlpredmatrix,test_label)
  
  return(list("valiperformance" = valiperformance, "testiperformance" = testiperformance))
}

evaluate.ensemble.sPLSDA.map <- function(predictlabels,truelabels){
  
  ## initial predict matirx
  Prelabel_matrix = matrix(0, nrow = length(truelabels), ncol = dim(predictlabels)[2])
  colnames(Prelabel_matrix) <- colnames(predictlabels)
  
  for(i in 1:length(truelabels)){
    for(j in 1:dim(predictlabels)[2]){
      Prelabel_matrix[i,j] <- as.integer(as.character(predictlabels[i,j]))
    }
  }
  
  ## voting
  prelable_final = matrix(0, nrow = length(truelabels), ncol = 1)
  for(i in 1:length(truelabels)){
    freq_table <- as.data.frame(table(Prelabel_matrix[i,]))
    vote_index <- which.max(freq_table$Freq)
    prelable_final[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
  }
  
  # calculate AUC
  pred1 <- prediction(prelable_final, truelabels)
  perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred1, measure = "auc")
  auc <- auc@y.values[[1]]
  
  TN = sum((1 - prelable_final)*(1 - truelabels)) # A:TN
  FP = sum(prelable_final*(1 - truelabels)) # B:FP
  FN = sum((1 - prelable_final)*truelabels) # C:FN
  TP = sum(prelable_final*truelabels) # D:TP
  accuracy = (TP + TN)/(TN + TP + FP + FN)
  sensitivity = TP/(TP + FN)
  specificity = TN/(TN + FP)
  recall = TP/(TP + FP)
  perf <- c(accuracy,sensitivity,specificity,recall,auc)
  return(perf)
}

SMSPL.rank <- function(dev_decval, dev_labels, v_iter, lambda, gamma, View_num, num_sample) {
  
  # Initialize the loss function
  loss = matrix(0, nrow = length(dev_labels), ncol = View_num)
  v_iter_new = matrix(0, nrow = length(dev_labels), ncol = View_num)
  
  # Calculate the loss
  for(m in 1:View_num){
    loss[,m] = (dev_decval[,m] - dev_labels)^2    #squared error
  }
  
  # Update the weight of samples
  for(m in 1:View_num){
    for(i in 1:length(dev_labels)){
      if(loss[i,m] < lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m]) - gamma){
        v_iter_new[i,m] <- 1
      }else if(loss[i,m] > lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m])){
        v_iter_new[i,m] <- 0
      }else{
        v_iter_new[i,m] <- (lambda[m] - loss[i,m])/gamma + sum(v_iter[i,])-v_iter[i,m]
      }
    }
  }
  
  ## sort sample
  selectedposidx = list()
  selectednegidx = list()
  selecedidx = list()
  pos_lambda = matrix(0, nrow = 1, ncol = View_num)
  neg_lambda = matrix(0, nrow = 1, ncol = View_num)
  V_iter = matrix(0, nrow = length(dev_labels), ncol = View_num)
  posidx = which(dev_labels==1)	#postive id mapping
  negidx = which(dev_labels==0)	#negative id mapping
  
  for(i in 1:View_num){
    
    if(length(which(v_iter_new[posidx,i]!=0))!=0){  # a certain class has samples
      pos_lambda[i] = sort(loss[posidx,i][which(v_iter_new[posidx,i]!=0)])[min(length(posidx), 
                                                                               num_sample, 
                                                                               length(which(v_iter_new[posidx,i]!=0)))]
    }
    if(length(which(v_iter_new[negidx,i]!=0))!=0){
      neg_lambda[i] = sort(loss[negidx,i][which(v_iter_new[negidx,i]!=0)])[min(length(negidx), 
                                                                               num_sample, 
                                                                               length(which(v_iter_new[negidx,i]!=0)))]
    }
    
    if(length(unique(loss[posidx,i]))!=1 &&  (length(posidx[which(loss[posidx,i] <= pos_lambda[i])]) > 1)){  # select samples
      selectedposidx[[i]] <- intersect(posidx[which(v_iter_new[posidx,i] != 0)],   ## v_iter_new != 0 && loss is small
                                       posidx[which(loss[posidx,i] <= pos_lambda[i])])
    }else{
      selectedposidx[[i]] <- sample(posidx, size = min(num_sample, length(posidx)), replace = FALSE)
      v_iter_new[selectedposidx[[i]],i] <- 0.1   ## setting small weighting if all sample is 0
    }
    
    if(length(unique(loss[negidx,i]))!=1 && (length(negidx[which(loss[negidx,i] <= neg_lambda[i])]) > 1)){
      selectednegidx[[i]] <- intersect(negidx[which(v_iter_new[negidx,i] != 0)],
                                       negidx[which(loss[negidx,i] <= neg_lambda[i])])
    }else{
      selectednegidx[[i]] <- sample(negidx, size = min(num_sample, length(negidx)), replace = FALSE)
      v_iter_new[selectednegidx[[i]],i] <- 0.1   ## setting small weighting if all sample is 0
    }
    
    selecedidx[[i]] = c(selectedposidx[[i]], selectednegidx[[i]])
    
    V_iter[selecedidx[[i]],i] <-v_iter_new[selecedidx[[i]],i]
  }
  cat("The ",View_num, "-th modality select ",length(selecedidx[[i]])," samples.\n", sep = "")
  #return the result
  return(V_iter)
  
}

evaluate.performance.SMSPL <- function(pred_label, true_label, View_num){

  ## initial predict matirx
  Prelabel_matrix = matrix(0, nrow = length(true_label), ncol = View_num)

  for(i in 1:length(true_label)){
    for(j in 1:View_num){
      Prelabel_matrix[i,j] <- as.integer(pred_label[i,j]>0.5)
    }
  }
  
  ## voting
  prelable_final = matrix(0, nrow = length(true_label), ncol = 1)
  for(i in 1:length(true_label)){
    freq_table <- as.data.frame(table(Prelabel_matrix[i,]))
    vote_index <- which.max(freq_table$Freq)
    prelable_final[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
  }
  
  # calculate AUC
  pred1 <- prediction(prelable_final, true_label)
  perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
  auc <- performance(pred1, measure = "auc")
  auc <- auc@y.values[[1]]
  
  TN = sum((1 - prelable_final)*(1 - true_label)) # A:TN
  FP = sum(prelable_final*(1 - true_label)) # B:FP
  FN = sum((1 - prelable_final)*true_label) # C:FN
  TP = sum(prelable_final*true_label) # D:TP
  accuracy = (TP + TN)/(TN + TP + FP + FN)
  sensitivity = TP/(TP + FN)
  specificity = TN/(TN + FP)
  recall = TP/(TP + FP)
  perf <- c(accuracy,sensitivity,specificity,recall,auc)
  return(perf)
  
}
