evaluate.ensemble <- function(predictlabels,truelabels){
  
  prelable_final = matrix(0, nrow = length(truelabels), ncol = 1)
  for(i in 1:length(truelabels)){
    freq_table <- as.data.frame(table(predictlabels[i,]))
    vote_index <- which.max(freq_table$Freq)
    prelable_final[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
  }
  return(prelable_final)
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


mvselfpace.rank.multiclass <- function(dev_prob, true_label, lambda, 
                                       gamma, v_iter, View_id, View_num, 
                                       sample_select) {
  
  # initialize the loss function
  loss = matrix(0, nrow = length(true_label), ncol = View_num)
  label = matrix(1, nrow = length(true_label), ncol = 1)
  
  #calculate the loss
  for(m in 1: View_num){
    if(m != View_id){
      loss[,m] = (dev_prob[,m] - label)^2    #squared error
    }else{
      next;
    }
  }
  
  # Update v(View_num-j)
  for(m in 1:View_num){
    if(m != View_id){
      for(i in 1:length(true_label)){
        if(loss[i,m] < lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m])){
          v_iter[i,m] = 1
        }else{
          v_iter[i,m] = 0
        }
      }
    }
  }
  
  # Update vj
  loss[,View_id] = (dev_prob[,View_id] - label)^2
  for(i in 1:length(true_label)){
    if(loss[i,View_id] < lambda[View_id] + gamma * (sum(v_iter[i,])-v_iter[i,View_id])){
      v_iter[i,View_id] = 1
    }else{
      v_iter[i,View_id] = 0
    }
  }
  
  ## sort sample
  class.idx <- list()
  sample.thr <- list()
  selectedidx <- list()
  selectedsample <- list()
  V_iter = matrix(0, nrow = length(true_label), ncol = View_num)
  for(i in 1:length(unique(true_label))){
    sample.thr[[i]] <- matrix(0, nrow = 1, ncol = View_num)
  }
  
  for(i in 1:View_num){
    for(j in 1:length(unique(true_label))){
      
      class.idx <- which(true_label==j)
      
      if(length(which(v_iter[class.idx,i]==1))!=0){
        sample.thr[[j]][,i] <- sort(loss[class.idx,i][which(v_iter[class.idx,i]==1)])[min(length(class.idx), 
                                                                                          sample_select[[j]],
                                                                                          length(which(v_iter[class.idx,i]==1)))]
      }
      if(length(unique(loss[class.idx,i]))!=1){
        selectedidx[[j]] <- intersect(class.idx[which(v_iter[class.idx,i] == 1)],   ## v_iter = 1 && loss is small
                                      class.idx[which(loss[class.idx,i] <= sample.thr[[j]][,i])])
      }
    }
    selectedsample[[i]] <- unlist(selectedidx)
    V_iter[selectedsample[[i]],i] = 1
  }

  cat("The ",View_id, "-th modality select ",length(selectedsample[[i]])," samples.\n", sep = "")
  
  #return the result
  return(V_iter)
  
}

calculate.final.label <- function(pred_label, true_label){
  final_label <- matrix(0, nrow = length(true_label), ncol = 1)
  for(i in 1:length(true_label)){
    freq_table <- as.data.frame(table(pred_label[i,]))
    vote_index <- which.max(freq_table$Freq)
    final_label[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
  }
  return(final_label)
  
}

selfpaced.muliticlass <- function(X_train, Y_train, X_test, Y_test, lambda, uplambda, Iter_num) {
  
  #the list storing the result for each iteration
  valpredmatrix = list()
  valprobmatrix = list()
  evlpredmatrix = list()
  evlprobmatrix = list()
  coefmatrix = list()
  coefnummatrix = list()
  coefnamematrix = list()
  valmaps <- replicate(Iter_num,0)
  evlmaps <- replicate(Iter_num,0)
  label_train <- matrix(1, nrow = length(Y_train), ncol = 1)
  label_test <- matrix(1, nrow = length(Y_test), ncol = 1)
  
  #the starting values
  cvfit<-cv.glmnet(X_train,Y_train,alpha=1,family="multinomial",type.multinomial="grouped")
  val.pred <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
  val.prob <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
  
  v.idx = selfpace.rank.multiple(dev_decval = val.prob, 
                                 dev_labels = Y_train, 
                                 lambda = lambda)
  this.training.vidx = v.idx
  iters.vidx = list()	
  
  for(iter in 1:Iter_num) {
    if(length(this.training.vidx) == length(Y_train)){break}
    cat("Starting the ",iter,"-th iteration.\t", sep = "")
    # iters.vidx[[iter]] = this.training.vidx
    
    # glmnet (Lasso, Elastic net & L2)
    cvfit<-cv.glmnet(data.matrix(X_train[this.training.vidx,]),
                     Y_train[this.training.vidx],
                     alpha=1,
                     family="multinomial",
                     type.multinomial="grouped",
                     nfolds = 5, 
                     lambda = seq(0.06,0.11,by=0.01))
    
    valprediction <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
    valprobiction <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
    
    tsprediction <- as.numeric(predict(cvfit,X_test,type="class",s="lambda.min"))
    tsprobiction <- apply(predict(cvfit,X_test,type="response",s="lambda.min"),1,max)
    
    coefprediction <- predict(cvfit,type="coefficients",s="lambda.min")  # coef
    coef.idx <- which(coefprediction$`1`[-1]!=0)
    coef.name <- rownames(coefprediction$`1`)[coef.idx]
    coef.number <- length(coef.idx)
    
    #evaluate the training and test error
    val_loss <- sum((valprobiction - label_train)^2)
    evl_loss <- sum((tsprobiction - label_test)^2)
    
    #self-paced learning
    selectedidx = selfpace.rank.multiple(dev_decval = valprobiction, 
                                         dev_labels = Y_train, 
                                         lambda = lambda)
    this.training.vidx = selectedidx
    cat("Select ", length(selectedidx), " samples.\t", sep = "")
    
    #change the parameter accoding to the step size
    for(i in 1:length(lambda)){
      lambda[[i]] = lambda[[i]] + uplambda[[i]]
    }
    
    
    #store the prediction for this iteration
    valpredmatrix[[iter]] = valprediction
    evlpredmatrix[[iter]] = tsprediction
    coefmatrix[[iter]]= coefprediction
    coefnummatrix[[iter]] = coef.number
    coefnamematrix[[iter]] = coef.name
    valmaps[iter] <- val_loss
    evlmaps[iter] <- evl_loss
    
    cat("Finish the ",iter,"-th iteration.\n", sep = "")
  }
  
  results <- list("valpredmatrix" = valpredmatrix, 
                  "evlpredmatrix" = evlpredmatrix, 
                  "valmaps" = valmaps,
                  "evlmaps" = evlmaps,
                  "Coef" = coefmatrix, 
                  "NumbernonzeroCoef" = coefnummatrix,
                  "Coef.name" = coefnamematrix)
  return(results)
}


selfpace.rank.multiple <- function(dev_decval, dev_labels, lambda) {
  
  #calculate the loss
  label = matrix(1, nrow = length(dev_labels), ncol = 1)
  loss = (dev_decval-label)^2	
  
  class.idx <- list()
  sample.thr <- list()
  selectedposidx <- list()
  for(i in 1:length(unique(dev_labels))){
    class.idx[[i]] <- which(dev_labels==i)
    sample.thr[[i]] <- sort(loss[class.idx[[i]]])[min(length(class.idx[[i]]), lambda[[i]])]
    
    if(length(unique(loss[class.idx[[i]]]))!=1){
      selectedposidx[[i]] <- class.idx[[i]][which(loss[class.idx[[i]]] <= sample.thr[[i]])]
    }else{
      selectedposidx[[i]] <- sample(class.idx[[i]], 
                                    size = min(lambda[[i]], length(class.idx[[i]])), 
                                    replace = FALSE)
    }
  }
  
  selectedidx = unlist(selectedposidx)
  
  return(selectedidx)
}



selfpaced.EN.multiclass <- function(X_train, Y_train, lambda, uplambda, Iter_num) {
  
  #the list storing the result for each iteration
  cvfitmatrix = list()
  valmaps <- replicate(Iter_num,0)
  label_train <- matrix(1, nrow = length(Y_train), ncol = 1)
  
  #the starting values
  cvfit<-cv.glmnet(X_train,Y_train,alpha=1,family="multinomial",type.multinomial="grouped")
  valpred <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
  valprob <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
  
  this.training.vidx = selfpace.rank1.multiclass(dev_decval = valprob, 
                                                 dev_labels = Y_train, 
                                                 lambda = lambda)
  
  for(iter in 1:Iter_num){
    
    if(length(this.training.vidx) == length(Y_train)){break}
    cat("Starting the ",iter,"-th iteration.\t", sep = "")
    
    # glmnet (Lasso, Elastic net & L2)
    cvfit<-cv.glmnet(data.matrix(X_train[this.training.vidx,]),
                     Y_train[this.training.vidx],
                     alpha=1,
                     family="multinomial",
                     type.multinomial="grouped",
                     nfolds = 5, 
                     lambda = seq(0.06,0.11,by=0.01))
    
    valprediction <- as.numeric(predict(cvfit,X_train,type="class",s="lambda.min"))
    valprobiction <- apply(predict(cvfit,X_train,type="response",s="lambda.min"),1,max)
    
    #self-paced learning
    selectedidx = selfpace.rank1.multiclass(dev_decval = valprobiction, dev_labels = Y_train, lambda)
    this.training.vidx = selectedidx
    cat("Select ", length(selectedidx), " samples.\t", sep = "")
    
    #change the parameter accoding to the step size
    for(i in 1:length(lambda)){
      lambda[[i]] = lambda[[i]] + uplambda[[i]]
    }
    
    #store the prediction for this iteration
    cvfitmatrix[[iter]] <- cvfit
    
    #evaluate the training and test error
    val_loss <- sum((valprobiction - label_train)^2)
    valmaps[iter] <- val_loss
    
    cat("Finish the ",iter,"-th iteration.\n", sep = "")
  }
  
  ## best results ##
  best.iter <- which(valmaps == min(valmaps[1:length(cvfitmatrix)]))
  best.cvfit <- cvfitmatrix[[best.iter]]
  
  return(best.cvfit)
}


selfpace.rank1.multiclass <- function(dev_decval, dev_labels, lambda){
  #calculate the loss
  label = matrix(1, nrow = length(dev_labels), ncol = 1)
  loss = (dev_decval-label)^2	
  
  class.idx <- list()
  sample.thr <- list()
  selectedposidx <- list()
  for(i in 1:length(unique(dev_labels))){
    class.idx[[i]] <- which(dev_labels==i)
    sample.thr[[i]] <- sort(loss[class.idx[[i]]])[min(length(class.idx[[i]]), lambda[[i]])]
    
    if(length(unique(loss[class.idx[[i]]]))!=1){
      selectedposidx[[i]] <- class.idx[[i]][which(loss[class.idx[[i]]] <= sample.thr[[i]])]
    }else{
      selectedposidx[[i]] <- sample(class.idx[[i]], 
                                    size = min(lambda[[i]], length(class.idx[[i]])), 
                                    replace = FALSE)
    }
  }
  
  selectedidx = unlist(selectedposidx)
  
  return(selectedidx)
}


evaluate.ensemble.sPLSDA <- function(predictlabels,truelabels){
  prelable_final = matrix(0, nrow = length(truelabels), ncol = 1)
  
  Prelabel_matrix = matrix(0, nrow = length(truelabels), ncol = dim(predictlabels)[2])
  colnames(Prelabel_matrix) <- colnames(predictlabels)
  
  for(i in 1:length(truelabels)){
    for(j in 1:dim(predictlabels)[2]){
      Prelabel_matrix[i,j] <- as.integer(as.character(predictlabels[i,j]))
    }
  }
  
  for(i in 1:length(truelabels)){
    freq_table <- as.data.frame(table(Prelabel_matrix[i,]))
    vote_index <- which.max(freq_table$Freq)
    prelable_final[i] <- as.numeric(as.character(freq_table$Var1[vote_index]))
  }
  return(prelable_final)
}

SMSPL.multiclass.rank <- function(dev_prob, dev_labels, v_iter, lambda, gamma, View_num, num_sample) {
  
  # Initialize the loss function
  loss = matrix(0, nrow = length(dev_labels), ncol = View_num)
  label = matrix(1, nrow = length(dev_labels), ncol = 1)
  v_iter_new = matrix(0, nrow = length(dev_labels), ncol = View_num)
  
  # Calculate the loss
  for(m in 1:View_num){
    loss[,m] = (dev_prob[,m] - label)^2    #squared error
  }
  
  # Update the weight of samples
  for(m in 1:View_num){
    for(i in 1:length(dev_labels)){
      if(loss[i,m] < lambda[m] + gamma * (sum(v_iter[i,])-v_iter[i,m]) - gamma){
        v_iter_new[i,m] <- 1
      }else if(loss[i,m] > lambda[m] + gamma * (sum(v_iter[i,]) - v_iter[i,m])){
        v_iter_new[i,m] <- 0
      }else{
        v_iter_new[i,m] <- (lambda[m] - loss[i,m])/gamma + sum(v_iter[i,]) - v_iter[i,m]
      }
    }
  }
  
  ## sort sample
  class.idx <- list()
  sample.thr <- list()
  selectedidx <- list()
  selectedsample <- list()
  V_iter = matrix(0, nrow = length(dev_labels), ncol = View_num)
  for(i in 1:length(unique(dev_labels))){
    sample.thr[[i]] <- matrix(0, nrow = 1, ncol = View_num)
  }
  
  for(i in 1:View_num){
    for(j in 1:length(unique(dev_labels))){
      
      class.idx <- which(dev_labels==j)
      
      if(length(which(v_iter_new[class.idx,i] != 0)) != 0){
        sample.thr[[j]][,i] <- sort(loss[class.idx,i][which(v_iter_new[class.idx,i]!=0)])[min(length(class.idx), 
                                                                                          num_sample[[j]],
                                                                                          length(which(v_iter_new[class.idx,i]!=0)))]
      }
      
      if((length(unique(loss[class.idx,i]))!=1)  &&  (length(class.idx[which(loss[class.idx,i] <= sample.thr[[j]][,i])]) > 1)){
        selectedidx[[j]] <- intersect(class.idx[which(v_iter_new[class.idx,i] == 1)],   ## v_iter = 1 && loss is small
                                      class.idx[which(loss[class.idx,i] <= sample.thr[[j]][,i])])
      }
    }
    selectedsample[[i]] = unlist(selectedidx)
    V_iter[selectedsample[[i]],i] <-v_iter_new[selectedsample[[i]],i]
  }

  cat("Selected ",length(selectedsample[[i]])," samples.\n", sep = "")
  #return the result
  return(V_iter)
  
}
