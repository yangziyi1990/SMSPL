##--------------------------------
## Step 1. Generate simulation data
##--------------------------------
J <- 3 # number of omics
p <- c(2000,500,1500)  # number of features in each omics
sam_num <- 1500
train_num <- c(100,150,200) # sample numbers of training
test_num <- 100  # sample numbers of testing
para_correlation <- c(0)  # c(0,0.2,0.4,0.6,0.8)
para_noise <- c(0,0.4,0.8)
S <- 1  # index of train number  (100,150,200)
N <- 2  # index of noise number (0,0.4,0.8)

## Initalize ##
beta.matrix <- list()
data.train <- list()
data.test <- list()
sample.train <- list()
sample.test <- list()
y_train <- NULL
y_test <- NULL
label.matrix.train <- matrix(0, nrow = sam_num, ncol = J)
label.matrix.test <- matrix(0, nrow = sam_num, ncol = J)
prob.matrix.train <- matrix(0, nrow = sam_num, ncol = J)
prob.matrix.test <- matrix(0, nrow = sam_num, ncol = J)

for(i in 1:J){
  data.train[[i]] <- list()
  data.test[[i]] <- list()
  sample.train[[i]] <- list()
  sample.test[[i]] <- list()
}

## Inital beta matrix ##
beta1 <- c(1.5,-1.4,2,-1.6,2.3,2.1,-1.8,2.3,1.9,-1.6,1.5,-1.4,2,-1.6,2.3,2.1,-1.8,2.3,1.9,-1.6,rep(0,1980))
beta2 <- c(rep(0,5),1.4,-1.1,-1.8,1.5,2.1,rep(0,500-10))
beta3 <- c(rep(0,10),1.1,2,2.3,-1.5,-1.9,1.7,-1.2,-1.6,1.1,2,2.3,-1.5,-1.9,1.7,-1.2,rep(0,1500-25))

# beta1 <- c(1.5,-1.4,2,-1.6,2.3,2.1,-1.8,2.3,1.3,-1.6,rep(0,1990))
# beta2 <- c(rep(0,5),1.4,-1.1,-1.8,rep(0,500-8))
# beta3 <- c(rep(0,10),1.1,2,2.3,-1.5,-1.5,1.7,rep(0,1500-16))

# beta1 <- c(1.5,-1.4,2,-1.6,2.3,2.1,-1.8,2.3,1.3,-1.6,2.1,-1.8,2.3,1.3,-1.6,-1.2, rep(0,1984))
# beta2 <- c(rep(0,5),1.4,-1.1,-1.8,1.5, rep(0,500-9))
# beta3 <- c(rep(0,10),1.1,2,2.3,-1.5,-1.5,1.7,2,2.3,-1.5,1.6, rep(0,1500-20))

beta.matrix <- list("beta1" = beta1,
                    "beta2" = beta2,
                    "beta3" = beta3)

## Inital pool of X ##
for(i in 1:J){
  sample.train[[i]] <- matrix(rnorm(sam_num * p[i], mean = 0, sd = 1),
                              sam_num, p[i])
  sample.test[[i]] <- matrix(rnorm(sam_num * p[i], mean = 0, sd = 1),
                             sam_num, p[i])
}

## Select sample with the sample label in all omics ##
prob.train1.idx <- matrix(0, nrow = sam_num, ncol = J)
prob.train0.idx <- matrix(0, nrow = sam_num, ncol = J)
prob.test1.idx <- matrix(0, nrow = sam_num, ncol = J)
prob.test0.idx <- matrix(0, nrow = sam_num, ncol = J)

for(i in 1:J){
  ## Pool of training
  noise_rnorm <- rnorm(sam_num,mean = 0, sd = 4)
  l.train <- sample.train[[i]] %*% beta.matrix[[i]] + para_noise[N] * noise_rnorm
  prob.train <- exp(l.train)/(1+exp(l.train))
  
  r1.train <- runif(length(l.train), min = 0, max = 0)
  r1.train[which(prob.train>0.5)] <- 1
  prob.train1.idx[,i] <- r1.train
  
  r0.train <- runif(length(l.train), min = 0, max = 0)
  r0.train[which(prob.train<=0.5)] <- 1
  prob.train0.idx[,i] <- r0.train
  
  prob.matrix.train[,i] <- prob.train
  prob.train[which(prob.train>0.5)] <- 1
  prob.train[which(prob.train<=0.5)] <- 0
  label.matrix.train[,i] <- prob.train
  
  
  ## Pool of testing
  l.test <- sample.test[[i]] %*% beta.matrix[[i]]
  prob.test <- exp(l.test)/(1+exp(l.test))
  
  r1.test <- runif(length(l.test), min = 0, max = 0)
  r1.test[which(prob.test>0.5)] <- 1
  prob.test1.idx[,i] <- r1.test
  
  r0.test <- runif(length(l.test), min = 0, max = 0)
  r0.test[which(prob.test<=0.5)] <- 1
  prob.test0.idx[,i] <- r0.test
  
  prob.matrix.test[,i] <- prob.test
  prob.test[which(prob.test>0.5)] <- 1
  prob.test[which(prob.test<=0.5)] <- 0
  label.matrix.test[,i] <- prob.test
  
}
##-----------------------------------------
## setting the confidence of the samples
##-----------------------------------------
## Training dataset
train1.rank.list <- list()
train0.rank.list <- list()
high1.idx <- list()
medium1.idx <- list()
low1.idx <- list()
high0.idx <- list()
medium0.idx <- list()
low0.idx <- list()


for(i in 1:J){
  train1.rank.list[[i]] <- order(prob.matrix.train[,i][which(prob.train1.idx[,i]==1)],decreasing=TRUE)
  train0.rank.list[[i]] <- order(prob.matrix.train[,i][which(prob.train0.idx[,i]==1)],decreasing=FALSE)
  
  len.train1 <- length(train1.rank.list[[i]])
  high1.idx[[i]] <- train1.rank.list[[i]][1:round(0.3*len.train1)]
  medium1.idx[[i]] <- train1.rank.list[[i]][(round(0.3*len.train1)+1) : round(0.7*len.train1)]
  low1.idx[[i]] <- train1.rank.list[[i]][(round(0.7*len.train1)+1) : round(len.train1)]
  
  len.train0 <- length(train0.rank.list[[i]])
  high0.idx[[i]] <- train0.rank.list[[i]][1:round(0.3*len.train0)]
  medium0.idx[[i]] <- train0.rank.list[[i]][(round(0.3*len.train0)+1) : round(0.7*len.train0)]
  low0.idx[[i]] <- train0.rank.list[[i]][(round(0.7*len.train0)+1) : round(len.train0)]
  
  
  select1.train.idx.high <- which(prob.train1.idx[,i]==1)[high1.idx[[i]][1:round(0.15*train_num[S])]]
  select1.train.idx.medium <- which(prob.train1.idx[,i]==1)[medium1.idx[[i]][1:round(0.2*train_num[S])]]
  select1.train.idx.low <- which(prob.train1.idx[,i]==1)[low1.idx[[i]][1:round(0.15*train_num[S])]]
  
  select0.train.idx.high <- which(prob.train0.idx[,i]==1)[high0.idx[[i]][1:round(0.15*train_num[S])]]
  select0.train.idx.medium <- which(prob.train0.idx[,i]==1)[medium0.idx[[i]][1:round(0.2*train_num[S])]]
  select0.train.idx.low <- which(prob.train0.idx[,i]==1)[low0.idx[[i]][1:round(0.15*train_num[S])]]
  
  select.train.idx.high <- c(select0.train.idx.high, select1.train.idx.high)
  select.train.idx.medium <- c(select0.train.idx.medium, select1.train.idx.medium)
  select.train.idx.low <- c(select0.train.idx.low, select1.train.idx.low)
  
  confident.train <- c(rep("high",length(select.train.idx.high)), 
                       rep("medium",length(select.train.idx.medium)),
                       rep("low",length(select.train.idx.low)))
  
}


## Test dataset
test1.rank.list <- list()
test0.rank.list <- list()
high1.idx <- list()
medium1.idx <- list()
low1.idx <- list()
high0.idx <- list()
medium0.idx <- list()
low0.idx <- list()


for(i in 1:J){
  test1.rank.list[[i]] <- order(prob.matrix.test[,i][which(prob.test1.idx[,i]==1)],decreasing=TRUE)
  test0.rank.list[[i]] <- order(prob.matrix.test[,i][which(prob.test0.idx[,i]==1)],decreasing=FALSE)
  
  len.test1 <- length(test1.rank.list[[i]])
  high1.idx[[i]] <- test1.rank.list[[i]][1:round(0.3*len.test1)]
  medium1.idx[[i]] <- test1.rank.list[[i]][(round(0.3*len.test1)+1) : round(0.7*len.test1)]
  low1.idx[[i]] <- test1.rank.list[[i]][(round(0.7*len.test1)+1) : round(len.test1)]
  
  len.test0 <- length(test0.rank.list[[i]])
  high0.idx[[i]] <- test0.rank.list[[i]][1:round(0.3*len.test0)]
  medium0.idx[[i]] <- test0.rank.list[[i]][(round(0.3*len.test0)+1) : round(0.7*len.test0)]
  low0.idx[[i]] <- test0.rank.list[[i]][(round(0.7*len.test0)+1) : round(len.test0)]
  
  
  select1.test.idx.high <- which(prob.test1.idx[,i]==1)[high1.idx[[i]][1:round(0.15*test_num)]]
  select1.test.idx.medium <- which(prob.test1.idx[,i]==1)[medium1.idx[[i]][1:round(0.2*test_num)]]
  select1.test.idx.low <- which(prob.test1.idx[,i]==1)[low1.idx[[i]][1:round(0.15*test_num)]]
  
  select0.test.idx.high <- which(prob.test0.idx[,i]==1)[high0.idx[[i]][1:round(0.15*test_num)]]
  select0.test.idx.medium <- which(prob.test0.idx[,i]==1)[medium0.idx[[i]][1:round(0.2*test_num)]]
  select0.test.idx.low <- which(prob.test0.idx[,i]==1)[low0.idx[[i]][1:round(0.15*test_num)]]
  
  select.test.idx.high <- c(select0.test.idx.high, select1.test.idx.high)
  select.test.idx.medium <- c(select0.test.idx.medium, select1.test.idx.medium)
  select.test.idx.low <- c(select0.test.idx.low, select1.test.idx.low)
  
  confident.test <- c(rep("high",length(select.test.idx.high)), 
                       rep("medium",length(select.test.idx.medium)),
                       rep("low",length(select.test.idx.low)))
  
}

##------------------------------
## generate simulation data
##------------------------------
for(i in 1:J){
  data.train[[i]] <- sample.train[[i]][c(select.train.idx.high, select.train.idx.medium, select.train.idx.low),]
  y_train <- c(rep(0,length(select.train.idx.high)/2),
               rep(1,length(select.train.idx.high)/2),
               rep(0,length(select.train.idx.medium)/2),
               rep(1,length(select.train.idx.medium)/2),
               rep(0,length(select.train.idx.low)/2),
               rep(1,length(select.train.idx.low)/2))

  data.test[[i]] <- sample.test[[i]][c(select.test.idx.high, select.test.idx.medium, select.test.idx.low),]
  y_test <- c(rep(0,length(select.test.idx.high)/2),
               rep(1,length(select.test.idx.high)/2),
               rep(0,length(select.test.idx.medium)/2),
               rep(1,length(select.test.idx.medium)/2),
               rep(0,length(select.test.idx.low)/2),
               rep(1,length(select.test.idx.low)/2))
}

