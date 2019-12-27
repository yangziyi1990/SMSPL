##--------------------------------
## Step 1. Generate simulation data
##--------------------------------
J <- 3                           # number of omics
p <- c(2000,500,1500)            # number of features in each omics
sam_num <- 1500
train_num <- c(100,150,200)      # sample numbers of training
test_num <- 100                  # sample numbers of testing
para_correlation <- c(0)         # c(0,0.2,0.4,0.6,0.8)
para_noise <- c(0,0.4,0.8)
S <- 3                           # index of train number  (100,150,200)
N <- 1                           # index of noise number (0,0.4,0.8)

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

for(i in 1:J){
  ## Pool of training
  noise_rnorm <- rnorm(sam_num,mean = 0, sd = 4)
  l.train <- sample.train[[i]] %*% beta.matrix[[i]] + para_noise[N] * noise_rnorm
  prob.train <- exp(l.train)/(1+exp(l.train))
  prob.matrix.train[,i] <- prob.train
  prob.train[which(prob.train>0.5)] <- 1
  prob.train[which(prob.train<=0.5)] <- 0
  label.matrix.train[,i] <- prob.train
  ## Pool of testing
  l.test <- sample.test[[i]] %*% beta.matrix[[i]]
  prob.test <- exp(l.test)/(1+exp(l.test))
  prob.matrix.test[,i] <- prob.test
  prob.test[which(prob.test>0.5)] <- 1
  prob.test[which(prob.test<=0.5)] <- 0
  label.matrix.test[,i] <- prob.test
  
}

sam0_train.idx <- which(apply(label.matrix.train, 1, sum)==0)
sam1_train.idx <- which(apply(label.matrix.train, 1, sum)==3)
sam0_test.idx <- which(apply(label.matrix.test, 1, sum)==0)
sam1_test.idx <- which(apply(label.matrix.test, 1, sum)==3)

## select confidence samples ##
select0.train.idx <- sample(x = sam0_train.idx, size = train_num[S]/2)
select1.train.idx <- sample(x = sam1_train.idx, size = train_num[S]/2)
select0.test.idx <- sample(x = sam0_test.idx, size = test_num/2)
select1.test.idx <- sample(x = sam1_test.idx, size = test_num/2)

select.train.idx <- c(select0.train.idx,select1.train.idx)
select.test.idx <- c(select0.test.idx,select1.test.idx)

##------------------------------
## generate simulation data
##------------------------------

for(i in 1:J){
  data.train[[i]] <- sample.train[[i]][select.train.idx,]
  y_train <- c(rep(0,train_num[S]/2),rep(1,train_num[S]/2))
  data.test[[i]] <- sample.test[[i]][select.test.idx,]
  y_test <- c(rep(0,test_num/2),rep(1,test_num/2))
}
