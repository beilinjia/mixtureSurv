########################################################
##### This file is an example for simulation study #####
##### including data generation                    #####
########################################################
##### Author: Beilin Jia
source("./functions.R")

N=1000
lambda1=exp(0);lambda2=exp(1.5)
k1=1;k2=3
ini_lambda1=0.5;ini_lambda2=1;ini_k1=0.5;ini_k2=1
expLambda=0.04
true_beta=c(0.4,0.2,-0.6,-0.3,rep(0,7))
initial_beta=matrix(rep(0,11),11,1)
nrep=1000;n.test=10000;dx=10

#### generate training and testing datasets
Train1000_2g <- list()
set.seed(364)
for(i in 1:nrep){
  Train1000_2g[[i]] <- genrDatwb_2g(N=N,dx=dx,lambda1=lambda1,lambda2=lambda2,k1=k1,k2=k2,
                                    expLambda=expLambda,true_beta=true_beta)
}

set.seed(652)
Test1000_2g <- list()
for(i in 1:nrep){
  Test1000_2g[[i]] <- genrDatwb_2g(N=n.test,dx=dx,lambda1=lambda1,lambda2=lambda2,k1=k1,k2=k2,
                                   expLambda=expLambda,true_beta=true_beta)
}

####################################
#### determine number of latent groups
####################################

logLik_n1000true2g1 <- c()
BIC_n1000true2g1 <- c()
for(i in 1:1000){
  weibull_fit <- survreg(Surv(y, delta) ~ 1, data = Train1000_2g[[i]], dist="weibull")
  logLik_n1000true2g1 <- c(logLik_n1000true2g1,weibull_fit$loglik[1])
  BIC_n1000true2g1 <- c(BIC_n1000true2g1,BIC_cal(logL=weibull_fit$loglik[1],p=2,n=N))
}

logLik_n1000true2g2 <- c()
BIC_n1000true2g2 <- c()
for(i in 1:1000){
  weibull_fit <- logLikCal(dataX=train[,1:10],Y=train$y,delta=train$delta,
                           ini_surv=matrix(c(0.5,0.5,1,1),ncol=1),dx=10,numGrp=2)
  if(weibull_fit$converge==1){
    logLik_n1000true2g2 <- c(logLik_n1000true2g2, weibull_fit$logLik)
    BIC_n1000true2g2 <- c(BIC_n1000true2g2,BIC_cal(logL=weibull_fit$logLik,p=15,n=N))
  }else{
    logLik_n1000true2g2 <- c(logLik_n1000true2g2, NA)
    BIC_n1000true2g2 <- c(BIC_n1000true2g2, NA)
  }
}

temp <- list()
logLik_n1000true2g3 <- c()
BIC_n1000true2g3 <- c()
for(i in 1:1000){
  temp <- logLikCal(dataX=train[,1:10],Y=train$y,delta=train$delta,ini_surv=matrix(c(0.5,0.5,1,1,1.3,1.3),ncol=1),
                    dx=10,numGrp=3)
  if(temp$converge==1){
    logLik_n1000true2g3 <- c(logLik_n1000true2g3, temp$logLik)
    BIC_n1000true2g3 <- c(BIC_n1000true2g3,BIC_cal(logL=temp$logLik,p=28,n=N))
  }else{
    logLik_n1000true2g3 <- c(logLik_n1000true2g3, NA)
    BIC_n1000true2g3 <- c(BIC_n1000true2g3, NA)
  }
}



####################################
#### estimate beta coefficients and variable selection
####################################

set.seed(554)
beta2gwb_EMnr_n1000 <- list()
res2gwb_EMnr_n1000 <- list()
beta2gwb_al_EMnr_n1000 <- list()
cv2gwb_al_EMnr_n1000 <- list()
for (i in 1:nrep){
  train <- Train1000_2g[[i]]
  beta2gwb_EMnr_n1000[[i]] <- betaOptim2g(dataX=train[1:10],Y=train$y,delta=train$delta,initial_beta=matrix(rep(0,11),ncol=1),tau=0.0001,dx=10,
                                          ini_lambda1=ini_lambda1,ini_lambda2=ini_lambda2,ini_k1=ini_k1,ini_k2=ini_k2)
  test <- Test1000_2g[[i]]
  B.new <- test_res2g(beta=beta2gwb_EMnr_n1000[[i]]$beta, testX=test[,1:10], dx=10)
  res2gwb_EMnr_n1000[[i]] <- cbind(B.new,test$B)
  cv2gwb_al_EMnr_n1000[[i]] <- cv_glmalasso2g(dataX=train[,1:10],Y=train$y,delta=train$delta,beta_EM=beta2gwb_EMnr_n1000[[i]]$beta,nfolds=5,
                                              lambda1_EM=beta2gwb_EMnr_n1000[[i]]$lambda1,lambda2_EM=beta2gwb_EMnr_n1000[[i]]$lambda2,
                                              k1_EM=beta2gwb_EMnr_n1000[[i]]$k1,k2_EM=beta2gwb_EMnr_n1000[[i]]$k2,dx=10)
  beta2gwb_al_EMnr_n1000[[i]] <- cv2gwb_al_EMnr_n1000[[i]]$beta
}

############################
#### post selection
############################

beta2gwb_ps_aln1000 <- list()
test2gwb_ps_aln1000 <- list()
set.seed(882)
for(i in 1:nrep){
  train <- Train1000_2g[[i]]
  beta2gwb_ps_aln1000[[i]] <- betaOptim2g_postsele(dataX=train[,1:10],Y=train$y,delta=train$delta,beta_al=beta2gwb_al_EMnr_n1000[[i]],
                                                     tau=0.0001,ini_lambda1=ini_lambda1,ini_lambda2=ini_lambda2,ini_k1=ini_k1,ini_k2=ini_k2)
  test <- Test1000_2g[[i]]
  B.new <- test2g_res_postsele(beta_ps=beta2gwb_ps_aln1000[[i]]$beta,testX=test[,1:10],
                               x2_nonzero=beta2gwb_ps_aln1000[[i]]$x2_nonzero)
  test2gwb_ps_aln1000[[i]] <- cbind(B.new, test$B)
}


