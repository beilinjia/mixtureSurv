##### This file contains functions for simulation study in "Inferring Latent Heterogeneity Using High-Dimensional Feature
##### Variables Supervised by Survival Outcome"
##### Author: Beilin Jia

###########################################
#### function to calculate log-likelihood, for BIC calculation, for different number of latent subgroups assumed in the data
###########################################
# @param: ini_surv: initial value of Weibull distribution, e.g.,ini_surv=(k1, lmabda1, k2, lambda2), where ki is the shape parameter for ith latent subgroup, and lambdai is the scale parameter for ith latent subgroup
# dx: number of covariates in the data
# numGrp: number of latent subgroups assumed in the data
logLikCal <- function(dataX,Y,delta,ini_surv,dx,numGrp){
  initial_beta <- matrix(0,ncol=1,nrow=(numGrp-1)*(dx+1))
  
  n <- nrow(dataX)
  X <- as.matrix(cbind(rep(1,n), dataX))
  
  diffv <- c()
  diff <- 100
  iter <- 0
  beta.curr <- initial_beta
  theta.curr <- ini_surv
  param.curr <- matrix(c(theta.curr,beta.curr),(dx+1)*(numGrp-1)+2*numGrp,1)
  k1.curr <- ini_surv[1,];lambda1.curr <- ini_surv[2,]
  k2.curr <- ini_surv[3,];lambda2.curr <- ini_surv[4,]
  if(numGrp==3){ k3.curr <- ini_surv[5,]
    lambda3.curr <- ini_surv[6,]}
  if(numGrp==4){ 
    k3.curr <- ini_surv[5,]
    lambda3.curr <- ini_surv[6,]
    k4.curr <- ini_surv[7,]
    lambda4.curr <- ini_surv[8,]}
  logLik.curr <- logLik_iter(dataX=dataX,Y=Y,delta=delta,surv=matrix(param.curr[1:(2*numGrp),],ncol=1),beta=matrix(param.curr[-(1:(2*numGrp)),],ncol=1),numGrp=numGrp)
  logLikv <- logLik.curr
  
  while(diff>1 & iter < 1000){
    E1 <- matrix(1,nrow=n,ncol=1)
    E2 <- exp(X%*%beta.curr[1:(dx+1),])
    if(numGrp==3){E3 <- exp(X%*%beta.curr[(dx+2):(2*dx+2),])}
    if(numGrp==4){
      E3 <- exp(X%*%beta.curr[(dx+2):(2*dx+2),])
      E4 <- exp(X%*%beta.curr[(2*dx+3):(3*(dx+1)),])}
    if(numGrp==3){
      Denom <- E1 + E2 + E3
    }else if(numGrp==4){
      Denom <- E1 + E2 + E3 + E4
    }else{Denom <- E1 + E2}
    
    S1 <- exp(-(Y/lambda1.curr)^k1.curr)
    S2 <- exp(-(Y/lambda2.curr)^k2.curr)
    f1 <- k1.curr*lambda1.curr^(-k1.curr)*Y^(k1.curr-1)*exp(-(Y/lambda1.curr)^k1.curr)
    f2 <- k2.curr*lambda2.curr^(-k2.curr)*Y^(k2.curr-1)*exp(-(Y/lambda2.curr)^k2.curr)
    Q1 <- (delta*f1+(1-delta)*S1)*E1
    Q2 <- (delta*f2+(1-delta)*S2)*E2
    if(numGrp==3){
      S3 <- exp(-(Y/lambda3.curr)^k3.curr)
      f3 <- k3.curr*lambda3.curr^(-k3.curr)*Y^(k3.curr-1)*exp(-(Y/lambda3.curr)^k3.curr)
      Q3 <- (delta*f3+(1-delta)*S3)*E3
    }
    if(numGrp==4){
      S3 <- exp(-(Y/lambda3.curr)^k3.curr)
      f3 <- k3.curr*lambda3.curr^(-k3.curr)*Y^(k3.curr-1)*exp(-(Y/lambda3.curr)^k3.curr)
      Q3 <- (delta*f3+(1-delta)*S3)*E3
      S4 <- exp(-(Y/lambda4.curr)^k4.curr)
      f4 <- k4.curr*lambda4.curr^(-k4.curr)*Y^(k4.curr-1)*exp(-(Y/lambda4.curr)^k4.curr)
      Q4 <- (delta*f4+(1-delta)*S4)*E4
    }
    if(numGrp==3){
      SQ <- Q1+Q2+Q3
      Q1 <- Q1/SQ
      Q2 <- Q2/SQ 
      Q3 <- Q3/SQ 
    }else if(numGrp==4){
      SQ <- Q1+Q2+Q3+Q4
      Q1 <- Q1/SQ
      Q2 <- Q2/SQ 
      Q3 <- Q3/SQ
      Q4 <- Q4/SQ
    }else{SQ <- Q1+Q2
      Q1 <- Q1/SQ
      Q2 <- Q2/SQ}
    
    parDerv1_1beta <- t(X) %*% (Q2-E2/Denom)
    Derv1beta <- parDerv1_1beta
    
    parDerv1_k1 <- sum(Q1*(delta*(1/k1.curr + log(Y/lambda1.curr)) - (Y/lambda1.curr)^k1.curr * log(Y/lambda1.curr)))
    parDerv1_lambda1 <- sum(Q1*(delta*(-k1.curr/lambda1.curr) + k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)))
    parDerv1_k2 <- sum(Q2*(delta*(1/k2.curr + log(Y/lambda2.curr)) - (Y/lambda2.curr)^k2.curr * log(Y/lambda2.curr)))
    parDerv1_lambda2 <- sum(Q2*(delta*(-k2.curr/lambda2.curr) + k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    parDerv1_surv <- matrix(c(parDerv1_k1,parDerv1_lambda1,parDerv1_k2,parDerv1_lambda2),4,1)
    
    parDerv2_beta <- - t(X) %*% diag(as.vector(E2*(Denom-E2)/Denom^2),nrow=n,ncol=n) %*% X
    Derv2beta <- parDerv2_beta
    
    parDerv2_k1k1 <- sum(Q1*(delta*(-1/k1.curr^2) - (Y/lambda1.curr)^k1.curr*(log(Y/lambda1.curr))^2))
    parDerv2_k1lambda1 <- sum(Q1*(delta*(-1/lambda1.curr) + Y^k1.curr*lambda1.curr^(-k1.curr-1)*(k1.curr*log(Y/lambda1.curr)+1)))
    parDerv2_lambda1lambda1 <- sum(Q1*(delta*(k1.curr/lambda1.curr^2) - k1.curr*(k1.curr+1)*Y^k1.curr*lambda1.curr^(-k1.curr-2)))
    parDerv2_k2k2 <- sum(Q2*(delta*(-1/k2.curr^2) - (Y/lambda2.curr)^k2.curr*(log(Y/lambda2.curr))^2))
    parDerv2_k2lambda2 <- sum(Q2*(delta*(-1/lambda2.curr)+Y^k2.curr*lambda2.curr^(-k2.curr-1)*(k2.curr*log(Y/lambda2.curr)+1)))
    parDerv2_lambda2lambda2 <- sum(Q2*(delta*(k2.curr/lambda2.curr^2) - k2.curr*(k2.curr+1)*Y^k2.curr*lambda2.curr^(-k2.curr-2)))
    parDerv2_surv <- matrix(c(parDerv2_k1k1,parDerv2_k1lambda1,0,0,
                              parDerv2_k1lambda1,parDerv2_lambda1lambda1,0,0,
                              0,0,parDerv2_k2k2,parDerv2_k2lambda2,
                              0,0,parDerv2_k2lambda2,parDerv2_lambda2lambda2),4,4)
    if(numGrp==3){
      parDerv1_2beta <- t(X) %*% (Q3-E3/Denom)
      Derv1beta <- rbind(parDerv1_1beta,parDerv1_2beta)
      
      parDerv1_k3 <- sum(Q3*(delta*(1/k3.curr + log(Y/lambda3.curr)) - (Y/lambda3.curr)^k3.curr * log(Y/lambda3.curr)))
      parDerv1_lambda3 <- sum(Q3*(delta*(-k3.curr/lambda3.curr) + k3.curr*Y^k3.curr*lambda3.curr^(-k3.curr-1)))
      parDerv1_surv <- matrix(c(parDerv1_k1,parDerv1_lambda1,parDerv1_k2,parDerv1_lambda2,parDerv1_k3,parDerv1_lambda3),6,1)
      
      parDerv2_11beta <- - t(X) %*% diag(as.vector(E2*(1+E3)/(Denom^2)),nrow=n,ncol=n) %*% X
      parDerv2_12beta <- t(X) %*% diag(as.vector((E2*E3)/(Denom^2)),nrow=n,ncol=n) %*% X
      parDerv2_21beta <- parDerv2_12beta
      parDerv2_22beta <- - t(X) %*% diag(as.vector(E3*(1+E2)/(Denom^2)),nrow=n,ncol=n) %*% X
      Derv2beta <- rbind(cbind(parDerv2_11beta,parDerv2_12beta), cbind(parDerv2_21beta,parDerv2_22beta))
      
      parDerv2_k3k3 <- sum(Q3*(delta*(-1/k3.curr^2) - (Y/lambda3.curr)^k3.curr*(log(Y/lambda3.curr))^2))
      parDerv2_k3lambda3 <- sum(Q3*(delta*(-1/lambda3.curr) + Y^k3.curr*lambda3.curr^(-k3.curr-1)*(k3.curr*log(Y/lambda3.curr)+1)))
      parDerv2_lambda3lambda3 <- sum(Q3*(delta*(k3.curr/lambda3.curr^2) - k3.curr*(k3.curr+1)*Y^k3.curr*lambda3.curr^(-k3.curr-2)))
      parDerv2_surv <- matrix(c(parDerv2_k1k1,parDerv2_k1lambda1,0,0,0,0,
                                parDerv2_k1lambda1,parDerv2_lambda1lambda1,0,0,0,0,
                                0,0,parDerv2_k2k2,parDerv2_k2lambda2,0,0,
                                0,0,parDerv2_k2lambda2,parDerv2_lambda2lambda2,0,0,
                                0,0,0,0,parDerv2_k3k3,parDerv2_k3lambda3,
                                0,0,0,0,parDerv2_k3lambda3,parDerv2_lambda3lambda3),6,6)
    }
     
    if(numGrp==4){
      parDerv1_2beta <- t(X) %*% (Q3-E3/Denom)
      parDerv1_3beta <- t(X) %*% (Q4-E4/Denom)
      Derv1beta <- rbind(parDerv1_1beta,parDerv1_2beta,parDerv1_3beta)
      
      parDerv1_k3 <- sum(Q3*(delta*(1/k3.curr + log(Y/lambda3.curr)) - (Y/lambda3.curr)^k3.curr * log(Y/lambda3.curr)))
      parDerv1_lambda3 <- sum(Q3*(delta*(-k3.curr/lambda3.curr) + k3.curr*Y^k3.curr*lambda3.curr^(-k3.curr-1)))
      parDerv1_k4 <- sum(Q4*(delta*(1/k4.curr + log(Y/lambda4.curr)) - (Y/lambda4.curr)^k4.curr * log(Y/lambda4.curr)))
      parDerv1_lambda4 <- sum(Q4*(delta*(-k4.curr/lambda4.curr) + k4.curr*Y^k4.curr*lambda4.curr^(-k4.curr-1)))
      
      parDerv1_surv <- matrix(c(parDerv1_k1,parDerv1_lambda1,parDerv1_k2,parDerv1_lambda2,parDerv1_k3,parDerv1_lambda3,parDerv1_k4,parDerv1_lambda4),8,1)
      
      parDerv2_11beta <- - t(X) %*% diag(as.vector(E2*(1+E3+E4)/(Denom^2)),nrow=n,ncol=n) %*% X
      parDerv2_12beta <- t(X) %*% diag(as.vector((E2*E3)/(Denom^2)),nrow=n,ncol=n) %*% X
      parDerv2_21beta <- parDerv2_12beta
      parDerv2_13beta <- t(X) %*% diag(as.vector((E2*E4)/(Denom^2)),nrow=n,ncol=n) %*% X
      parDerv2_31beta <- parDerv2_13beta
      parDerv2_22beta <- - t(X) %*% diag(as.vector(E3*(1+E2+E4)/(Denom^2)),nrow=n,ncol=n) %*% X
      parDerv2_23beta <- t(X) %*% diag(as.vector((E3*E4)/(Denom^2)),nrow=n,ncol=n) %*% X
      parDerv2_32beta <- parDerv2_23beta
      parDerv2_33beta <- - t(X) %*% diag(as.vector(E4*(1+E2+E3)/(Denom^2)),nrow=n,ncol=n) %*% X
      
      Derv2beta <- rbind(cbind(parDerv2_11beta,parDerv2_12beta,parDerv2_13beta), 
                         cbind(parDerv2_21beta,parDerv2_22beta,parDerv2_23beta),
                         cbind(parDerv2_31beta,parDerv2_32beta,parDerv2_33beta))
      
      parDerv2_k3k3 <- sum(Q3*(delta*(-1/k3.curr^2) - (Y/lambda3.curr)^k3.curr*(log(Y/lambda3.curr))^2))
      parDerv2_k3lambda3 <- sum(Q3*(delta*(-1/lambda3.curr) + Y^k3.curr*lambda3.curr^(-k3.curr-1)*(k3.curr*log(Y/lambda3.curr)+1)))
      parDerv2_lambda3lambda3 <- sum(Q3*(delta*(k3.curr/lambda3.curr^2) - k3.curr*(k3.curr+1)*Y^k3.curr*lambda3.curr^(-k3.curr-2)))
      parDerv2_k4k4 <- sum(Q4*(delta*(-1/k4.curr^2) - (Y/lambda4.curr)^k4.curr*(log(Y/lambda4.curr))^2))
      parDerv2_k4lambda4 <- sum(Q4*(delta*(-1/lambda4.curr) + Y^k4.curr*lambda4.curr^(-k4.curr-1)*(k4.curr*log(Y/lambda4.curr)+1)))
      parDerv2_lambda4lambda4 <- sum(Q4*(delta*(k4.curr/lambda4.curr^2) - k4.curr*(k4.curr+1)*Y^k4.curr*lambda4.curr^(-k4.curr-2)))
      
      parDerv2_surv <- matrix(c(parDerv2_k1k1,parDerv2_k1lambda1,0,0,0,0,0,0,
                                parDerv2_k1lambda1,parDerv2_lambda1lambda1,0,0,0,0,0,0,
                                0,0,parDerv2_k2k2,parDerv2_k2lambda2,0,0,0,0,
                                0,0,parDerv2_k2lambda2,parDerv2_lambda2lambda2,0,0,0,0,
                                0,0,0,0,parDerv2_k3k3,parDerv2_k3lambda3,0,0,
                                0,0,0,0,parDerv2_k3lambda3,parDerv2_lambda3lambda3,0,0,
                                0,0,0,0,0,0,parDerv2_k4k4,parDerv2_k4lambda4,
                                0,0,0,0,0,0,parDerv2_k4lambda4,parDerv2_lambda4lambda4),8,8)
      
    }
    sizeb <- 0
    beta.new <- beta.curr - solve(Derv2beta)%*%Derv1beta
    while(sum((beta.new>2)|(beta.new<(-2)))>0){
      sizeb <- sizeb-1
      beta.new <- beta.curr - (2^sizeb*solve(Derv2beta))%*%Derv1beta
    }
    sizet <- 0
    if(min(abs(diag(parDerv2_surv)))<10^(-10)){break}
    theta.new <- theta.curr - (2^sizet*solve(parDerv2_surv))%*%parDerv1_surv
    while(sum(theta.new<0)>0){
      sizet <- sizet-1
      theta.new <- theta.curr - (2^sizet*solve(parDerv2_surv))%*%parDerv1_surv
    }
    param.new <- matrix(c(theta.new,beta.new),(numGrp-1)*(dx+1)+2*numGrp,1)
    logLik.new <- logLik_iter(dataX=dataX,Y=Y,delta=delta,surv=matrix(param.new[1:(2*numGrp),],ncol=1),beta=matrix(param.new[-(1:(2*numGrp)),],ncol=1),numGrp=numGrp)
    diff <- abs(logLik.new - logLik.curr)
    diffv <- c(diffv,diff)
    theta.curr <- theta.new
    beta.curr <- beta.new
    param.curr <- param.new
    logLik.curr <- logLik.new
    logLikv <- c(logLikv,logLik.curr)
    k1.curr <- theta.new[1,];lambda1.curr <- theta.new[2,]
    k2.curr <- theta.new[3,];lambda2.curr <- theta.new[4,]
    if(numGrp==3){
      k3.curr <- theta.new[5,]
      lambda3.curr <- theta.new[6,]
    }
    if(numGrp==4){
      k3.curr <- theta.new[5,]
      lambda3.curr <- theta.new[6,]
      k4.curr <- theta.new[7,]
      lambda4.curr <- theta.new[8,]
    }
    iter <- iter+1
  }
  if (diff<1){converge=1}else{converge=0}
  return(list(logLikv=logLikv,param=param.curr,diffv=diffv,
              logLik=logLik.curr,converge=converge))
}

###########################################
#### function to calculate log-Likelihood at each iteration of EM algorithm, for different number of subgroups assumed in the data
###########################################
# @param: surv: estimated parameters for survival distributions
# beta: estimated beta coefficients
# numGrp: number of latent subgroups assumed in the data, taking values of 2,3,4
logLik_iter <- function(dataX,Y,delta,surv,beta,numGrp){
  n <- nrow(dataX)
  X <- as.matrix(cbind(rep(1,n), dataX))
  
  k1 <- surv[1,];lambda1 <- surv[2,]
  k2 <- surv[3,];lambda2 <- surv[4,]
  if(numGrp==3){k3 <- surv[5,]
    lambda3 <- surv[6,]}
  if(numGrp==4){
    k3 <- surv[5,]
    lambda3 <- surv[6,]
    k4 <- surv[7,]
    lambda4 <- surv[8,]}
  E1 <- matrix(1,nrow=n,ncol=1)
  E2 <- exp(X%*%beta[1:(dx+1),])
  if(numGrp==3){E3 <- exp(X%*%beta[(dx+2):(2*dx+2),])}
  if(numGrp==4){
    E3 <- exp(X%*%beta[(dx+2):(2*dx+2),])
    E4 <- exp(X%*%beta[(2*dx+3):(3*(dx+1)),])}
  
  if(numGrp==3){
    Denom <- E1 + E2 + E3
  }else if(numGrp==4){
    Denom <- E1 + E2 + E3 + E4
  }else{Denom <- E1 + E2}
  
  S1 <- exp(-(Y/lambda1)^k1)
  S2 <- exp(-(Y/lambda2)^k2)
  f1 <- k1*lambda1^(-k1)*Y^(k1-1)*exp(-(Y/lambda1)^k1)
  f2 <- k2*lambda2^(-k2)*Y^(k2-1)*exp(-(Y/lambda2)^k2)
  if(numGrp==3){S3 <- exp(-(Y/lambda3)^k3)
  f3 <- k3*lambda3^(-k3)*Y^(k3-1)*exp(-(Y/lambda3)^k3)}
  if(numGrp==4){
    S3 <- exp(-(Y/lambda3)^k3)
    f3 <- k3*lambda3^(-k3)*Y^(k3-1)*exp(-(Y/lambda3)^k3)
    S4 <- exp(-(Y/lambda4)^k4)
    f4 <- k4*lambda3^(-k4)*Y^(k4-1)*exp(-(Y/lambda4)^k4)}
  
  if(numGrp==3){
    res <- sum(delta*log(f1*E1/Denom + f2*E2/Denom + f3*E3/Denom) + (1-delta)*log(S1*E1/Denom + S2*E2/Denom + S3*E3/Denom))
  }else if(numGrp==4){
    res <- sum(delta*log(f1*E1/Denom + f2*E2/Denom + f3*E3/Denom + f4*E4/Denom) + (1-delta)*log(S1*E1/Denom + S2*E2/Denom + S3*E3/Denom + S4*E4/Denom))
  }else{res <- sum(delta*log(f1*E1/Denom + f2*E2/Denom) + (1-delta)*log(S1*E1/Denom + S2*E2/Denom))}
  return(res)
}

###########################################
#### function to calculate BIC
###########################################
# @param: logL: log-likelihood
# p: number of parameters
# n: number of observations in the data
BIC_cal <- function(logL,p,n){
  res = -2*logL + p*log(n)
  return(res)
}


###########################################
#### function to compute \hat{\beta} based on EM algorithm, for two latent subgroups
###########################################
# @param: dataX: covariate, input of a n*dx matrix. (sample size = n, number of covariates = dx)
# Y: observation time, input of a n*1 matrix. (sample size = n)
# delta: event status, input of a n*1 matrix, 1 = event, 0 = censored.
# initial_beta: initial estimate of beta, input of a (dx+1)*1 matrix. (number of covariates = dx)
# tau: tolerance paramter for EM algorithm
# lambda1, k1: parameters of weibull distribution for group 1. (scale param: lambda1, shape param: k1)
# lambda2, k2: parameters of weibull distribution for group 2. (scale param: lambda2, shape param: k2)
# dx: number of covariates 
## returns: beta: a dx*1 matrix of estimated beta coefficients.
# iter: number of iterations
# diffv: a vector of length iter, consisting of changes in beta coefficients in each iteration
# converge: 1=converged, 0=not converged
betaOptim2g <- function(dataX, Y,delta,initial_beta,tau,ini_lambda1,ini_lambda2,ini_k1,ini_k2,dx){
    n <- nrow(dataX)
    X <- cbind(rep(1,n), dataX)
    X <- as.matrix(X)
    
    diffv <- c()
    diff <- 1
    iter <- 0
    beta.curr <- initial_beta
    theta.curr <- matrix(c(ini_k1,ini_lambda1,ini_k2,ini_lambda2),4,1)
    param.curr <- matrix(c(theta.curr,beta.curr),dx+5,1)
    k1.curr <- ini_k1;k2.curr <- ini_k2;lambda1.curr <- ini_lambda1;lambda2.curr <- ini_lambda2
    
    while(diff>tau & diff < 10){
      E1 <- matrix(1,nrow=n,ncol=1)
      E2 <- exp(X%*%beta.curr[1:(dx+1),])
      Denom <- E1 + E2
      
      S1 <- exp(-(Y/lambda1.curr)^k1.curr)
      S2 <- exp(-(Y/lambda2.curr)^k2.curr)
      f1 <- k1.curr*lambda1.curr^(-k1.curr)*Y^(k1.curr-1)*exp(-(Y/lambda1.curr)^k1.curr)
      f2 <- k2.curr*lambda2.curr^(-k2.curr)*Y^(k2.curr-1)*exp(-(Y/lambda2.curr)^k2.curr)
      
      Q1 <- (delta*f1+(1-delta)*S1)*E1
      Q2 <- (delta*f2+(1-delta)*S2)*E2
      SQ <- Q1+Q2
      Q1 <- Q1/SQ
      Q2 <- Q2/SQ
      
      parDerv1_beta <- t(X) %*% (Q2-E2/Denom)
      
      parDerv1_k1 <- sum(Q1*(delta*(1/k1.curr + log(Y/lambda1.curr)) - (Y/lambda1.curr)^k1.curr * log(Y/lambda1.curr)))
      parDerv1_lambda1 <- sum(Q1*(delta*(-k1.curr/lambda1.curr) + k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)))
      parDerv1_k2 <- sum(Q2*(delta*(1/k2.curr + log(Y/lambda2.curr)) - (Y/lambda2.curr)^k2.curr * log(Y/lambda2.curr)))
      parDerv1_lambda2 <- sum(Q2*(delta*(-k2.curr/lambda2.curr) + k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
      parDerv1_surv <- matrix(c(parDerv1_k1,parDerv1_lambda1,parDerv1_k2,parDerv1_lambda2),4,1)
      
      parDerv2_beta <- - t(X) %*% diag(as.vector(E2*(Denom-E2)/Denom^2),nrow=n,ncol=n) %*% X
      
      parDerv2_k1k1 <- sum(Q1*(delta*(-1/k1.curr^2) - (Y/lambda1.curr)^k1.curr*(log(Y/lambda1.curr))^2))
      parDerv2_k1lambda1 <- sum(Q1*(delta*(-1/lambda1.curr) + Y^k1.curr*lambda1.curr^(-k1.curr-1)*(k1.curr*log(Y/lambda1.curr)+1)))
      parDerv2_lambda1lambda1 <- sum(Q1*(delta*(k1.curr/lambda1.curr^2) - k1.curr*(k1.curr+1)*Y^k1.curr*lambda1.curr^(-k1.curr-2)))
      parDerv2_k2k2 <- sum(Q2*(delta*(-1/k2.curr^2) - (Y/lambda2.curr)^k2.curr*(log(Y/lambda2.curr))^2))
      parDerv2_k2lambda2 <- sum(Q2*(delta*(-1/lambda2.curr)+Y^k2.curr*lambda2.curr^(-k2.curr-1)*(k2.curr*log(Y/lambda2.curr)+1)))
      parDerv2_lambda2lambda2 <- sum(Q2*(delta*(k2.curr/lambda2.curr^2) - k2.curr*(k2.curr+1)*Y^k2.curr*lambda2.curr^(-k2.curr-2)))
      parDerv2_surv <- matrix(c(parDerv2_k1k1,parDerv2_k1lambda1,0,0,
                                parDerv2_k1lambda1,parDerv2_lambda1lambda1,0,0,
                                0,0,parDerv2_k2k2,parDerv2_k2lambda2,
                                0,0,parDerv2_k2lambda2,parDerv2_lambda2lambda2),4,4)
      
      beta.new <- beta.curr - solve(parDerv2_beta)%*%parDerv1_beta
      size <- 0
      theta.new <- theta.curr - (2^size*solve(parDerv2_surv))%*%parDerv1_surv
      while(sum(theta.new<0)>0){
        size <- size-1
        theta.new <- theta.curr - (2^size*solve(parDerv2_surv))%*%parDerv1_surv
      }
      param.new <- matrix(c(theta.new,beta.new),dx+5,1)
      diff <- sqrt(t(param.new-param.curr)%*%(param.new-param.curr))
      diffv <- c(diffv,diff)
      theta.curr <- theta.new
      beta.curr <- beta.new
      param.curr <- param.new
      k1.curr <- theta.new[1,];lambda1.curr <- theta.new[2,]
      k2.curr <- theta.new[3,];lambda2.curr <- theta.new[4,]
      iter <- iter+1
    }
    if (diff>10) {converge=0} 
    if (diff<tau) {converge=1}
    return(list(beta=beta.curr,theta=theta.curr,k1=k1.curr,k2=k2.curr,lambda1=lambda1.curr,
                lambda2=lambda2.curr, iter=iter,diffv=diffv, converge=converge))
}


###########################################
###### function for variable selection, for two latent subgroups
###########################################
## @param: dataX: covariate, input of a n*dx matrix. (sample size = n, number of covariates = dx)
# Y: observation time, input of a n*1 matrix. (sample size = n)
# delta: event status, input of a n*1 matrix, 1 = event, 0 = censored.
# beta_EM: beta estimate from EM algorithm, output of function betaOptim2g. 
# nfolds: number of folds in cross validation. 
# lambda1, k1: parameters of weibull distribution for group 1. (scale param: lambda1, shape param: k1)
# lambda2, k2: parameters of weibull distribution for group 2. (scale param: lambda2, shape param: k2)
# dx: number of covariates 
## returns: beta: beta coefficients based on optimal tunning parameter lambda from cross validation
# cv.alasso: a glmnet object of cross validation fit
library(tidyr)
library(glmnet)
cv_glmalasso2g <- function(dataX,Y,delta,beta_EM,nfolds=5,lambda1_EM,lambda2_EM,k1_EM,k2_EM,dx){
  n <- nrow(dataX)
  X <- as.matrix(cbind(rep(1,n),dataX))
  S1 <- exp(-(Y/lambda1_EM)^k1_EM)
  S2 <- exp(-(Y/lambda2_EM)^k2_EM)
  f1 <- k1_EM*lambda1_EM^(-k1_EM)*Y^(k1_EM-1)*exp(-(Y/lambda1_EM)^k1_EM)
  f2 <- k2_EM*lambda2_EM^(-k2_EM)*Y^(k2_EM-1)*exp(-(Y/lambda2_EM)^k2_EM)
  
  E1 <- matrix(1,nrow=n,ncol=1)
  E2 <- exp(X%*%beta_EM[1:(dx+1),])
  Denom <- E1 + E2
  
  Q1 <- (delta*f1+(1-delta)*S1)*E1
  Q2 <- (delta*f2+(1-delta)*S2)*E2
  SQ <- Q1+Q2
  Q1 <- Q1/SQ
  Q2 <- Q2/SQ
  
  dataX$q1 <- Q1
  dataX$q2 <- Q2
  data.l <- gather(dataX,key=posterior, value=q,q1:q2)
  data.l$b[data.l$posterior=="q1"] <- 1
  data.l$b[data.l$posterior=="q2"] <- 2
  data.l$b <- as.factor(data.l$b)
  glmfit <- glmnet(x=as.matrix(data.l[,1:dx]), y=data.l$b, family="binomial", weights=data.l$q, alpha=0, lambda=0)
  beta_EM <- as.numeric(coef(glmfit))
  cv.alasso <- cv.glmnet(x=as.matrix(data.l[,1:dx]), y=data.l$b, family="binomial", weights=data.l$q, alpha=1,
                         nfolds=nfolds, penalty.factor=1/abs(as.numeric(beta_EM)),lambda=2^(seq(-16,16,0.5)))
  coef_alasso <- coef(cv.alasso,lambda=cv.alasso$lambda.1se)
  coef <- matrix(as.numeric(coef_alasso), ncol=1)
  
  return(list(beta=coef, cv.alasso=cv.alasso))
}

####################################################
#### function for post selection prediction, for two latent subgroups
####################################################
## @param: dataX: covariate, input of a n*dx matrix. (sample size = n, number of covariates = dx)
# Y: observation time, input of a n*1 matrix. (sample size = n)
# delta: event status, input of a n*1 matrix, 1 = event, 0 = censored.
# beta_al: beta estimate from adaptive lasso, output of function cv_glmalasso2g.
# tau: tolerance paramter for EM algorithm
# lambda1, k1: parameters of weibull distribution for group 1. (scale param: lambda1, shape param: k1)
# lambda2, k2: parameters of weibull distribution for group 2. (scale param: lambda2, shape param: k2)
## returns beta: a dx*1 matrix of estimated beta coefficients.
# iter: number of iterations
# diffv: a vector of length iter, consisting of changes in beta coefficients in each iteration
# x_nonzero: a vector of indice of nonzero covariates
# obsInfoM: observed information matrix based on Louis formula
# converge: 1=converged, 0=not converged
betaOptim2g_postsele <- function(dataX,Y,delta,beta_al,tau,ini_lambda1,ini_lambda2,ini_k1,ini_k2){
  n <- nrow(dataX)
  x_nonzero <- which(beta_al[,1]!=0)-1
  data.2 <- cbind(dataX[,x_nonzero[-1]],Y)
  if(is.null(ncol(data.2))){
    numbeta.2 <- 0
    X.2 <- as.matrix(cbind(rep(1,n)))
  } else {
    numbeta.2 <- ncol(data.2)-1# data.2 contains y and selected x for comparing group 2 vs. 1
    X.2 <- as.matrix(cbind(rep(1,n),data.2[,-ncol(data.2)]))
  }
  
  initial_beta <- matrix(rep(0,1+numbeta.2),ncol=1)
  diffv <- c()
  diffbv <- c()
  diffsv <- c()
  diff <- 1
  iter <- 0
  beta.curr <- initial_beta
  theta.curr <- matrix(c(ini_k1,ini_lambda1,ini_k2,ini_lambda2),4,1)
  param.curr <- matrix(c(theta.curr,beta.curr),5+numbeta.2,1)
  k1.curr <- ini_k1;k2.curr <- ini_k2;lambda1.curr <- ini_lambda1;lambda2.curr <- ini_lambda2
  logLik <- c()
  
  while(diff>tau & diff < 10){
    E1 <- matrix(1,nrow=n,ncol=1)
    E2 <- exp(X.2%*%beta.curr)
    Denom <- E1 + E2
    
    S1 <- exp(-(Y/lambda1.curr)^k1.curr)
    S2 <- exp(-(Y/lambda2.curr)^k2.curr)
    f1 <- k1.curr*lambda1.curr^(-k1.curr)*Y^(k1.curr-1)*exp(-(Y/lambda1.curr)^k1.curr)
    f2 <- k2.curr*lambda2.curr^(-k2.curr)*Y^(k2.curr-1)*exp(-(Y/lambda2.curr)^k2.curr)
    
    Q1 <- (delta*f1+(1-delta)*S1)*E1
    Q2 <- (delta*f2+(1-delta)*S2)*E2
    SQ <- Q1+Q2
    Q1 <- Q1/SQ
    Q2 <- Q2/SQ
    
    parDerv1_beta <- t(X.2) %*% (Q2-E2/Denom)
    
    parDerv1_k1 <- sum(Q1*(delta*(1/k1.curr + log(Y/lambda1.curr)) - (Y/lambda1.curr)^k1.curr * log(Y/lambda1.curr)))
    parDerv1_lambda1 <- sum(Q1*(delta*(-k1.curr/lambda1.curr) + k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)))
    parDerv1_k2 <- sum(Q2*(delta*(1/k2.curr + log(Y/lambda2.curr)) - (Y/lambda2.curr)^k2.curr * log(Y/lambda2.curr)))
    parDerv1_lambda2 <- sum(Q2*(delta*(-k2.curr/lambda2.curr) + k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    parDerv1_surv <- matrix(c(parDerv1_k1,parDerv1_lambda1,parDerv1_k2,parDerv1_lambda2),4,1)
    
    parDerv2_beta <- - t(X.2) %*% diag(as.vector(E2*(Denom-E2)/Denom^2),nrow=n,ncol=n) %*% X.2
    
    parDerv2_k1k1 <- sum(Q1*(delta*(-1/k1.curr^2) - (Y/lambda1.curr)^k1.curr*(log(Y/lambda1.curr))^2))
    parDerv2_k1lambda1 <- sum(Q1*(delta*(-1/lambda1.curr) + Y^k1.curr*lambda1.curr^(-k1.curr-1)*(k1.curr*log(Y/lambda1.curr)+1)))
    parDerv2_lambda1lambda1 <- sum(Q1*(delta*(k1.curr/lambda1.curr^2) - k1.curr*(k1.curr+1)*Y^k1.curr*lambda1.curr^(-k1.curr-2)))
    parDerv2_k2k2 <- sum(Q2*(delta*(-1/k2.curr^2) - (Y/lambda2.curr)^k2.curr*(log(Y/lambda2.curr))^2))
    parDerv2_k2lambda2 <- sum(Q2*(delta*(-1/lambda2.curr)+Y^k2.curr*lambda2.curr^(-k2.curr-1)*(k2.curr*log(Y/lambda2.curr)+1)))
    parDerv2_lambda2lambda2 <- sum(Q2*(delta*(k2.curr/lambda2.curr^2) - k2.curr*(k2.curr+1)*Y^k2.curr*lambda2.curr^(-k2.curr-2)))
    parDerv2_surv <- matrix(c(parDerv2_k1k1,parDerv2_k1lambda1,0,0,
                              parDerv2_k1lambda1,parDerv2_lambda1lambda1,0,0,
                              0,0,parDerv2_k2k2,parDerv2_k2lambda2,
                              0,0,parDerv2_k2lambda2,parDerv2_lambda2lambda2),4,4)
    
    beta.new <- beta.curr - solve(parDerv2_beta)%*%parDerv1_beta
    size <- 0
    theta.new <- theta.curr - (2^size*solve(parDerv2_surv))%*%parDerv1_surv
    while(sum(theta.new<0)>0){
      size <- size-1
      theta.new <- theta.curr - (2^size*solve(parDerv2_surv))%*%parDerv1_surv
    }
    param.new <- matrix(c(theta.new,beta.new),5+numbeta.2,1)
    diffb <- sqrt(t(beta.new-beta.curr)%*%(beta.new-beta.curr))
    diffs <- sqrt(t(theta.new-theta.curr)%*%(theta.new-theta.curr)) 
    diff <- sqrt(t(param.new-param.curr)%*%(param.new-param.curr))
    diffv <- c(diffv,diff)
    diffbv <- c(diffbv,diffb)
    diffsv <- c(diffsv,diffs)
    theta.curr <- theta.new
    beta.curr <- beta.new
    param.curr <- param.new
    k1.curr <- theta.new[1,];lambda1.curr <- theta.new[2,]
    k2.curr <- theta.new[3,];lambda2.curr <- theta.new[4,]
    iter <- iter+1
  }
  if(diff>10){converge = 0}
  if(diff<tau){converge = 1}
  if(diff<tau){
    E1 <- matrix(1,nrow=n,ncol=1)
    E2 <- exp(X.2%*%beta.new)
    Denom <- E1 + E2
    
    S1 <- exp(-(Y/lambda1.curr)^k1.curr)
    S2 <- exp(-(Y/lambda2.curr)^k2.curr)
    f1 <- k1.curr*lambda1.curr^(-k1.curr)*Y^(k1.curr-1)*exp(-(Y/lambda1.curr)^k1.curr)
    f2 <- k2.curr*lambda2.curr^(-k2.curr)*Y^(k2.curr-1)*exp(-(Y/lambda2.curr)^k2.curr)
    
    Q1 <- (delta*f1+(1-delta)*S1)*E1
    Q2 <- (delta*f2+(1-delta)*S2)*E2
    SQ <- Q1+Q2
    Q1 <- Q1/SQ
    Q2 <- Q2/SQ
    
    # expectation of - second derivative of complete likelihood
    expBk1k1 <- sum(Q1*(delta/k1.curr^2 + (Y/lambda1.curr)^k1.curr*(log(Y/lambda1.curr))^2))
    expBk1lambda1 <- sum(Q1*(delta/lambda1.curr - Y^k1.curr*lambda1.curr^(-k1.curr-1)*(k1.curr*log(Y/lambda1.curr)+1)))
    expBlambda1lambda1 <- sum(Q1*(delta*(-k1.curr/lambda1.curr^2) + (k1.curr+1)*k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-2)))
    expBk2k2 <- sum(Q2*(delta/k2.curr^2 + (Y/lambda2.curr)^k2.curr*(log(Y/lambda2.curr))^2))
    expBk2lambda2 <- sum(Q2*(delta/lambda2.curr - Y^k2.curr*lambda2.curr^(-k2.curr-1)*(k2.curr*log(Y/lambda2.curr)+1)))
    expBlambda2lambda2 <- sum(Q2*(delta*(-k2.curr/lambda2.curr^2) + (k2.curr+1)*k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-2)))
    
    expBbeta <- t(X.2) %*% diag(as.vector(E2/Denom^2),nrow=n,ncol=n) %*% X.2
    expB1 <- matrix(c(expBk1k1,expBk1lambda1,0,0,rep(0,1+numbeta.2),
                      expBk1lambda1,expBlambda1lambda1,0,0,rep(0,1+numbeta.2),
                      0,0,expBk2k2,expBk2lambda2,rep(0,1+numbeta.2),
                      0,0,expBk2lambda2,expBlambda2lambda2,rep(0,1+numbeta.2)),nrow=4,ncol=5+numbeta.2,byrow=TRUE)
    expB2 <- cbind(matrix(0,nrow=1+numbeta.2,ncol=4),expBbeta)
    expB <- rbind(expB1,expB2)
    # expectation of SS^T of complete likelihood, S: gradient vector
    expSStk1k1 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - 
                            (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr))^2)
    expSStk1lambda1 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * 
                             (delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)))
    expSStk1k2 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * 
                        Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)))
    expSStk1lambda2 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - 
                                 (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * 
                             Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    expSStk1beta <- matrix(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - 
                                 (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * (E1/Denom - 1),nrow=1) %*% X.2
    expSStlambda1lambda1 <- sum(Q1*(delta*(-k1.curr/lambda1.curr)+
                                      k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1))^2)
    expSStlambda1k2 <- sum(Q1*(delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)) * 
                             Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)))
    expSStlambda1lambda2 <- sum(Q1*(delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)) * 
                                  Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    expSStlambda1beta <- matrix(Q1*(delta*(-k1.curr/lambda1.curr)+
                                      k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)) * (E1/Denom - 1),nrow=1) %*% X.2 
    expSStk2k2 <- sum(Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr))^2)
    expSStk2lambda2 <- sum(Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - 
                                 (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)) * 
                             (delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    expSStk2beta <- matrix(Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - 
                                 (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)) * (E1/Denom),nrow=1) %*% X.2
    expSStlambda2lambda2 <- sum(Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1))^2)
    expSStlambda2beta <- matrix(Q2*(delta*(-k2.curr/lambda2.curr)+
                                      k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)) * (E1/Denom),nrow=1) %*% X.2 
    expSStbeta <- t(X.2) %*% diag(as.vector((Q1*E2^2 + Q2)/Denom^2),nrow=n,ncol=n) %*% X.2
    
    expSSt1 <- matrix(c(expSStk1k1,expSStk1lambda1,expSStk1k2,expSStk1lambda2,expSStk1beta,
                        expSStk1lambda1,expSStlambda1lambda1,expSStlambda1k2,expSStlambda1lambda2,expSStlambda2beta,
                        expSStk1k2,expSStlambda1k2,expSStk2k2,expSStk2lambda2,expSStk2beta,
                        expSStk1lambda2,expSStlambda1lambda2,expSStk2lambda2,expSStlambda2lambda2,expSStlambda2beta),nrow=4,ncol=5+numbeta.2,byrow=TRUE)
    expSSt21 <- matrix(c(expSStk1beta,expSStlambda1beta,expSStk2beta,expSStlambda2beta),nrow=1+numbeta.2,ncol=4,byrow=FALSE)
    expSSt2 <- cbind(expSSt21,expSStbeta)
    expSSt <- rbind(expSSt1,expSSt2)
    # (expectation of S)^2
    expSk1 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)))
    expSlambda1 <- sum(Q1*(delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)))
    expSk2 <- sum(Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)))
    expSlambda2 <- sum(Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    expSbeta <- t(t(X.2)%*%(E1/Denom - Q1))
    
    expS2k1k1 <- sum(Q1^2*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - 
                             (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr))^2)
    expS2k1lambda1 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * 
                            Q1*(delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)))
    expS2k1k2 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * 
                       Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)))
    expS2k1lambda2 <- sum(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * 
                            Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    expS2k1beta <- matrix(Q1*(delta*(1/k1.curr+log(Y)-log(lambda1.curr)) - 
                                (Y/lambda1.curr)^k1.curr*log(Y/lambda1.curr)) * (E1/Denom - Q1),nrow=1) %*% X.2
    expS2lambda1lambda1 <- sum(Q1^2*(delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1))^2)
    expS2lambda1k2 <- sum(Q1*(delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)) * 
                            Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)))
    expS2lambda1lambda2 <- sum(Q1*(delta*(-k1.curr/lambda1.curr)+k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)) * 
                                 Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    expS2lambda1beta <- matrix(Q1*(delta*(-k1.curr/lambda1.curr)+
                                     k1.curr*Y^k1.curr*lambda1.curr^(-k1.curr-1)) * (E1/Denom - Q1),nrow=1) %*% X.2
    expS2k2k2 <- sum(Q2^2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr))^2)
    expS2k2lambda2 <- sum(Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)) * 
                            Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)))
    expS2k2beta <- matrix(Q2*(delta*(1/k2.curr+log(Y)-log(lambda2.curr)) - 
                                (Y/lambda2.curr)^k2.curr*log(Y/lambda2.curr)) * (E1/Denom - Q1),nrow=1) %*% X.2
    expS2lambda2lambda2 <- sum(Q2^2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1))^2)
    expS2lambda2beta <- matrix(Q2*(delta*(-k2.curr/lambda2.curr)+k2.curr*Y^k2.curr*lambda2.curr^(-k2.curr-1)) * (E1/Denom - Q1), nrow=1) %*% X.2
    expS2beta <- t(X.2) %*% diag(as.vector(((-Q1*E2 + Q2)/Denom)^2),nrow=n,ncol=n) %*% X.2
    
    expS21 <- matrix(c(expS2k1k1,expS2k1lambda1,expS2k1k2,expS2k1lambda2,expS2k1beta,
                       expS2k1lambda1,expS2lambda1lambda1,expS2lambda1k2,expS2lambda1lambda2,expS2lambda1beta,
                       expS2k1k2,expS2lambda1k2,expS2k2k2,expS2k2lambda2,expS2k2beta,
                       expS2k1lambda2,expS2lambda1lambda2,expS2k2lambda2,expS2lambda2lambda2,expS2lambda2beta),nrow=4,ncol=5+numbeta.2,byrow=TRUE)
    expS221 <- matrix(c(expS2k1beta,expS2lambda1beta,expS2k2beta,expS2lambda2beta),nrow=1+numbeta.2,ncol=4,byrow=FALSE)
    expS22 <- cbind(expS221,expS2beta)
    expS2 <- rbind(expS21,expS22)
    
    obsInfoM <- expB - expSSt + expS2
  }else{obsInfoM <- NA}
  
  return(list(beta=beta.curr,k1=k1.curr,k2=k2.curr,lambda1=lambda1.curr,
              lambda2=lambda2.curr, iter=iter,diffv=diffv, converge=converge, 
              x_nonzero=x_nonzero, obsInfoM=obsInfoM))
}

########################################
##### function to generate test results, for two latent subgroups
########################################
## @param: beta: beta estimate from training set, output of function betaOptim2g.
# testX: covariates of testing dataset, input of a n*dx matrix. (sample size = n, number of covariates = dx)
# dx: number of variables (or number of important variables)
## returns: B.new: a vector of length n with estimated group membership
test_res2g <- function(beta, testX, dx){
  n <- nrow(testX)
  X.new <- as.matrix(cbind(rep(1,n),testX))
  prob2 <- exp(X.new%*%beta)/(1+exp(X.new%*%beta))
  prob1 <- 1/(1+exp(X.new%*%beta))
  maxv <- pmax(prob1,prob2)
  B.new <- rep(0,n)
  B.new[(maxv==prob1)] <- 1
  B.new[(maxv==prob2)] <- 2
  return(B.new)
}

##################################################
##### function to generate test results based on post selection prediction, for two latent subgroups
##################################################
# @param: beta_ps: estimated beta from post selection prediction
# testX: covariates of testing dataset, an input of n*q, (sample size = n, number of nonzero covariates = q)
# x_nonzero: a vector of indices of nonzero covariates, output of function betaOptim2g_postsele
# returns: B.new: a vector of length n with estimated group membership
test2g_res_postsele <- function(beta_ps, testX, x_nonzero){
  n <- nrow(testX)
  X.new.2 <- as.matrix(cbind(rep(1,n),testX[,x_nonzero[-1]]))
  beta <- beta_ps
  number.2 <- length(x_nonzero)
  prob2 <- exp(X.new.2%*%matrix(beta[1:number.2,]))/(1+exp(X.new.2%*%matrix(beta[1:number.2,])))
  prob1 <- 1/(1+exp(X.new.2%*%matrix(beta[1:number.2,])))
  prob <- cbind(prob1,prob2)
  max <- apply(prob,1,max)
  B.new <- rep(0,n)
  B.new[(max==prob[,1])] <- 1
  B.new[(max==prob[,2])] <- 2
  return(B.new)
}

#################################
#### function to evaluate variable selection results. 
#### Compute correct/incorrect zeros and proportion of nonzero/selected variables
#################################
## @param: beta_al: a list of length nrep, each component is a matrix of estimated beta coefficients from adaptive lasso (an output of function cv_glmalasso2g)
# dx: number of variables
# dxx: number of important variables
## returns: propselec12: a nrep*dx matrix, each row represents a vector of indicators whether corresponding variable is zero (1=nonzero, 0=zero). 
# corr0_12: a vector of length nrep, each elements represents number of correct zeros of one replication. 
# incorr0_12: a vector of length nrep, each elements represents number of incorrect zeros of one replication. 
prop_selec2g <- function(beta_al, dx, dxx){
  coef <- list()
  nrep <- length(beta_al)
  corr0_12 = rep(0,nrep)
  incorr0_12 = rep(0,nrep)
  propselec12 = matrix(0, ncol=dx, nrow=nrep)
  for(i in 1:nrep){
    coefi <- beta_al[[i]]
    incorr0_12[i] <- sum(coefi[2:(dxx+1)]==0)
    corr0_12[i] <- sum(coefi[(dxx+2):(dx+1)]==0)
    for (j in 1:dx){
      propselec12[i,j] <- (coefi[(j+1)]!=0)
    }
  }

  return(list(propselec12=propselec12,
              corr0_12=corr0_12,
              incorr0_12=incorr0_12))
}

