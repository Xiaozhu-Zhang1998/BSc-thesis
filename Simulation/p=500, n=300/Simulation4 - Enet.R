## Load the Packages
rm(list = ls())
library(MASS)
library(Matrix)
library(APML0)
library(Hmisc)


## ==============================
## PART I: Foundation
## ==============================
q <- 500
n <- 300
mu <- rep(0,q)
Sigma <- matrix(0, nrow = q, ncol = q)

## ## Specify rho ##
rho <- 0.9
## ## Specify epsilon ##
l_sigma <- 1

## The Covariance Matrix
for(i in 1:q)
{
  for(j in 1:q)
  {
    Sigma[i,j] <- rho ^ abs(i-j)
  }
}


## The Coefficients beta
b1 <- c(1.5, 1.4, 0.7, 0.8, 0.9)
b2 <- c(-1.2, -1.0, -0.6, -0.5, -0.7)
beta <- c()
for(i in 1:5)
{
  beta <- c(beta, b1[i])
  beta <- c(beta, rep(0,49))
  beta <- c(beta, b2[i])
  beta <- c(beta, rep(0,49))
}


## The new example x0
x0 <- read.csv("x0_rho9.csv", header = F, sep=',')
x0 <- x0$V1

## set threashold c for sqrt(n-3)*z
c <- qnorm(1-0.001, 0, 1, lower.tail = T, log.p = F)
r <- (exp(2*c/sqrt(n-3))-1)/(exp(2*c/sqrt(n-3))+1)





## ==============================
## PART VIII: Enet Encapsulation
## ==============================

## Enet: Enet(0,0,0,0,0)
Enet <- function(TP, FP, CI, MSE, RISK)
{
  ## --------------------
  ## (1.1) Generate Data
  ## --------------------
  ## The Design Matrix X
  X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  ## epsilon
  e <- mvrnorm(n = n, mu = 0, Sigma = l_sigma)
  ## Survival Time
  Time <- exp(X %*% beta + e)
  ## Status
  status <- rbinom(n, 1, 0.6)
  
  ## --------------------
  ## (1.2) Conduct the Enet
  ## --------------------
  y <- cbind(Time, status)
  colnames(y) <- c("time", "status")
  alpha <- seq(0.1, 0.9, 0.1)
  cv <- c()
  for(i in 1:9)
  {
    fit.Enet <- APML0(X, y, family = "cox",
                      penalty = "Enet", alpha = alpha[i],
                      lambda = NULL, nfolds = 5, keep.beta = FALSE)
    cv <- as.numeric(c(cv, fit.Enet$fit0[2]))
  }
  cv.min <- min(cv)
  i <- c(1:9)[cv==cv.min]
  fit.Enet <- APML0(X, y, family = "cox",
                    penalty = "Enet", alpha = alpha[i],
                    lambda = NULL, nfolds = 5, keep.beta = FALSE)
  
  ## index 1: TP - Number of true positive covariates
  tp <- sum(-fit.Enet$Beta0 > 0 & beta > 0)
  TP <- c(TP, tp)
  ## index 2: FP - Number of false positive covariates
  fp <- sum(-fit.Enet$Beta0 > 0 & beta <= 0)
  FP <- c(FP, fp)
  ## index 3: CI - Correspondence Index
  CI <- c(CI, as.numeric(1-rcorr.cens(X %*% fit.Enet$Beta0,Surv(Time, status))[1]))
  ## index 4: MSE - Mean squared error
  V <-(t(X) %*% X)/n
  mse <- t(-fit.Enet$Beta0-beta) %*% V %*% (-fit.Enet$Beta0-beta)
  MSE <- c(MSE, mse)
  ## index 5: RISK - predicted lambda
  RISK <- c(RISK, exp(t(x0) %*% (-fit.Enet$Beta0 - beta)) )
  return(list(TP = TP, FP = FP, CI = CI, MSE = MSE, RISK = RISK))
}




## ==============================
## PART IX: Replication of Enet
## ==============================
N <- 50
TP <- c()
FP <- c()
CI <- c()
MSE <- c()
RISK<- c()
for(K in 1:N)
{
  fun.Enet <- Enet(TP, FP, CI, MSE, RISK)
  TP <- fun.Enet$TP
  FP <- fun.Enet$FP
  CI <- fun.Enet$CI
  MSE <- fun.Enet$MSE
  RISK <- fun.Enet$RISK
}
## Record the CVM, TP, TN, CI, NN NOW!
mean(TP)
mean(FP)
mean(CI)
mean(MSE)
mean(RISK)


