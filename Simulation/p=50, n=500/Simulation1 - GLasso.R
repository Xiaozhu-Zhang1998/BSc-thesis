## Load the Packages
rm(list = ls())
library(MASS)
library(Matrix)
library(APML0)
library(Hmisc)


## ==============================
## PART I: Foundation
## ==============================
q <- 50
n <- 500
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
  beta <- c(beta, 0, 0, 0, 0)
  beta <- c(beta, b2[i])
  beta <- c(beta, 0, 0, 0, 0)
}




## The new example x0
x0 <- read.csv("x0_rho9.csv", header = F, sep=',')
x0 <- x0$V1

## The weight
beta.pre <- read.csv("1_wbeta_rho9_60.csv", header = F, sep = ',')
beta.pre <- beta.pre$V1
wbeta <- abs(1/beta.pre)

## set threashold c for sqrt(n-3)*z
c <- qnorm(1-0.001, 0, 1, lower.tail = T, log.p = F)
r <- (exp(2*c/sqrt(n-3))-1)/(exp(2*c/sqrt(n-3))+1)




## ==============================
## PART II: GLASSO Encapsulation
## ==============================

## GLASSO: GLASSO(0,0,0,0,0)
GLASSO <- function(wbeta, TP, FP, CI, MSE, RISK)
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
  ## (1.2) The Graph Laplacian 'Omega'
  ## --------------------
  ## (a) 0-1 Adjacency
  # Corr <- cor(X)
  # Omega <- matrix(0,q,q)
  # Omega[abs(Corr) > r] <- 1
  # diag(Omega) <- rep(0, q)
  
  ## (b) Correlation
  # Omega <- cor(X)
  # Omega[Omega < 0] <- 0
  
  ## (c) Covariance
  Omega <- cov(X)^3
  
  ## --------------------
  ## (1.3) Conduct the GLasso
  ## --------------------
  y <- cbind(Time, status)
  colnames(y) <- c("time", "status")
  alpha <- seq(0.1, 0.9, 0.1)
  cv <- c()
  for(i in 1:9)
  {
    fit.gLasso <- APML0(X, y, family = "cox",
                        penalty = "Net", Omega = Omega, alpha = alpha[i], 
                        lambda = NULL, nlambda = 50, rlambda = 0.0001,
                        wbeta = wbeta,
                        nfolds = 5, keep.beta = FALSE)
    cv <- as.numeric(c(cv, fit.gLasso$fit0[2]))
  }
  cv.min <- min(cv)
  i <- c(1:9)[cv==cv.min]
  fit.gLasso <- APML0(X, y, family = "cox",
                      penalty = "Net", Omega = Omega, alpha = alpha[i], 
                      lambda = NULL, nlambda = 50, rlambda = 0.0001,
                      wbeta = wbeta,
                      nfolds = 5, keep.beta = FALSE)
  
  ## --------------------
  ## (1.4) Store the Performance
  ## --------------------
  ## index 1: TP - Number of true positive covariates
  tp <- sum(-fit.gLasso$Beta0 > 0 & beta > 0)
  TP <- c(TP, tp)
  ## index 2: FP - Number of false positive covariates
  fp <- sum(-fit.gLasso$Beta0 > 0 & beta <= 0)
  FP <- c(FP, fp)
  ## index 3: CI - Correspondence Index
  CI <- c(CI, as.numeric(1-rcorr.cens(X %*% fit.gLasso$Beta0,Surv(Time, status))[1]))
  ## index 4: MSE - Mean squared error
  V <-(t(X) %*% X)/n
  mse <- t(-fit.gLasso$Beta0-beta) %*% V %*% (-fit.gLasso$Beta0-beta)
  MSE <- c(MSE, mse)
  ## index 5: RISK - predicted lambda
  RISK <- c(RISK, exp(t(x0) %*% (-fit.gLasso$Beta0-beta)) )
  return(list(TP = TP, FP = FP, CI = CI, MSE = MSE, RISK = RISK))
}




## ==============================
## PART III: Replication of GLASSO
## ==============================
## Start!
N <- 50
TP <- c()
FP <- c()
CI <- c()
MSE <- c()
RISK<- c()
for(K in 1:N)
{
  fun.glasso <- GLASSO(wbeta, TP, FP, CI, MSE, RISK)
  TP <- fun.glasso$TP
  FP <- fun.glasso$FP
  CI <- fun.glasso$CI
  MSE <- fun.glasso$MSE
  RISK <- fun.glasso$RISK
}
## Record the TP, FP, CI, MSE, RISK NOW!
mean(TP)
mean(FP)
mean(CI)
mean(MSE)
mean(RISK)




