## Load the Packages
rm(list = ls())
library(MASS)
library(Matrix)
library(APML0)
library(Hmisc)

## ====================
## PART I: Load Data
## ====================
## --------------------
## (1.1) The Foundation
## --------------------
setwd("C:\\Users\\小竹子\\Desktop\\Thesis\\Codes\\Codes - Empirical Study")
raw.dat <- read.csv("Company.csv")
X <- raw.dat[,-c(1,2,3,4,5,6)]
q <- ncol(X)
n <- nrow(X)
y <- cbind(raw.dat$time, raw.dat$status)
colnames(y) <- c("time", "status")




## --------------------
## (1.2) The Graph Laplacian 'Omega'
## --------------------
Corr <- cor(X)
## n=243, sqrt(n-3)*z ~ N(0,1)
## set threashold c for sqrt(n-3)*z
c <- qnorm(1-0.001, 0, 1, lower.tail = T, log.p = F)
r <- (exp(2*c/sqrt(n-3))-1)/(exp(2*c/sqrt(n-3))+1)
Omega <- matrix(0,q,q)
Omega[abs(Corr) > r] <- 1
diag(Omega) <- rep(0, q)




## ====================
## PART II: Regression
## ====================
## --------------------
## (2.1) Conduct the Preliminary Ridge Regression
## --------------------
set.seed(1213)
fit.pre <- APML0(as.matrix(X), y, family = "cox", 
                 penalty = "Enet", alpha = 0)
beta.pre <- fit.pre$Beta[,50]
wbeta <- abs(1/beta.pre)
sgn <- rep(-1, q)
sgn[beta.pre < 0] <- 1

## --------------------
## (2.2) Conduct the aGLasso
## --------------------
alpha <- seq(0.1, 0.9, 0.01)
cv <- c()
for(i in 1:81)
{
  set.seed(1213)
  fit.apLasso <- APML0(as.matrix(X), y, family = "cox",
                       penalty = "Net", Omega = Omega, alpha = alpha[i],
                       lambda = NULL, nlambda = 50, rlambda = 0.0001,
                       wbeta = wbeta, sgn = sgn,
                       nfolds = 5, keep.beta = FALSE)
  cv <- as.numeric(c(cv, fit.apLasso$fit0[2]))
}
cv.min <- min(cv)
i <- c(1:81)[cv==cv.min]
# [1] 7 
set.seed(1213)
fit.apLasso <- APML0(as.matrix(X), y, family = "cox",
                     penalty = "Net", Omega = Omega, alpha = alpha[i],
                     lambda = NULL, nlambda = 50, rlambda = 0.0001,
                     wbeta = wbeta, sgn = sgn,
                     nfolds = 5, keep.beta = FALSE)




## ====================
## PART III: Results Analysis
## ====================
## --------------------
## (3.1) The Coefficients
## --------------------
fit.apLasso$Beta0
# [1]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
# [6]  0.000000000  0.000000000  0.000000000 -0.124483045 -0.598824516
# [11]  0.000000000  0.000000000  1.786282945  0.000000000 -0.297380847
# [16]  0.000000000  0.000000000  0.000000000 -0.031651957  0.000000000
# [21]  0.000000000  0.000000000  0.004866881  0.000000000  0.000000000
# [26] -0.026771392 -0.322309567 -0.220922417  0.000000000  0.140258796

## Hazard Ratio
exp(fit.apLasso$Beta0)
# [1] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
# [8] 1.0000000 0.8829532 0.5494571 1.0000000 1.0000000 5.9672307 1.0000000
# [15] 0.7427611 1.0000000 1.0000000 1.0000000 0.9688437 1.0000000 1.0000000
# [22] 1.0000000 1.0048787 1.0000000 1.0000000 0.9735838 0.7244739 0.8017789
# [29] 1.0000000 1.1505715

## --------------------
## (3.2) Other indexes
## --------------------
## lambda
fit.apLasso$fit0[1]
#        lambda
# 1 0.003536778

## alpha
alpha[i]
# [1] 0.16

## CVM
fit.apLasso$fit0[2]
#       cvm
# 1 5.187291

## C-index
1-rcorr.cens(as.matrix(X) %*% fit.apLasso$Beta0, Surv(raw.dat$time, raw.dat$status))[1]
#  C Index 
# 0.701375 

## #(Non-Zero Variables)
sum(fit.apLasso$Beta0!=0)
# [1] 10


## --------------------
## (3.3) Significance Test
## --------------------
## Preparation
index <- (1:q)[fit.apLasso$Beta0!=0]
len <- length(index)
beta.n <- fit.apLasso$Beta0[fit.apLasso$Beta0!=0]
Xn <- as.matrix(X[, index])


## Global Wald Test
I <- matrix(rep(0, len * len), len)
for(u in 1:len)
{
  for(v in 1:len)
  {
    ## working on I[u,v]
    for(i in 1:n)
    {
      if(status[i]==1)
      {
        nu <- sum((Time >= Time[i]) * Xn[,u] * Xn[,v] * exp(Xn %*% beta.n))
        de <- sum((Time >= Time[i]) * exp(Xn %*% beta.n))
        I[u,v] <- I[u,v] + nu/de
      }
    }
  }
}
## Evaluate
qChi <- as.numeric(t(beta.n) %*% I %*% beta.n)
qChi
pchisq(qChi, len, ncp=0, lower.tail = F, log.p = F)


## Local Wald Tests
## Write a function
Wald.test <- function(k)
{
  ## Prepare
  beta.k <- beta.n
  beta.k[k] <- 0
  ## Start!
  I <- matrix(rep(0, len * len), len)
  for(u in 1:len)
  {
    for(v in 1:len)
    {
      ## working on I[u,v]
      for(i in 1:n)
      {
        if(status[i]==1)
        {
          nu <- sum((Time >= Time[i]) * Xn[,u] * Xn[,v] * exp(Xn %*% beta.k))
          de <- sum((Time >= Time[i]) * exp(Xn %*% beta.k))
          I[u,v] <- I[u,v] + nu/de
        }
      }
    }
  }
  ## Evaluate
  qChi <- t(beta.n[k]-0) %*% solve(solve(I)[k,k]) %*% (beta.n[k]-0)
  p.value <- pchisq(qChi, 1, ncp=0, lower.tail = F, log.p = F)
  ## return
  return(list(Chisq = qChi, p.value = p.value))
}

for(k in 1:10)
{
  Tst <- Wald.test(k)
  print(Tst$p.value)
}




## ====================
## PART IV: Network Analysis
## ====================
## (4.1) Correlation Network
Corr.gephi <- abs(Corr)
Corr.gephi[Corr.gephi==1] <- 0
write.table(Corr.gephi,"Corr.csv", row.names=FALSE, col.names=FALSE,sep=",")

## (4.2) Adfacency Matrix
write.table(Omega,"Omega.csv", row.names=FALSE, col.names=FALSE,sep=",")

## (4.3) Boxes
## Now start loop!
alpha <- c(0.2, 0.5, 0.7, 0.9)
lambda <- c(0.1, 0.01, 0.001, 0.0001)
i=4
j=1
fit.gephi <- APML0(as.matrix(X), y, family = "cox",
                   penalty = "Net", Omega = Omega, wbeta = wbeta, sgn = sgn,
                   alpha = alpha[i], lambda = lambda[j])
(1:q)[as.numeric(fit.gephi$Beta)!=0]




## ====================
## PART V: Dynamic Projection
## ====================
## Calculate H0(t)
Time <- raw.dat$time
status <- raw.dat$status
t <- sort(Time)
t <- unique(t)
H0 <- c()
for(i in 1:length(t))  
{
  term <- 0
  for(j in 1:n)
  {
    if(Time[j] <= t[i] & status[j]==1)  # Ti[j]=='ti', t[i]=='t'
    {
      term <- term + 1/sum(exp(as.matrix(X) %*% fit.apLasso$Beta0) * as.numeric(Time>=Time[j]))
    }
  }
  H0 <- c(H0, term)
}
## Calculate S0(t)=exp(-H0(t))
S0 <- exp(-H0)
## Save S0
write.table(data.frame(t,S0),"S0.csv", row.names=FALSE, col.names=FALSE,sep=",")

## Calcualte S(t|X)
S.t.X <- c()
for(i in 1:length(S0))
{
  S.t.X <- cbind(S.t.X, S0[i]^as.numeric(exp(as.matrix(X) %*% fit.apLasso$Beta0)))
} 

S.t.X.df <- data.frame(rbind(t, S.t.X))
row.names(S.t.X.df) <- c("time", as.character(raw.dat$证券代码))
write.table(S.t.X.df, "S_t_X.csv", row.names=TRUE, col.names=FALSE,sep=",")



## Reconstruct the dataset
S1 <- t(S.t.X)
S2 <- c(as.numeric(S1[,74]), as.numeric(S1[,84]), as.numeric(S1[,191]), 
        as.numeric(S1[,233]), as.numeric(S1[,241]))
S2 <- cbind(rep(t, 5), S2)
S2 <- cbind(S2, c(rep(1, length(t)), rep(2, length(t)), rep(3, length(t)),
                  rep(4, length(t)), rep(5, length(t)) ))
S2 <- data.frame(S2)
colnames(S2) <- c("Time", "SP", "Company")
S2$Company <- as.factor(S2$Company)
levels(S2$Company) <- c('中国一重','贝因美','美的集团','惠达卫浴','中宠股份')

## Plot
ggplot(S2, aes(Time, SP, color = Company)) + 
  geom_step(size = 1.5, linetype = 1) +
  labs(title="Survival Curves", y = "Survival Probability") + 
  scale_x_continuous(breaks = c(1000,2000,3000,4000,5000,6000,7000)) + 
  scale_y_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line(size = 1, colour="black"))



## Calculate S(t|bar(X))
a=t(apply(X,2,mean)) %*% fit.apLasso$Beta0
