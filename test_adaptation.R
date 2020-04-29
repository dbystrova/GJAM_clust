#test the extensions with adaptative MH

rm(list=ls())
setwd("~/Phd/Master/Code/gjam 4/")
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(Rcpp)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions.R")
source("R/gjam.R")

d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
source_data(d)
xdata <- forestTraits$xdata[,c(1,2,8)]


formula <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

y  <- gjamReZero(forestTraits$treesDeZero)  # extract y
treeYdata  <- gjamTrimY(y,10)$y             # at least 10 plots

ng=1000
burnin=500
#old version
#rl <- list(r = 8, N = 20, rate=10,shape=10)
#new version
rl <- list(r = 8, N = 20) #this is totally default version which should call basic

#rl1 <- list(r = 8, N = 20, DRtype="1",rate=10,shape=10)

rl1 <- list(r = 8, N = 20, DRtype="1") #No Kn identified, in this case default is ceiling(0.2*S)
rl1 <- list(r = 8, N = 20, DRtype="1", K=10) #Kn identified


rl2  <- list(r = 8, N = 20,DRtype="2") #here to modify N


### Need to modify to have N inside
N_eps<-floor(.compute_tau_mean(0.3,2,0.1) + 2*.compute_tau_var(0.3,2,0.1))
rl3   <- list(r = 8, N = N_eps, DRtype="3", sigma_py=0.25, alpha=2)
N_eps<-floor(.compute_tau_mean(0.5,10,0.1) + 2*.compute_tau_var(0.5,10,0.1))
rl4   <- list(r = 8, N = N_eps, DRtype="4",ro.disc=0.5) #here to modify N

ml4   <- list(ng = ng, burnin = burnin, typeNames = 'DA', reductList = rl4) #change ml
ml3   <- list(ng = ng, burnin = burnin, typeNames = 'DA', reductList = rl3) #change ml
ml2   <- list(ng = ng, burnin = burnin, typeNames = 'DA', reductList = rl2) #change ml
ml1   <- list(ng = ng, burnin = burnin, typeNames = 'DA', reductList = rl1) #change ml
ml   <- list(ng = ng, burnin = burnin, typeNames = 'DA', reductList = rl) #change ml

form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

fit <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml)
fit1 <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml1)
fit2 <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml2)
fit3 <- gjam(form,xdata =xdata, ydata = treeYdata, modelList = ml3)
fit4 <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml4)



#formula=form, xdata=xdata, ydata= treeYdata, modelList = ml4



alpha<-fit1$chains$alpha.DP_g
alpha<-mcmc(alpha)
x11()
plot(alpha)
acfplot(alpha)
cumuplot(alpha)
dev.off()

alpha<-fit2$chains$alpha.DP_g
alpha<-mcmc(alpha)
x11()
plot(alpha)
acfplot(alpha)
cumuplot(alpha)
dev.off()

#compute acceptance ratio to validate adaptation
(length(rle(fit2$chains$alpha.DP_g)[[1]])-1)/ng


alpha<-fit4$chains$alpha.PY_g
alpha<-mcmc(alpha)
x11()
plot(alpha)
acfplot(alpha)
cumuplot(alpha)
dev.off()

#compute acceptance ratio to validate adaptation
(length(rle(fit4$chains$alpha.PY_g)[[1]])-1)/ng


discount<-mcmc(fit4$chains$discount.PY_g)
x11()
plot(discount)
acfplot(discount)
cumuplot(discount)
dev.off()

last_weight <- fit1$chains$pk_g[ml1$ng,]
last_weight <- fit2$chains$pk_g[ml2$ng,]
last_weight <- fit3$chains$pk_g[ml3$ng,]
last_weight <- fit4$chains$pk_g[ml4$ng,]




#############




.setupReduct <- function(modelList, S, Q, n){
  
  REDUCT <- F
  N <- r <- NULL
  
  rl <- NULL
  
  if( 'REDUCT' %in% names(modelList) ){
    rl <- list(N = NULL, r = NULL )
    if(!modelList$REDUCT)return( rl )
  }
  
  npar <- (S+1)/2 + Q
  
  ratio <- 1/5
  N <- min(c(5, S))
  r <- N - 1
  if(npar/n > ratio){
    N <- ceiling( ( ratio*n - Q )/5 )
    if(N > 25)N <- 25
    if(N < 4)N  <- 4
    r <- ceiling( N/2 )
  }
  
  if( 'reductList' %in% names(modelList) ){
    REDUCT <- T
    rl <- modelList$reductList
    N <- rl$N
    r <- rl$r
    if(N >= S){
      N <- S - 1
      warning(' dimension reduction requires reductList$N < no. responses ')
    }
  }
  
  if( !'reductList' %in% names(modelList) ){
    rl <- list(r = r, N = N, alpha.DP = S)
  }
  rl
}




Q <- ncol(x)
n <- nrow(x)

