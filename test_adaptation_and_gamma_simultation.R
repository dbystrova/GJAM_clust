#Test the extensions with adaptative MH
rm(list=ls())
#setwd("~/Phd/Master/Code/gjam 4/")
setwd("~/Documents/GitHub/GJAM_clust")
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

## Get the Forest data 
d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
source_data(d)
xdata <- forestTraits$xdata[,c(1,2,8)]

### Define the formula for covariates
formula <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

### Delete most rare species
y  <- gjamReZero(forestTraits$treesDeZero)  # extract y
treeYdata  <- gjamTrimY(y,10)$y             # at least 10 plots


### Set the parameters for fitting the MCMC [number of iterations(ng), and burnin ]
ng=500
burnin=100


## We define the list, which define the type of dimension reduction we want to use
rl <- list(r = 8, N = 20) #this is totally default version which should call basic model
##DR type 1, the dimension reduction with the truncation version of DP
rl1 <- list(r = 8, N = 20, DRtype="1", K=10) #K prior identified

##DR type 2, the dimension reduction with the truncation version of DP
##DR type 2, the dimension reduction with the truncation version of DP
rl2  <- list(r = 8, N = 20,DRtype="2", K=15) #here to modify N 

##DR type 3, the dimension reduction with the truncation version of PY
rl3   <- list(r = 8, DRtype="3", sigma_py=0.25, K=15)


##DR type 4, the dimension reduction with the truncation version of PY
rl4   <- list(r = 8, DRtype="4",ro.disc=0.5, K=15) #here to modify N
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



#formula=form; xdata=xdata; ydata= treeYdata; modelList = ml1



alpha<-fit1$chains$alpha.DP_g
alpha<-mcmc(alpha)
x11()
plot(alpha)
acfplot(alpha)
cumuplot(alpha)
#dev.off()

alpha<-fit2$chains$alpha.DP_g
alpha<-mcmc(alpha)
x11()
plot(alpha)
acfplot(alpha)
cumuplot(alpha)
dev.off()

ng
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

###########################################################################################






.setupReduct <- function(modelList, S, Q, n){
  
  N <- r <- rl <- NULL
  
  if( 'REDUCT' %in% names(modelList) | 'reductList' %in% names(modelList) ){
    rl <- list(N = NULL, r = NULL )
    if(!modelList$REDUCT)return( rl )
  }
  
  if(n < 2*S | S > 200){
    N  <- round( S/3 )
    if(N > 25)N <- 25
    if(N <= 4)N <- 4
    r  <- ceiling( N/2 )
    rl <- list(r = r, N = N, alpha.DP = S)
    warning( 'dimension reduction' )
  }
  rl
}




a <- 2
b <- NULL

!is.null(a) && !is.null(b) && a==b

is.null(b) ||  a==b

##########################################################################################







