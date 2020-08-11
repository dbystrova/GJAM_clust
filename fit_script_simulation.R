## Fitting the models
#rm(list=ls())
#setwd("~/Documents/GitHub/GJAM_clust")
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(Rcpp)
library(ggplot2)
library(AUC)
library(formattable)
library(mcclust.ext)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(factoextra)
library(Hmsc)
library(knitr)
library(tidyverse)
library(corrplot)
library(rootSolve)
library(FactoMineR)
library(ggsci)
library(viridis)
library(rust)
library(gtools)
library(CryptRndTest)
Rcpp::sourceCpp('src/cppFns.cpp')
Rcpp::sourceCpp('src/user_fns.cpp')

source("R/gjamHfunctions.R")
source("R/gjam.R")
source("BNP_functions.R")
source("rlaptrans.r")

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}



###### Simulated data

generate_data<-function(Sp=50,nsamples=500,qval=20,Ktrue=4){
  S<-Sp   #species number
  n<- nsamples  #number of samples
  q<- qval      #number of columns in true A matrix
  env<-runif(-50,50,n=n)
  X<-cbind(1,poly(env,2)) #nxK
  # idx<-sample(S)
  # B_0<-seq(0,100,length.out=S)[idx]
  # B_1<-seq(0,100,length.out=S)[idx]
  # B_2<-seq(0,100,length.out=S)[idx]
  # B<-cbind(B_0,B_1,B_2) #SxK
  x0<- c(rep(1, floor(S/3)),rep(0, floor(S/3)))
  x1<- c(rep(0, floor(S/3)),rep(1, floor(S/3)))
  x2<- c(rep(0, floor(S/3)),rep(0, floor(S/3)))
  x0f<- c(x0,rep(0, S-length(x0)))
  x1f<- c(x1,rep(0, S-length(x1)))
  x2f<- c(x2,rep(1, S-length(x2)))
  B<- cbind(x0f,x1f,x2f)
  
  L<-X%*%t(B) #nxS
  
  K_t<- Ktrue
  cat("True number of clusters : ",K_t,"\n")
  A<-matrix(NA,nrow=ceiling(K_t),ncol =q)
  for(i in 1:ceiling(K_t)){
    A[i,]<-mvrnorm(n = 1,rep(0,q), Sigma=3*diag(q)) #Nxq short and skinny
  }
  idx<-sample((1:ceiling(K_t)),S,replace=T)  # true clustering
  Lambda<-A[idx,] #Sxr tall and skinny
  Sigma<-Lambda%*%t(Lambda)+0.1*diag(S) #SxS
  Sigma_true<-Sigma
  Y<-mvrnorm(n = n, mu=rep(0,S), Sigma=Sigma)
  #change here for Y_new
  e<- rnorm(n, 0,1)
  Y_new<- L+Y+e
  xdata<-as.data.frame(X[,-1])
  colnames(xdata)<-c("env1","env2")
  return(list(xdata=xdata, Y=Y_new,idx=idx,S_true=Sigma_true))
}

data_set<- generate_data(Sp=112,nsamples=500,qval=5,Ktrue=10)

xdata<-data_set$xdata

Y<-data_set$Y
#input clustering
idx<- data_set$idx
#input Sigma 
Sigma_true<- data_set$S_true
# given formaul
formula<-as.formula(~env1+env2)


iterations= 2000
burn_period= 500

rl <- list(r =5, N = 100)
ml   <- list(ng = iterations, burnin = burn_period, typeNames = 'CON', reductList = rl,PREDICTX = F)
fit_gjam<-gjam(formula, xdata = xdata, ydata = Y, modelList = ml)


K_chain<- apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))


load("IJulia_part/Cnk_mat_112_05.Rdata")
load("IJulia_part/Cnk_mat_112_025.Rdata")
load("IJulia_part/Cnk_mat_112_H05.Rdata")
load("IJulia_part/Cnk_mat_112_H025.Rdata")
par = compute_alpha_PYM(H=112,n=112,sigma=0.25,Mat_prior= Cnk_112_112_H025, K=10)
rl2   <- list(r = 5, DRtype="2" ,N=112, alpha_py=par,sigma_py=0.25,K=10, Precomp_mat=Cnk_112_112_025)
ml2   <- list(ng = iterations, burnin = burn_period, typeNames = 'CON', reductList = rl2,PREDICTX = F)
fit_gjamPY1<-gjam(formula, xdata = xdata, ydata = Y, modelList = ml2)

arandi(fit_gjamPY1$chains$kgibbs[2000,], data_set$idx)

K_chain_2<- apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))
plot(1:2000, K_chain_2)
####################

##### PCA data 
set.seed(123)
PA_pdata<- load_object("Bauges_dataset/PA_data_clean_PCA.RData")
train_ind <- load_object( "Bauges_dataset/PCAtrain_ind.Rds")

y<- PA_pdata[,7:(ncol(PA_pdata)-2)]
Ydata<- gjamTrimY(y,20)$y 
xdata_train <- PA_pdata[train_ind, (ncol(PA_pdata)-1):(ncol(PA_pdata))]
xdata_test <- PA_pdata[-train_ind, (ncol(PA_pdata)-1):(ncol(PA_pdata))]
Ydata_train<- Ydata[train_ind,]
Ydata_test<- Ydata[-train_ind,]
S<- ncol(Ydata_train)
S_prev <- colSums(Ydata, na.rm = TRUE, dims = 1)
p_w<- S_prev[1:(length(S_prev))]/sum(S_prev[1:(length(S_prev))])

formula <- as.formula( ~   PC1  + PC2 + I(PC1^2) + I(PC2^2))

iterations=80000
burn_period=30000
K_prior=16
r_reduct = 5

folderpath="PCA_analysis/test/"


rl <- list(r =r_reduct, N = S)
ml   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl,PREDICTX = F)
fit_gjam<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml)

save(fit_gjam, file = paste0(folderpath,"fit_gjam.Rdata"))

##################################################################################################

rl1 <- list(r = r_reduct, N = S,DRtype="1", K=K_prior) #prior is the number of plant functional groups
ml1   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl1,PREDICTX = F) #change ml
fit_gjamDP1<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml1)

save(fit_gjamDP1, file =paste0(folderpath,"fit_gjamDP1.Rdata"))

##################################################################################################
# PYM
load("IJulia_part/Cnk_mat_112_025.Rdata")
load("IJulia_part/Cnk_mat_112_H025.Rdata")

load("IJulia_part/Cnk_mat_112_05.Rdata")
load("IJulia_part/Cnk_mat_112_H05.Rdata")
par = compute_alpha_PYM(H=112,n=112,sigma=0.5,Mat_prior= Cnk_112_112_H05, K=K_prior)
rl2   <- list(r = r_reduct, DRtype="2" ,N=112, alpha_py=par,sigma_py=0.5,K=K_prior, Precomp_mat=Cnk_112_112_05)
ml2   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl2,PREDICTX = F)
fit_gjamPY1<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml2)

save(fit_gjamPY1, file =paste0(folderpath,"fit_gjamPY1.Rdata"))

##################################################################################################
# PY_SB
#rl3   <- list(r = r_reduct, DRtype="3" ,sigma_py=0.25,K=K_prior)
#ml3   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl3,PREDICTX = F)
#fit_gjamPY2<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml3)


#save(fit_gjamPY2, file =paste0(folderpath,"fit_gjamPY2.Rdata"))

########






