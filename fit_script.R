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
source("R/BNP_functions.R")
source("R/rlaptrans.r")

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

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
K_prior=56
r_reduct = 5

#folderpath="PCA_analysis/test/"
folderpath="PCA_analysis/r_wp/wp_56/models/chain_2/"


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

load("IJulia_part/C_nk_matrix/Cnk_mat_112_05.Rdata")
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H05.Rdata")
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H08.Rdata")
load("IJulia_part/C_nk_matrix/Cnk_mat_112_08.Rdata")

par = compute_alpha_PYM(H=112,n=112,sigma=0.8,Mat_prior= Cnk_112_112_H08, K=K_prior)
rl2   <- list(r = r_reduct, DRtype="2" ,N=112, alpha_py=par,sigma_py=0.8,K=K_prior, Precomp_mat=Cnk_112_112_08)
ml2   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl2,PREDICTX = F)
fit_gjamPY1<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml2)

save(fit_gjamPY1, file =paste0(folderpath,"fit_gjamPY1.Rdata"))

##################################################################################################