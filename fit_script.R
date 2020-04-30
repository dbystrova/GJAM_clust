## Fitting the models
#rm(list=ls())
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

Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions.R")
source("R/gjam.R")

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

iterations=30000
burn_period=12000
K_prior=16
r_reduct = 5

folderpath="PCA_analysis/r5_25/"

##conditional prediction
columns<-1:ncol(Ydata_train)
ycs<- sample(columns, 10)
y_c_p <-columns[ !columns %in% ycs]
##### Model spceification

rl <- list(r =r_reduct, N = S)
ml   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl,PREDICTX = F)
fit_gjam<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml)
save(fit_gjam, file = paste0(folderpath,"fit_gjam.Rdata"))


#rl1  <- list(r = r_reduct, N = 150, DRtype="1", K=K_prior) #prior is N
#ml1   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl1,PREDICTX = F) #change ml
#fit_gjamDP1<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml1)
#save(fit_gjamDP1, file = paste0(folderpath,"fit_gjamDP1.Rdata"))

rl2  <- list(r = r_reduct, N = S,DRtype="2", K=K_prior) #prior is Number of plant functional groups
ml2   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl2,PREDICTX = F) #change ml
fit_gjamDP2<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml2)
save(fit_gjamDP2, file =paste0(folderpath,"fit_gjamDP2.Rdata"))

##### PY1
rl3   <- list(r = r_reduct, DRtype="3" ,sigma_py=0.25,K=K_prior)
ml3   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl3,PREDICTX = F)
fit_gjamPY1<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml3)
save(fit_gjamPY1, file = paste0(folderpath,"fit_gjamPY1.Rdata"))
##### PY2
rl4   <- list(r = r_reduct,DRtype="4",ro.disc=0.5, sigma_py=0.25, K=K_prior)
ml4   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl4,PREDICTX = F)
fit_gjamPY2<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml4)
save(fit_gjamPY2, file = paste0(folderpath,"fit_gjamPY2.Rdata"))

