## test file
library(repmis)
library(gjam)
library(MASS)
library(truncnorm)
library(coda)
library(RcppArmadillo)
library(arm)
library(Rcpp)
library(raster)
library(ggplot2)
library(rgdal)
library(biomod2)
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
library("FactoMineR")
library("ggsci")

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

iterations=1000
burn_period=300
K_prior=16
r_reduct = 5

folderpath="PCA_analysis/test"

##conditional prediction
columns<-1:ncol(Ydata_train)
ycs<- sample(columns, 10)
y_c_p <-columns[ !columns %in% ycs]

##### Model spceification

rl <- list(r =r_reduct, N = S)
ml   <- list(ng = iterations, burnin = burn_period, typeNames = 'PA', reductList = rl,PREDICTX = F)
fit_gjam<-gjam(formula, xdata = xdata_train, ydata = Ydata_train, modelList = ml)

save(fit_gjam, file = paste0(folderpath,"fit_gjam.Rdata"))



