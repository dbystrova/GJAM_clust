rm(list=ls())
setwd("~/GitHub/GJAMF/gjam 4")
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

rl <- list(r = 8, N = 20, rate=10,shape=10,V=1)
rl1 <- list(r = 8, N = 20, DRtype="1",rate=10,shape=10)
rl2  <- list(r = 8, N = 20,DRtype="2",rate=10,shape=10,V=1) #here to modify N
N_eps<-floor(.compute_tau_mean(0.3,2,0.1) + 2*.compute_tau_var(0.3,2,0.1))
rl3   <- list(r = 8, N = N_eps, DRtype="3", sigma_py=0.3, alpha=2)
N_eps<-floor(.compute_tau_mean(0.5,10,0.1) + 2*.compute_tau_var(0.5,10,0.1))
rl4   <- list(r = 8, N = N_eps, DRtype="4", rate=10,shape=10,V1=1,ro.disc=0.5) #here to modify N

ml4   <- list(ng = 200, burnin = 100, typeNames = 'DA', reductList = rl4) #change ml
ml3   <- list(ng = 200, burnin = 100, typeNames = 'DA', reductList = rl3) #change ml
ml2   <- list(ng = 200, burnin = 100, typeNames = 'DA', reductList = rl2) #change ml
ml1   <- list(ng = 200, burnin = 100, typeNames = 'DA', reductList = rl1) #change ml
ml   <- list(ng = 200, burnin = 100, typeNames = 'DA', reductList = rl) #change ml

form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )

fit <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml)
fit1 <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml1)
fit2 <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml2)
fit3 <- gjam(form,xdata =xdata, ydata = treeYdata, modelList = ml3)
fit4 <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml4)

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

alpha<-fit4$chains$alpha.PY_g
alpha<-mcmc(alpha)
x11()
plot(alpha)
acfplot(alpha)
cumuplot(alpha)
dev.off()

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

