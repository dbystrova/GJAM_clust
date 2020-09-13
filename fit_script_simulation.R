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



###### Simulated data

generate_data<-function(Sp=50,nsamples=500,qval=20,Ktrue=4){
  S<-Sp   #species number
  n<- nsamples  #number of samples
  q<- qval      #number of columns in true A matrix
  env<-runif(-50,50,n=n)
  X<-cbind(1,poly(env,2)) #nxK
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
# given formula
formula<-as.formula(~env1+env2)


folderpath="Simulation/"

iterations= 5000
burn_period= 2500

rl <- list(r =5, N = 112)
ml   <- list(ng = iterations, burnin = burn_period, typeNames = 'CON', reductList = rl,PREDICTX = F)
fit_gjam<-gjam(formula, xdata = xdata, ydata = Y, modelList = ml)


K_chain<- apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))
arandi(fit_gjam$chains$kgibbs[2000,], data_set$idx)
DP_clust <- gjamCluster(fit_gjam, K=10, prior_clust =data_set$idx  )

arandi(DP_clust$VI_est[[1]], data_set$idx)

## DP1


rl1 <- list(r = 5, N = 112,DRtype="1", K=10) #prior is the number of plant functional groups
ml1   <- list(ng = iterations, burnin = burn_period, typeNames = 'CON', reductList = rl1,PREDICTX = F) #change ml
fit_gjamDP1<-gjam(formula, xdata = xdata, ydata = Y, modelList = ml1)


K_chain_DP1<- apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x)))
arandi(fit_gjamDP1$chains$kgibbs[2000,], data_set$idx)
DP1_clust <- gjamCluster(fit_gjamDP1, K=10, prior_clust =data_set$idx  )

arandi(DP1_clust$VI_est[[1]], data_set$idx)

#

load("IJulia_part/C_nk_matrix/Cnk_mat_112_05.Rdata") ## Cnk matrix() for use in the model fitting
load("IJulia_part/C_nk_matrix/Cnk_mat_112_025.Rdata")
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H05.Rdata") ##Cnk matrix for computing the alpha parameter
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H025.Rdata")
par = compute_alpha_PYM(H=112,n=112,sigma=0.25,Mat_prior= Cnk_112_112_H025, K=30)
rl2   <- list(r = 5, DRtype="2" ,N=112, alpha_py=par,sigma_py=0.25,K=10, Precomp_mat=Cnk_112_112_025)
ml2   <- list(ng = iterations, burnin = burn_period, typeNames = 'CON', reductList = rl2,PREDICTX = F)
fit_gjamPY1<-gjam(formula, xdata = xdata, ydata = Y, modelList = ml2)

arandi(fit_gjamPY1$chains$kgibbs[2000,], data_set$idx)

K_chain_2<- apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))
plot(1:2000, K_chain_2)

PY_clust <- gjamCluster(fit_gjamPY1, K=10, prior_clust =data_set$idx  )

arandi(PY_clust$VI_est[[3]], data_set$idx)




###### Check convergence 
library(coda)

gjam_mc<- mcmc(fit_gjam$chains$sgibbs[2500:5000,])
hist(effectiveSize(gjam_mc), main="ess(sigma) gjam",lwd=2,col=gray(.6),breaks=10)
plot(gjam_mc[,100])


###### Check RMSE
###### Check 






simulation_fun_gjam<-function(data_set,Sp, Ntr, rval,nsamples=500, Ktrue,q=20, it=1000, burn=500){
  S<-Sp
  n<- nsamples
  r <- rval
  iterations<-it
  
  K=sum(S/(S+(1:S)-1)) #104, his prior number of clusters when alpha=S
  cat("Prior expected number of clusters : ",Ktrue,"\n")
  K_t= Ktrue
  xdata<-data_set$xdata
  Y<-data_set$Y
  idx<- data_set$idx
  Sigma_true<- data_set$S_true
  formula<-as.formula(~env1+env2)
 
  rl   <- list(r = r, N = N_eps,rate=rate,shape=shape,V1=1,ro.disc=ro.disc) #here to modify N
  ml<-list(ng=it,burnin=burn,typeNames='CON',reductList=rl)
  fit<-.gjam_4(formula,xdata,ydata=as.data.frame(Y),modelList = ml)
  alpha.chains<-fit$chains$alpha.PY_g
  sigma.chains<-fit$chains$discount.PY_g
  pk_chains<- fit$chains$pk_g
  Ntr<-N_eps+1
  alpha.DP<-alpha.PY
  trace<-apply(fit$chains$kgibbs,1,function(x) length(unique(x)))
  ind_trace<- seq(1,it,by=1)
  trace_short<- trace[ind_trace]
  df<-as.data.frame(trace)
  df$iter<-1:it
  #####Alpha plot
  df_alpha <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
  df_alpha$alpha<- alpha.chains[-c(1:burn)]
  df_alpha$type<- "posterior"
  df_alpha_prior <- data.frame(matrix(NA, nrow =it-burn, ncol =1))
  #df_alpha_prior$alpha<- rgamma(it-burn, shape, rate)
  alpha_seq= seq(min(alpha.chains[-c(1:burn)]),max(alpha.chains[-c(1:burn)]),length=it-burn)
  df_alpha_prior$alpha <- dgamma(alpha_seq,rate,shape)
  
  df_alpha_prior$type<- "prior"
  df_alpha_all<- rbind(df_alpha[-1,],df_alpha_prior[-1,])
  ###Compute mean
  mu <- ddply(df_alpha_all, "type", summarise, grp.mean=mean(alpha))
  mu$grp.mean[which(mu$type=="prior")]=alpha.DP
  p_alpha_2<- ggplot(df_alpha, aes(x=alpha)) + geom_vline(data=mu, aes(xintercept=grp.mean, color=type),linetype="dashed")+
    geom_density(color="red")+labs(title=paste0("Posterior distribution for alpha"), caption=paste0("Number of iterations: ",it," burnin: ",burn," number of samples: ",nsamples," S=",S," ,r=",r," true gr K=",K_t, " ,N=",Ntr))+
    theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
    scale_color_manual(name = c("Legend"), values = c("prior"="#9999FF", "posterior"= "#FF6666"), labels=c("posterior mean","prior mean"))
  plot(p_alpha_2)
  
 
  N_dim<-(it-burn)
  sigma<-array(dim=c(Sp,Sp,N_dim))
  for(j in 1:N_dim){
    K<-fit$chains$kgibbs[j,]
    Z  <- matrix(fit$chains$sgibbs[j,],Ntr-1,r)
    sigma[,,j] <- .expandSigma(fit$chains$sigErrGibbs[j], Sp, Z = Z, fit$chains$kgibbs[j,], REDUCT = T) #sigma
  }
  sigma_mean<-apply(sigma,c(1,2),mean)
  err<-sum((sigma_mean-Sigma_true)^2)/(Sp*Sp)
  rmspe<-fit$fit$rmspeAll
  return(list(trace=trace_short,
              idx=idx,K=fit$chains$kgibbs[it,],
              alpha=alpha.DP,alpha.chains=alpha.chains,
              coeff_t=Sigma_true,coeff_f=sigma_mean,
              err=err,fit=rmspe))
  
}
####### Just one possible test case
#sim<-simulation_fun(Sp=50, Ntr=150, rval=3,nsamples=500, Ktrue=4,it=1000,burn=200)
# plot(as.vector(sim$coeff_t),as.vector(sim$coeff_f))
# x11()
# heatmap(sim$coeff_f)
# x11()
# heatmap(sim$coeff_t)
# plot(sim$trace)
# plot(sim$idx,sim$K)
#possible parameters to add in the loop:
# - type (T) of gjam to be fitted
# - n, the number of simulated normals
# - N, the truncation level
# - S, the number of species
# - K_t, the true number of clusters

###########Simulation for Continous data case :small S K=4###################################################



#####################################Simulation 2 K=10#######################################

####Small S, N==S, n=500
list5=list()
list0=list()
data_list=list()
lk<-list()
S_vec<-c(100,200)
r_vec<-5
#n_vec<-c(10)
n_samples<- 500
k<-1
it<-5000
burn<-2500
Ktr<-10
q<-20

path<- "~/Documents/GitHub/gjamed/sigma_post/"
for(i in 1:length(S_vec)){
  data_list=list()
  k=1
  for(l in (1:1)){
    data_list<- list.append(data_list,generate_data(Sp=S_vec[i],nsamples=n_samples,qval=q,Ktrue=Ktr))
    names(data_list)[[l]]<-paste0("S_",S_vec[i],"_q_",q,"n_",n_samples,"_K_",Ktr,"_l",l)
  }
  ########gjam 4  model list########################    
  l5<-list()
  l5<- lapply(data_list,simulation_fun_gjam4,Sp=S_vec[i], Ntr=150, q=q,rval=r_vec,nsamples=n_samples, Ktrue=Ktr,it=it,burn=burn)
  list5<-list.append(list5,assign(paste0("S_",S_vec[i],"_r_5_N_n_500_K",Ktr),l5))
  names(list5)[[k]]<-paste0("S_",S_vec[i],"_r_5_N_150_n_",n_samples,"_K",Ktr)
  save(list5, file = paste0(path,"Sigma_mod",S_vec[i],"K",Ktr,"_type4.Rda"))
  k=k+1
}

