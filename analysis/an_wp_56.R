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
library(rlist)
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
library(latex2exp)
library(cowplot)
library(rootSolve)
library(FactoMineR)
library(ggsci)
library(viridis)
library(rlist)
library(latex2exp)
library(GreedyEPL)
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions.R")
source("R/gjam.R")
source("R/BNP_functions.R")
source('analysis/analysis_functions.R')
load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

##### PCA data 

set.seed(123)
PA_pdata<- load_object("Bauges_dataset/PA_data_clean_PCA.RData")
#Bauges_plant<- cbind(PA_pdata$PC1,PA_pdata$PC2, PA_pdata[,7:(ncol(PA_pdata)-2)])
#names(Bauges_plant) = c("PC1", "PC2", names(PA_pdata[,7:(ncol(PA_pdata)-2)]))
#save(Bauges_plant, file = "Bauges_plant.Rdata")

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
Colnames_Y<- tibble(CN = 1:112, CODE_CBNA=colnames(Ydata))
formula <- as.formula( ~   PC1  + PC2 + I(PC1^2) + I(PC2^2))

folderpath="PCA_analysis/r_wp/wp_56/"
folderpath_mod="PCA_analysis/r_wp/wp_56/models/chain_1/"
##################################Species names and Functional groups #####################################
#Preparing species functional groups 
Species_names_groups<- read.csv("Bauges/PFG_Bauges_Description_2017.csv", sep="\t")
true_names<- as.data.frame(names(table(Species_names_groups$PFG)))
names(true_names)<- c("PFG")
true_names$K_n<- 1:16
Species_names_groups_num<- merge(Species_names_groups,true_names, by="PFG" )
#Species_names_groups_num
True_clustering<- tibble( CODE_CBNA =  colnames(Ydata)[1: (ncol(Ydata))])
True_clust<- merge(Species_names_groups_num,True_clustering, by="CODE_CBNA")
####################################################################################
Colnames_Y<- merge(Colnames_Y,Species_names_groups_num [,c(2,3,5)], by ="CODE_CBNA" )
Colnames_Y$species<- as.character(Colnames_Y$species)
Colnames_Y$species<- strtrim(Colnames_Y$species, 20)

## Load models 1st run
fit_gjamDP1<- load_object(paste0(folderpath_mod,"fit_gjamDP1.Rdata"))
## Load models 2nd run
fit_gjamPY1<- load_object(paste0(folderpath_mod,"fit_gjamPY1.Rdata"))
#Parameters
Pars_wp_56<- list(DP1 = fit_gjamDP1$modelList$reductList, PY=  fit_gjamPY1$modelList$reductList)
### Check the prior
#### Add prior/posterior plot
prior_nu1 <- fit_gjamDP1$modelList$reductList$otherpar$shape
prior_nu2 <- fit_gjamDP1$modelList$reductList$otherpar$rate
alpha_vec<- rgamma(50000, prior_nu1,prior_nu2)
x<- sapply(alpha_vec, functionDPM,n=112,N=112)
mean(x)

# pdf(file = "Plots/DP1_alpha_posterior_K56.pdf", width=2., height =3.2)
# aDP1 =tibble(Prior =alpha_vec,
#              Posterior =fit_gjamDP1$chains$alpha.DP_g[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng])%>%
#   gather(Distribution, Probability, Prior:Posterior)%>%
#   ggplot(aes(x=Probability, fill=Distribution)) +scale_color_viridis()+
#   geom_density(adjust = 2, alpha=0.5)+
#   #ggtitle(c) +
#   xlab(TeX(sprintf('$\\alpha$')))+
#   ylab("")+
#   theme_bw() +
#   #theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
#   theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "none", plot.title = element_text(hjust = 0.5))
# #aDP1
# plot(aDP1)
# dev.off()


aDP1_56 =tibble(Prior =alpha_vec,Posterior =fit_gjamDP1$chains$alpha.DP_g[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng])%>% gather(Distribution, Probability, Prior:Posterior)
aDP1_56$K_pr <- "K = 56"
#save(aDP1_56, file =  paste0(folderpath,"aDP1_56.Rdata"))
###################################################################################################


PY_prior<- function(k,H, n, alpha, sigma, Cnk_mat){
  n_vec<- 0:(n-2)
  fal_fact <- prod(alpha +1 +n_vec)
  coef = exp(log(factorial(H)) -  log(factorial(H - k)) - log(fal_fact))
  sum<- 0
  for (l in (k:n)){
    if (k==l){
      val0 = exp( lgamma(alpha/sigma +l ) - lgamma(alpha/sigma + 1) - l*log(H)  + Cnk_mat[n,l])
    }
    else{
      val0 = exp( lgamma(alpha/sigma +l ) - lgamma(alpha/sigma + 1) - l*log(H) + Strlng2(l, k, log = TRUE) + Cnk_mat[n,l])
    }
    sum<- sum + val0
  }
  return( (coef*sum)/sigma)
}

#PY_pars[[1]]$Precomp_mat
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H05.Rdata")
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H025.Rdata")
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H025.Rdata")
load("IJulia_part/C_nk_matrix/Cnk_mat_112_H08.Rdata")
library(CryptRndTest)

x_vec<- 1:112
pks<- sapply(x_vec,PY_prior,Pars_wp_56$PY$N,Pars_wp_56$PY$N,  Pars_wp_56$PY$alpha_py, Pars_wp_56$PY$sigma_py, Cnk_112_112_H08)
plot(x_vec, pks)
exp<-sum(x_vec*pks)
exp

gjam<- load_object("PCA_analysis/r5_models/chain_1/fit_gjam.Rdata")

trace_chain_DP1 = apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x)))
trace_chain_PY = apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))
trace_chain_DP = apply(gjam$chains$kgibbs,1,function(x) length(unique(x)))

df_trace_wp_56<- tibble(it= 1: length(apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))),
                      DP1=trace_chain_DP1,
                      PY1 = trace_chain_PY,
                      DP= trace_chain_DP)
#save(df_trace_wp_56, file =  paste0(folderpath,"df_trace_wp_56.Rdata"))

#pdf(file = "Plots/Kn_posterior_trace_56.pdf", width= 8.27, height = 9.69)
df_trace_wp_56 %>%
  gather(Chains, trace, DP1:DP)%>%
  ggplot(aes(x=it,y=trace,col=Chains))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))
#dev.off()

## Dsitribution

prob_clust_DP1<-vector("numeric", length = 112)
p_ks <-table(df_trace_wp_56$DP1)
prob_clust_DP1[as.numeric(c(names(p_ks)))]<- p_ks/sum(p_ks)
prob_clust_DP1[-as.numeric(c(names(p_ks)))]<- 0
plot(1:112,prob_clust_DP1, col="red", type="l")


prob_clust_DP<-vector("numeric", length = 112)
p_ks <-table(df_trace_wp_56$DP)
prob_clust_DP[as.numeric(c(names(p_ks)))]<- p_ks/sum(p_ks)
prob_clust_DP[-as.numeric(c(names(p_ks)))]<- 0
plot(1:112,prob_clust_DP, col="black", type="l")


prob_clust_PY<-vector("numeric", length = 112)
p_ks <-table(df_trace_wp_56$PY1)
prob_clust_PY[as.numeric(c(names(p_ks)))]<- p_ks/sum(p_ks)
prob_clust_PY[-as.numeric(c(names(p_ks)))]<- 0
plot(1:112,prob_clust_PY, col="red", type="l")
lines(1:112,prob_clust_DP1, col="blue")
lines(1:112,prob_clust_DP, col="black")

DP1_post_clust_56<- prob_clust_DP1
#save(DP1_post_clust_56, file= paste0(folderpath,"DP1_post_clust_56.Rdata"))

prob_clust_PY_56<- prob_clust_PY
#save(prob_clust_PY_56, file=paste0(folderpath,"prob_clust_PY_56.Rdata"))

#prob_clustDP<- prob_clust_DP
#save(prob_clustDP, file=paste0(folderpath, "prob_clustDP.Rdata"))

#### Add prior/posterior plot
prior_nu1 <- fit_gjamDP1$modelList$reductList$otherpar$shape
prior_nu2 <- fit_gjamDP1$modelList$reductList$otherpar$rate
alpha_vec<- rgamma(80000, prior_nu1,prior_nu2)
x<- sapply(alpha_vec, functionDPM,n=112,N=112)
mean(x)



##### Posterior clustering 

## CLuster estimates

DP1_clust_full_56<- gjamClust2_full(model= fit_gjamDP1,K= 16,true_clust =True_clust$K_n )
#save(DP1_clust_full_56, file =  paste0(folderpath,"Clustering_gre_diff_sp_DP1_56.Rdata"))
length(unique(DP1_clust_full_56$VI_est[[1]]))


PY_clust_full_56<- gjamClust2_full(model= fit_gjamPY1,K= 16,true_clust =True_clust$K_n )
#save(PY_clust_full_56, file =  paste0(folderpath,"Clustering_gre_diff_sp_PY_56.Rdata"))
length(unique(PY_clust_full_56$VI_est[[3]]))


Clust2<- load_object("~/Documents/GitHub/GJAM_clust/PCA_analysis/r5/Clusters/Clusters_all_2.Rdata")

Clust2$CBN <- as.numeric(Clust2$CODE_CBNA)
Clust2_sorted= Clust2[order(Clust2$CBN),]

arandi(DP1_clust_full_56$VI_est[[3]],Clust2_sorted$ClustDP1)
arandi(PY_clust_full_56$VI_est[[3]], Clust2_sorted$ClustPY1)

WP_56_DP1<-tibble(Model="DP1",K=length(unique(DP1_clust_full_56$VI_est[[3]])),Dist_PFG = arandi(DP1_clust_full_56$VI_est[[3]], True_clust$K_n), Dist_CM = arandi(DP1_clust_full_56$VI_est[[3]], Clust2_sorted$ClustDP1))
#save(WP_56_DP1, file =  paste0(folderpath,"WP_56_DP1.Rdata"))
WP_56_PY<- tibble(Model="PY",K=length(unique(PY_clust_full_56$VI_est[[3]])),Dist_PFG = arandi(PY_clust_full_56$VI_est[[3]], True_clust$K_n), Dist_CM = arandi(PY_clust_full_56$VI_est[[3]], Clust2_sorted$ClustPY1))
#save(WP_56_PY, file =  paste0(folderpath,"WP_56_PY.Rdata"))
##################################################################