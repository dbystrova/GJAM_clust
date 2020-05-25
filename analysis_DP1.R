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
source("BNP_functions.R")
source('analysis_functions.R')
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


K_prior=16
r_reduct = 5

folderpath="PCA_analysis/r5/"
folderpath2="PCA_analysis/r5_2/"
folderpath3="PCA_analysis/r5_3/"
folderpath4="PCA_analysis/r5_4/"
folderpath5="PCA_analysis/r5_5/"
folderpath6="PCA_analysis/r5_6/"

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

Colnames_Y<- merge(Colnames_Y,Species_names_groups_num [,c(2,3,5)], by ="CODE_CBNA" )
Colnames_Y$species<- as.character(Colnames_Y$species)
Colnames_Y$species<- strtrim(Colnames_Y$species, 20)


## Load models 1st run
fit_gjamDP1<- load_object(paste0(folderpath3,"fit_gjamDP1.Rdata"))
## Load models 2nd run
fit_gjamDP1_2<- load_object(paste0(folderpath4,"fit_gjamDP1.Rdata"))

########################################Prediction########################################################
#Prediction out-of sample on xtest
new_out <- list(xdata =xdata_test,  nsim = 5000) # effort unchanged 
#Prediction in-sample
new_in <- list(xdata =xdata_train,  nsim = 5000) # effort unchanged 
#Conditional prediction 

### Two chains combined 
fit_DP1_comb =fit_gjamDP1
fit_DP1_comb$modelList$ng=fit_DP1_comb$modelList$ng + (fit_DP1_comb$modelList$ng - fit_DP1_comb$modelList$burnin)
fit_DP1_comb$chains$bFacGibbs = rbind(fit_gjamDP1$chains$bFacGibbs, fit_gjamDP1_2$chains$bFacGibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,])
fit_DP1_comb$chains$bgibbs = rbind(fit_gjamDP1$chains$bgibbs, fit_gjamDP1_2$chains$bgibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,])
fit_DP1_comb$chains$bgibbsUn = rbind(fit_gjamDP1$chains$bgibbsUn, fit_gjamDP1_2$chains$bgibbsUn[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,])
fit_DP1_comb$chains$fSensGibbs = rbind(fit_gjamDP1$chains$bFacGibbs, fit_gjamDP1_2$chains$bFacGibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,])
fit_DP1_comb$chains$kgibbs = rbind(fit_gjamDP1$chains$kgibbs, fit_gjamDP1_2$chains$kgibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,])
fit_DP1_comb$chains$sgibbs = rbind(fit_gjamDP1$chains$sgibbs, fit_gjamDP1_2$chains$sgibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,])
fit_DP1_comb$chains$sigErrGibbs= c(fit_gjamDP1$chains$sigErrGibbs, fit_gjamDP1_2$chains$sigErrGibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng])
fit_DP1_comb$parameters$betaMu =  (fit_gjamDP1$parameters$betaMu +fit_gjamDP1_2$parameters$betaMu)/2

gjamDP1<- model_prediction_summary(model=fit_DP1_comb, list_out_s=new_out,list_in_s=new_in, Ytest=Ydata_test, Ytrain=Ydata_train)

#save(gjamDP1, file =  paste0(folderpath,"DP1_prediction.Rdata"))
#gjamDP1<- load_object( paste0(folderpath,"DP1_prediction.Rdata"))

gjamDP1_1<- model_prediction_summary(model=fit_gjamDP1, list_out_s=new_out,list_in_s=new_in, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamDP1_2<- model_prediction_summary(model=fit_gjamDP1_2, list_out_s=new_out,list_in_s=new_in, Ytest=Ydata_test, Ytrain=Ydata_train)

##Correlation list
Cor_DP1<- list()

#Convergence of parameters
Cor_DP1<- list.append(Cor_DP1, cor(gjamDP1_1$AUC_out, gjamDP1_2$AUC_out), cor(gjamDP1_1$AUC_in, gjamDP1_2$AUC_in))
#save(Cor_DP1, file =  paste0(folderpath,"CorDP1_prediction.Rdata"))

#### Save prediction object 
################################################################################################################

### Convergence

chain_1 <- fit_gjamDP1$chains$sgibbs[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng,]
chain_2 <- fit_gjamDP1_2$chains$sgibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,]
chains  = mcmc.list(mcmc(chain_1), mcmc(chain_2))
GR_value =gelman.diag(chains)
summary(GR_value$psrf[,1])
hist(GR_value$psrf[,1])
GR_value$psrf[GR_value$psrf[,1]>1.1,1]
names<- names(GR_value$psrf[GR_value$psrf[,1]>1.1,1])

##### Ploting chains with large psf
plots<- list()
for(key in names) {
  non_c_chains_1 <- fit_gjamDP1$chains$sgibbs[(fit_gjamDP1$modelList$burnin +1):fit_gjamDP1$modelList$ng,c(key)]
  non_c_chains_2 <- fit_gjamDP1_2$chains$sgibbs[(fit_gjamDP1_2$modelList$burnin +1):fit_gjamDP1_2$modelList$ng,c(key)]
  non_c_chain  = mcmc.list(mcmc(non_c_chains_1), mcmc(non_c_chains_2))
  x =gelman.diag(non_c_chain)
  x$psrf[1,]
  gelman.plot(non_c_chain)
  p = tibble(iterations = 1:(fit_gjamDP1_2$modelList$ng - fit_gjamDP1_2$modelList$burnin),
             chain_1 =non_c_chains_1,
             chain_2 = non_c_chains_2 ) %>%
    gather(Chains, trace, chain_1:chain_2)%>%
    ggplot(aes(x=iterations,y=trace,col=Chains))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
    xlab("iterations")+ylab("Sigma coeffcient") +theme_bw()
  plots<- list.append(plots,p)
}

prow <- plot_grid(
  plots[[1]] + theme(legend.position='none'),
  plots[[2]]+ theme(legend.position='none'),
  nrow = 2, ncol=1)

legend_b <- get_legend(plots[[2]]+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))

#pdf(file = "Plots/DP1_Rhat_large_sigma.pdf", width= 8.27, height = 9.69)
plot(p)
#dev.off()
###########################################################################

df_sigma<- tibble(ES=effectiveSize(mcmc(rbind(chain_1, chain_2))))
df_sigma%>% ggplot(aes(ES, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.5) 
#summary(df_sigma)
#min =100
#save(df_sigma, file =  paste0(folderpath,"Conv_sigma_DP1.Rdata"))
###### Convergence for beta coefficients
other_values<- c("other_PC1", "other_PC2","other_I(PC1^2)","other_I(PC2^2)","other_intercept")
chain_beta_1 <- fit_gjamDP1$chains$bgibbs[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng,1:(dim(fit_gjamDP1$chains$bgibbs)[2]-5)]
chain_beta_2 <- fit_gjamDP1_2$chains$bgibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,1:(dim(fit_gjamDP1_2$chains$bgibbs)[2]-5)]
beta_chains  = mcmc.list(mcmc(chain_beta_1), mcmc(chain_beta_2))
GR_value_beta =gelman.diag(beta_chains)
summary(GR_value_beta$psrf[,1])
hist(GR_value_beta$psrf[,1])
GR_value_beta$psrf[GR_value_beta$psrf[,1]>1.1,1]
names<- names(GR_value_beta$psrf[GR_value_beta$psrf[,1]>1.1,1])

plots<- list()
for(key in names) {
  non_c_chains_1 <- fit_gjamDP1$chains$bgibbs[(fit_gjamDP1$modelList$burnin +1):fit_gjamDP1$modelList$ng,c(key)]
  non_c_chains_2 <- fit_gjamDP1_2$chains$bgibbs[(fit_gjamDP1_2$modelList$burnin +1):fit_gjamDP1_2$modelList$ng,c(key)]
  non_c_chain  = mcmc.list(mcmc(non_c_chains_1), mcmc(non_c_chains_2))
  x =gelman.diag(non_c_chain)
  print(x$psrf[1,])
  gelman.plot(non_c_chain)
  p = tibble(iterations = 1:(fit_gjamDP1_2$modelList$ng - fit_gjamDP1_2$modelList$burnin),
             chain_1 =non_c_chains_1,
             chain_2 = non_c_chains_2 ) %>%
    gather(Chains, trace, chain_1:chain_2)%>%
    ggplot(aes(x=iterations,y=trace,col=Chains))+geom_line(alpha=0.5)+ scale_color_viridis(discrete=TRUE)+
    xlab("iterations")+ylab("Sigma coeffcient") +theme_bw()
  plots<- list.append(plots,p)
}

prow <- plot_grid(
  plots[[1]] + theme(legend.position='none'),
  plots[[2]]+ theme(legend.position='none'),
  plots[[3]]+ theme(legend.position='none'),
  nrow = 3, ncol=1)
legend_b <- get_legend(plots[[3]]+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))

#pdf(file = "Plots/DP1_Rhat_large_beta.pdf", width= 8.27, height = 9.69)
plot(p)
#dev.off()

df_beta<- tibble(ES= effectiveSize(mcmc(rbind(chain_beta_1, chain_beta_2))))
df_beta %>%ggplot(aes(ES))+ geom_histogram( position="identity", alpha=0.5) 
#save(df_beta, file =  paste0(folderpath,"Conv_beta_DP1.Rdata"))

###########################################################################################################
############################Trace plot##########################################################################
trace_chain = fit_DP1_comb$chains$kgibbs
trace_chain_1 = apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x)))
trace_chain_2 = apply(fit_gjamDP1_2$chains$kgibbs,1,function(x) length(unique(x)))

df_trace_DP1<- tibble(it= 1: length(apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x)))),
                      DP1=trace_chain_1,
                      DP1_2 = trace_chain_2)
#save(df_trace_DP1, file =  paste0(folderpath,"DP1_k_chains.Rdata"))

df_trace_DP1 %>%
  gather(Chains, trace, DP1:DP1_2)%>%
  ggplot(aes(x=it,y=trace,col=Chains))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))

######################################### alpha chains ########################################
### convergence

alpha_ch_1<- fit_gjamDP1$chains$alpha.DP_g[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng]
alpha_ch_2<- fit_gjamDP1_2$chains$alpha.DP_g[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng]

df_alpha_trace_DP1<- tibble(it= 1: (fit_gjamDP1$modelList$ng - fit_gjamDP1$modelList$burnin),
                      chain1=alpha_ch_1,
                      chain2 = alpha_ch_2)

df_alpha_trace_DP1 %>%
  gather(Chains, trace, chain1:chain2)%>%
  ggplot(aes(x=it,y=trace,col=Chains))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title=TeX(sprintf('Trace plot for DP1 parameter $\\alpha$')))+xlab("iterations")+ylab(TeX(sprintf('$\\alpha$'))) +theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15)) +theme(legend.text=element_text(size=15))


alpha_chains  = mcmc.list(mcmc(alpha_ch_1), mcmc(alpha_ch_2))
gelman.diag(alpha_chains)
#1.02       1.08
gelman.plot(alpha_chains)
lattice::xyplot(alpha_chains, autoburinin=FALSE)
effectiveSize(mcmc(alpha_ch_1))
#413.1239 
##### Prior/ Posterior distribution 



#### Add prior/posterior plot
prior_nu1 <- fit_gjamDP1_2$modelList$reductList$otherpar$shape
prior_nu2 <- fit_gjamDP1_2$modelList$reductList$otherpar$rate
alpha_vec<- rgamma(100000, prior_nu1,prior_nu2)
x<- sapply(alpha_vec, functionDPM,n=112,N=112)
mean(x)

###Posterior for alpha parameter DP2

aDP1 =tibble(Prior =alpha_vec,
             Posterior =c(alpha_ch_1,alpha_ch_2))%>%
  gather(Distribution, density, Prior:Posterior)%>%
  ggplot(aes(x=density, fill=Distribution)) +scale_color_viridis()+
  geom_density(adjust = 1, alpha=0.6)+ggtitle(TeX(sprintf('Prior and posterior distribution for DP1 $\\alpha$'))) + xlab("")+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))

aDP1
#prior_nu1 <- fit_gjamPY2$modelList$reductList$otherpar$shape
#prior_nu2 <- fit_gjamPY2$modelList$reductList$otherpar$rate
#alpha_vecPY<- rgamma(18000, prior_nu1,prior_nu2)
#x<- sapply(alpha_vecPY, functionPY,n=112, sigma_py=0.25)
#mean(x)


#pdf(file = "Plots/Posterior_for_DP1_parameters.pdf", width= 8.27, height = 9.69)
plot(aDP1)
#dev.off()


##################################Clustering #####################################
sp_num<- ncol(Ydata)-1
K_chain_1 <- fit_gjamDP1$chains$kgibbs[(fit_gjamDP1$modelList$burnin +1):fit_gjamDP1$modelList$ng,]
K_chain_2 <- fit_gjamDP1_2$chains$kgibbs[(fit_gjamDP1_2$modelList$burnin +1):fit_gjamDP1_2$modelList$ng,]
K_chain = rbind(K_chain_1,K_chain_2 )
# MatDP<- relabel_clust(K_chain[,1:sp_num])
# CM_DP<-  comp.psm(MatDP)
# mbind_DP.ext <- minbinder.ext(CM_DP,MatDP, method="all",include.greedy = TRUE)
# vi_DP.ext <- minVI(CM_DP,MatDP, method="all", include.greedy = TRUE)
# length(unique(vi_DP.ext$cl[1,]))
# DP_grEPL<- MinimiseEPL(MatDP, pars =list("VI"))

load("PCA_analysis/r5/Clusters_modells_1.Rdata")
load("PCA_analysis/r5/Clusters_modells_2.Rdata")
# 
# DP_grEPL$decision
# arandi(DP_grEPL$decision, Cluster_models_1$ClustPY1)
# arandi(vi_DP.ext$cl[1,], Cluster_models_1$ClustPY1)
# arandi(DP_grEPL$decision, Cluster_models_2$ClustPY1)
# arandi(vi_DP.ext$cl[1,], Cluster_models_2$ClustPY1)

DP1_clust<- gjamClust_full_test(model= fit_DP1_comb)
#save(DP1_clust, file =  paste0(folderpath,"Clustering_mcclust_diff_sp_DP1.Rdata"))
DP1_clust<- gjamClust(model= fit_DP1_comb)

DP1_clust_1<- gjamClust(model= fit_gjamDP1)
DP1_clust_2<- gjamClust(model= fit_gjamDP1_2)
arandi(DP1_clust_1$VI_est,DP1_clust_2$VI_est)
arandi(DP1_clust_2$VI_est,DP1_clust$VI_est)
arandi(DP1_clust$VI_est,Cluster_models_1$ClustPY1)
arandi(DP1_clust_2$VI_est,DP1_clust$VI_est)

#### Compare the obtained estimates with the PFG clusters

SW_fin_table_DP1<-tibble(Model="DP1",K=length(unique(DP1_clust$VI_est)),VI_dist = vi.dist(DP1_clust$VI_est, True_clust$K_n), AR_dist = arandi(DP1_clust$VI_est, True_clust$K_n))
#save(SW_fin_table_DP1, file =  paste0(folderpath,"SW_tab_DP1.Rdata"))
Cluster_DP1_1<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],ClustDP1=DP1_clust$VI_est)
#save(Cluster_DP1_1, file =  paste0(folderpath,"Cluster_DP1_1.Rdata"))

## Second clustering method
###########################

DP1_clust2_full<- gjamClust2_full(model= fit_DP1_comb)
#save(DP1_clust2_full, file =  paste0(folderpath,"Clustering_gre_diff_sp_DP1.Rdata"))

DP1_clust2_1<- gjamClust2(model= fit_gjamDP1)
DP1_clust2_2<- gjamClust2(model= fit_gjamDP1_2)

arandi(DP1_clust2_1$VI_est,DP1_clust2_2$VI_est)
arandi(DP1_clust2_2$VI_est, Cluster_models_1$ClustDP1)
arandi(DP1_clust2_1$VI_est, Cluster_models_2$ClustDP1)
arandi(DP1_clust2_full$VI_est[[3]], DP1_clust2$VI_est)

### using true clust initialization
DP1_clust2 <- gjamClust2(model= fit_DP1_comb,  pars =list(decision_init= True_clust$K_n, loss_type="VI"))

GRE_fin_table_DP1<-tibble(Model="PY1",K=length(unique(DP1_clust2$VI_est)),VI_dist = vi.dist(DP1_clust2$VI_est, True_clust$K_n), AR_dist = arandi(DP1_clust2$VI_est, True_clust$K_n))
#save(GRE_fin_table_DP1, file =  paste0(folderpath,"GRE_tab_DP1.Rdata"))
Cluster_DP1_2<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],ClustDP1=DP1_clust2$VI_est)
#save(Cluster_DP1_2, file =  paste0(folderpath,"Cluster_DP1_2.Rdata"))
##################################Covariance matrix######################################################################
#### Add cluster labels 
Colnames_Y_clust<- merge(Colnames_Y, Cluster_DP1_2, by ="CODE_CBNA")

### Covariance matrix for the mean 
#pdf(file = "Plots/Correlation_matrix_DP1.pdf", width= 8.27, height = 9.69)
MDP= (fit_gjamDP1$parameters$corMu + fit_gjamDP1_2$parameters$corMu)/2
Colnames_Y_clust$Sp_name_DP1 <- paste(Colnames_Y_clust$ClustDP1, Colnames_Y_clust$species,sep="_")
rownames(MDP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP1"],"other")
colnames(MDP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP1"],"other")
colors_vir=viridis(length(unique(Colnames_Y_clust$ClustDP1))+1, option = "magma")
LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_DP1),"ClustDP1"],length(unique(Colnames_Y_clust$ClustDP1))+1), function(x) colors_vir[x])
cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))


corrplot(MDP, diag = FALSE, order = "hclust", tl.cex = 0.45,tl.srt=45, method = "color", tl.col=LabelCol,col=cols(200),
         type = "full", title= "Correlation for the DP1 model (original)", mar=c(0,0,1,0))

corrplot(MDP, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45,  tl.col=LabelCol,
         method = "color",col=cols(200), type = "full",title= "Correlation for the DP1 model (DP1 groups)", mar=c(0,0,1,0))
#dev.off()



######## Precision matrix ################################################################################
library(qgraph)
library(matlib)
#inv_cor = solve(MDP)
inv_cor = inv(MDP)

PC_mat <- wi2net(inv_cor)
corrplot(PC_mat, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45,  tl.col=LabelCol,
         method = "color",col=cols(200), type = "full",title= "Correlation for the DP1 model (DP1 groups)", mar=c(0,0,1,0))

### Get precision matrix from the Correlation
PC_mat[abs(PC_mat)> 0.1] 
## values  are  super  small 
########################################################################################
#### Final Table
form<-c(formula)
Fin_tab_DP1<-as.data.frame(matrix(NA,nrow=10,ncol=6))
names(Fin_tab_DP1)<- c("Parameter","DP1","r", "iter", "burn","formula")
Fin_tab_DP1$iter<- fit_gjamDP1$modelList$ng
Fin_tab_DP1$burn<- fit_gjamDP1$modelList$burnin
Fin_tab_DP1$r<-fit_gjamDP1$modelList$reductList$r
Fin_tab_DP1$formula<-as.character(form)
Fin_tab_DP1[1,1]<- "DIC"
Fin_tab_DP1[1,2]<- fit_gjamDP1$fit$DIC/10000
Fin_tab_DP1[2,1]<- "mean AUC"
Fin_tab_DP1[2,2]<- mean(gjamDP1$AUC_out)
Fin_tab_DP1[3,1]<- "AUC in"
Fin_tab_DP1[3,2]<- mean(gjamDP1$AUC_in)
Fin_tab_DP1[4,1]<- "AR dist VI loss 1"
Fin_tab_DP1[4,2]<- SW_fin_table_DP1$AR_dist
Fin_tab_DP1[5,1]<- "VI dist VI loss"
Fin_tab_DP1[5,2]<- SW_fin_table_DP1$VI_dist
Fin_tab_DP1[6,1]<- "AR dist VI loss 2"
Fin_tab_DP1[6,2]<- GRE_fin_table_DP1$AR_dist
Fin_tab_DP1[7,1]<- "VI dist VI loss 2"
Fin_tab_DP1[7,2]<- GRE_fin_table_DP1$VI_dist
Fin_tab_DP1[8,1]<- "mean K"
Fin_tab_DP1[8,2]<- mean(apply(fit_DP1_comb$chains$kgibbs,1,function(x) length(unique(x)))[fit_DP1_comb$modelList$burnin:fit_PY_comb$modelList$ng])
Fin_tab_DP1[9,1]<- "K VI"
Fin_tab_DP1[9,2]<- SW_fin_table_DP1$K
Fin_tab_DP1[10,1]<- "K VI2"
Fin_tab_DP1[10,2]<- GRE_fin_table_DP1$K
Fin_tab_DP1[,2]<- round(Fin_tab_DP1[,2], 3)
#save(Fin_tab_DP1, file = paste0(folderpath,"Fin_tabDP1_r5.Rdata"))
#library("xlsx")
# Write the first data set in a new workbook
#write.xlsx(Fin_all, file = "Final_table.xlsx")
#########################################################################################
