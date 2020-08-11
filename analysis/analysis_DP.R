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



#iterations=1000
#burn_period=300
K_prior=16
r_reduct = 5

folderpath="PCA_analysis/r5/DP_analysis/"
folderpath2="PCA_analysis/r5_2/"
folderpath3="PCA_analysis/r5_3/"
folderpath4="PCA_analysis/r5_4/"

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
fit_gjamDP<- load_object("PCA_analysis/r5_3/fit_gjam.Rdata")
## Load models 2nd run
fit_gjamDP_2<- load_object("PCA_analysis/r5_4/fit_gjam.Rdata")

########################################Prediction########################################################
#Prediction out-of sample on xtest
new_out <- list(xdata =xdata_test,  nsim = 1000) # effort unchanged 
#Prediction in-sample
new_in <- list(xdata =xdata_train,  nsim = 1000) # effort unchanged 
#Conditional prediction 

### Two chains combined 
fit_DP_comb =fit_gjamDP
fit_DP_comb$modelList$ng=fit_DP_comb$modelList$ng + (fit_DP_comb$modelList$ng - fit_DP_comb$modelList$burnin)
fit_DP_comb$chains$bFacGibbs = rbind(fit_DP_comb$chains$bFacGibbs, fit_gjamDP_2$chains$bFacGibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,])
fit_DP_comb$chains$bgibbs = rbind(fit_DP_comb$chains$bgibbs, fit_gjamDP_2$chains$bgibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,])
fit_DP_comb$chains$bgibbsUn = rbind(fit_DP_comb$chains$bgibbsUn, fit_gjamDP_2$chains$bgibbsUn[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,])
fit_DP_comb$chains$fSensGibbs = rbind(fit_DP_comb$chains$bFacGibbs, fit_gjamDP_2$chains$bFacGibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,])
fit_DP_comb$chains$kgibbs = rbind(fit_DP_comb$chains$kgibbs, fit_gjamDP_2$chains$kgibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,])
fit_DP_comb$chains$sgibbs = rbind(fit_DP_comb$chains$sgibbs, fit_gjamDP_2$chains$sgibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,])
fit_DP_comb$chains$sigErrGibbs= c(fit_DP_comb$chains$sigErrGibbs, fit_gjamDP_2$chains$sigErrGibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng])
fit_DP_comb$parameters$betaMu =  (fit_gjamDP$parameters$betaMu +fit_gjamDP_2$parameters$betaMu)/2

gjamDP<- model_prediction_summary(model=fit_DP_comb, list_out_s=new_out,list_in_s=new_in, Ytest=Ydata_test, Ytrain=Ydata_train)

#save(gjamDP, file =  paste0(folderpath,"DP_prediction.Rdata"))


gjamDP_1<- model_prediction_summary(model=fit_gjamDP, list_out_s=new_out,list_in_s=new_in, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamDP_2<- model_prediction_summary(model=fit_gjamDP_2, list_out_s=new_out,list_in_s=new_in, Ytest=Ydata_test, Ytrain=Ydata_train)

##Correlation list
Cor_DP<- list()

#Convergence of parameters
Cor_DP<- list.append(Cor_DP, cor(gjamDP_1$AUC_out, gjamDP_2$AUC_out), cor(gjamDP_1$AUC_in, gjamDP_2$AUC_in))
#save(Cor_DP, file =  paste0(folderpath,"CorDP_prediction.Rdata"))

#### Save prediction object 
################################################################################################################

### Convergence

chain_1 <- fit_gjamDP$chains$sgibbs[(fit_gjamDP$modelList$burnin+1):fit_gjamDP$modelList$ng,]
chain_2 <- fit_gjamDP_2$chains$sgibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,]
chains  = mcmc.list(mcmc(chain_1), mcmc(chain_2))
GR_value =gelman.diag(chains)
summary(GR_value$psrf[,1])
hist(GR_value$psrf[,1])
#save(GR_value, file =  paste0(folderpath,"GR_value_sigma_DP.Rdata"))

###########################################################################

df_sigma<- tibble(ES=effectiveSize(mcmc(rbind(chain_1, chain_2))))
df_sigma%>% ggplot(aes(ES, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 
#summary(df_sigma)
#min =100
#save(df_sigma, file =  paste0(folderpath,"Conv_sigma_DP.Rdata"))
###### Convergence for beta coefficients
other_values<- c("other_PC1", "other_PC2","other_I(PC1^2)","other_I(PC2^2)","other_intercept")
chain_beta_1 <- fit_gjamDP$chains$bgibbs[(fit_gjamDP$modelList$burnin+1):fit_gjamDP$modelList$ng,1:(dim(fit_gjamDP$chains$bgibbs)[2]-5)]
chain_beta_2 <- fit_gjamDP_2$chains$bgibbs[(fit_gjamDP_2$modelList$burnin+1):fit_gjamDP_2$modelList$ng,1:(dim(fit_gjamDP_2$chains$bgibbs)[2]-5)]
beta_chains  = mcmc.list(mcmc(chain_beta_1), mcmc(chain_beta_2))
GR_value_beta =gelman.diag(beta_chains)
summary(GR_value_beta$psrf[,1])
hist(GR_value_beta$psrf[,1])
#save(GR_value_beta, file =  paste0(folderpath,"GR_value_beta_DP.Rdata"))

df_beta<- tibble(ES= effectiveSize(mcmc(rbind(chain_beta_1, chain_beta_2))))
df_beta %>%ggplot(aes(ES))+ geom_histogram( position="identity", alpha=0.5) 
save(df_beta, file =  paste0(folderpath,"Conv_beta_DP.Rdata"))

###########################################################################################################
############################Trace plot##########################################################################
#pdf(file = "Plots/Trace_plot_partitions.pdf", width= 8.27, height = 9.69)
trace_chain = fit_DP_comb$chains$kgibbs
trace_chain_1 = apply(fit_gjamDP$chains$kgibbs,1,function(x) length(unique(x)))
trace_chain_2 = apply(fit_gjamDP_2$chains$kgibbs,1,function(x) length(unique(x)))

df_trace_DP<- tibble(it= 1: length(apply(fit_gjamDP$chains$kgibbs,1,function(x) length(unique(x)))),
                      DP_1=trace_chain_1,
                      DP_2 = trace_chain_2)
#save(df_trace_DP, file =  paste0(folderpath,"DP_k_chains.Rdata"))

df_trace_DP %>%
  gather(Chains, trace, DP_1:DP_2)%>%
  ggplot(aes(x=it,y=trace,col=Chains))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))


##################################Clustering #####################################
sp_num<- ncol(Ydata)-1
K_chain_1 <- fit_gjamDP$chains$kgibbs[(fit_gjamDP$modelList$burnin +1):fit_gjamDP$modelList$ng,]
K_chain_2 <- fit_gjamDP_2$chains$kgibbs[(fit_gjamDP_2$modelList$burnin +1):fit_gjamDP_2$modelList$ng,]
K_chain = rbind(K_chain_1,K_chain_2 )
# MatDP<- relabel_clust(K_chain[,1:sp_num])
# CM_DP<-  comp.psm(MatDP)
# mbind_DP.ext <- minbinder.ext(CM_DP,MatDP, method="all",include.greedy = TRUE)
# vi_DP.ext <- minVI(CM_DP,MatDP, method="all", include.greedy = TRUE)
# length(unique(vi_DP.ext$cl[1,]))
# DP_grEPL<- MinimiseEPL(MatDP, pars =list("VI"))

# load("PCA_analysis/r5/Clusters_modells_1.Rdata")
# load("PCA_analysis/r5/Clusters_modells_2.Rdata")
# 
# DP_grEPL$decision
# arandi(DP_grEPL$decision, Cluster_models_1$ClustPY1)
# arandi(vi_DP.ext$cl[1,], Cluster_models_1$ClustPY1)
# arandi(DP_grEPL$decision, Cluster_models_2$ClustPY1)
# arandi(vi_DP.ext$cl[1,], Cluster_models_2$ClustPY1)

DP_clust_full<- gjamClust_full_test(model= fit_DP_comb)
#save(DP_clust_full, file =  paste0(folderpath,"Clustering_mcclust_diff_sp_DP.Rdata"))
DP_clust<- gjamClust(model= fit_DP_comb)

DP_clust_1<- gjamClust(model= fit_gjamDP)
DP_clust_2<- gjamClust(model= fit_gjamDP_2)
#arandi(PY1_clust_1$VI_est,PY1_clust_2$VI_est)

#arandi(PY1_clust_2$VI_est,Cluster_models_2$ClustPY1)
#### Compare the obtained estimates with the PFG clusters

SW_fin_table_DP<-tibble(Model="DP",K=length(unique(DP_clust$VI_est)),VI_dist = vi.dist(DP_clust$VI_est, True_clust$K_n), AR_dist = arandi(DP_clust$VI_est, True_clust$K_n))
#save(SW_fin_table_DP, file =  paste0(folderpath,"SW_tab_DP.Rdata"))
Cluster_DP_1<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],ClustDP=DP_clust$VI_est)
#save(Cluster_DP_1, file =  paste0(folderpath,"Cluster_DP_1.Rdata"))

## Second clustering method
###########################

DP_clust2_full<- gjamClust2_full(model= fit_DP_comb)
#save(DP_clust2_full, file =  paste0(folderpath,"Clustering_gre_diff_sp_DP.Rdata"))

DP_clust2_1<- gjamClust2(model= fit_gjamDP)
DP_clust2_2<- gjamClust2(model= fit_gjamDP_2)


arandi(DP_clust2_1$VI_est,DP_clust2_2$VI_est)
arandi(DP_clust2_2$VI_est, Cluster_models_1$ClustDP)
arandi(PY1_clust2_2$VI_est, Cluster_models_2$ClustDP)
arandi(DP_clust2_full$VI_est[[1]], DP_clust2_1$VI_est)

### using true clust initialization
DP_clust2 <- gjamClust2(model= fit_DP_comb,  pars =list(loss_type="VI"))


GRE_fin_table_DP<-tibble(Model="DP",K=length(unique(DP_clust2$VI_est)),VI_dist = vi.dist(DP_clust2$VI_est, True_clust$K_n), AR_dist = arandi(DP_clust2$VI_est, True_clust$K_n))
#save(GRE_fin_table_DP, file =  paste0(folderpath,"GRE_tab_DP.Rdata"))
Cluster_DP_2<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],ClustDP=DP_clust2$VI_est)
#save(Cluster_DP_2, file =  paste0(folderpath,"Cluster_DP_2.Rdata"))
Cluster_DP_2<- load_object( paste0(folderpath,"Cluster_DP_2.Rdata"))


K_chain <- apply(fit_DP_comb$chains$kgibbs[(fit_DP_comb$modelList$burnin +1):fit_DP_comb$modelList$ng,],1,function(x) length(unique(x)))

prob_clustDP<-vector("numeric", length = 112)
p_ks <-table(K_chain)
prob_clustDP[as.numeric(c(names(p_ks)))]<- p_ks/sum(p_ks)
prob_clustDP[-as.numeric(c(names(p_ks)))]<- 0
plot(1:112,prob_clustDP, col="red", type="l")
#save(prob_clustDP, file="DP_post_clust.Rdata")




##################################Covariance matrix######################################################################
#### Add cluster labels 
Colnames_Y_clust<- merge(Colnames_Y, Cluster_DP_2, by ="CODE_CBNA")

### Covariance matrix for the mean 
#pdf(file = "Plots/Correlation_matrix_DP.pdf", width= 8.27, height = 9.69)
# MDP= (fit_gjamDP$parameters$corMu + fit_gjamDP_2$parameters$corMu)/2
# Colnames_Y_clust$Sp_name_DP <- paste(Colnames_Y_clust$ClustDP, Colnames_Y_clust$species,sep="_")
# rownames(MDP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP"],"other")
# colnames(MDP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP"],"other")
# colors_vir=viridis(length(unique(Colnames_Y_clust$ClustDP))+1, option = "magma")
# LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_DP),"ClustDP"],length(unique(Colnames_Y_clust$ClustDP))+1), function(x) colors_vir[x])
# cols = colorRampPalette(c("dark blue","white","red"))
# col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
#                            "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
#                            "#67001F", "#B2182B", "#D6604D", "#F4A582"))
# 
# 
# 
# corrplot(MDP, diag = FALSE, order = "hclust", tl.cex = 0.45,tl.srt=45, method = "color", tl.col=LabelCol,col=cols(200),
#          type = "full", title= "Correlation for the DP model (original)", mar=c(0,0,1,0))
# 
# corrplot(MDP, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45,  tl.col=LabelCol,
#          method = "color",col=cols(200), type = "full",title= "Correlation for the DP model (DP groups)", mar=c(0,0,1,0))
# 
# #dev.off()


pdf(file = "Plots/Correlation_matrix_DP.pdf", width= 8.27, height = 9.69)
MDP= (fit_gjamDP$parameters$corMu + fit_gjamDP_2$parameters$corMu)/2
cor_matrix_DP = MDP[1:111,1:111]
for (i in 1:111){
  if (Colnames_Y_clust$ClustDP[i] > 9){
    Colnames_Y_clust$Sp_name_DP[i] <- paste("Group",Colnames_Y_clust$ClustDP[i], Colnames_Y_clust$species[i],sep=":")
  }
  else{
    Colnames_Y_clust$Sp_name_DP[i] <- paste(" Group",Colnames_Y_clust$ClustDP[i], Colnames_Y_clust$species[i],sep=":")
  }
}
#Colnames_Y_clust$Sp_name_DP1 <- paste("Group",Colnames_Y_clust$ClustDP1, Colnames_Y_clust$species,sep=":")
rownames(cor_matrix_DP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP"])
colnames(cor_matrix_DP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP"])
colors_vir=viridis(length(unique(Colnames_Y_clust$ClustDP)), option = "magma")
LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_DP),"ClustDP"],length(unique(Colnames_Y_clust$ClustDP))), function(x) colors_vir[x])
cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))


corrplot(cor_matrix_DP, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45,  tl.col=LabelCol,
         method = "color",col=cols(200), type = "lower",title= "", mar=c(0,0,1,0))

dev.off()

######## Precision matrix ################################################################################
library(qgraph)
library(matlib)
#inv_cor = solve(MDP)
inv_cor = inv(MDP)

PC_mat <- wi2net(inv_cor)

### Get precision matrix from the Correlation
PC_mat[abs(PC_mat)> 0.09] 
## values  are  super  small 
########################################################################################
#### Final Table
form<-c(formula)
Fin_tab_PY<-as.data.frame(matrix(NA,nrow=10,ncol=6))
names(Fin_tab_PY)<- c("Parameter","PY","r", "iter", "burn","formula")
Fin_tab_PY$iter<- fit_gjamPY1$modelList$ng
Fin_tab_PY$burn<- fit_gjamPY1$modelList$burnin
Fin_tab_PY$r<-fit_gjamPY1$modelList$reductList$r
Fin_tab_PY$formula<-as.character(form)
Fin_tab_PY[1,1]<- "DIC"
Fin_tab_PY[1,2]<- fit_gjamPY1$fit$DIC/10000
Fin_tab_PY[2,1]<- "mean AUC"
Fin_tab_PY[2,2]<- mean(gjamPY$AUC_out)
Fin_tab_PY[3,1]<- "AUC in"
Fin_tab_PY[3,2]<- mean(gjamPY$AUC_in)
Fin_tab_PY[4,1]<- "AR dist VI loss 1"
Fin_tab_PY[4,2]<- SW_fin_table_PY1$AR_dist
Fin_tab_PY[5,1]<- "VI dist VI loss"
Fin_tab_PY[5,2]<- SW_fin_table_PY1$VI_dist
Fin_tab_PY[6,1]<- "AR dist VI loss 2"
Fin_tab_PY[6,2]<- GRE_fin_table_PY1$AR_dist
Fin_tab_PY[7,1]<- "VI dist VI loss 2"
Fin_tab_PY[7,2]<- GRE_fin_table_PY1$VI_dist
Fin_tab_PY[8,1]<- "mean K"
Fin_tab_PY[8,2]<- mean(apply(fit_PY_comb$chains$kgibbs,1,function(x) length(unique(x)))[fit_PY_comb$modelList$burnin:fit_PY_comb$modelList$ng])
Fin_tab_PY[9,1]<- "K VI"
Fin_tab_PY[9,2]<- SW_fin_table_PY1$K
Fin_tab_PY[10,1]<- "K VI2"
Fin_tab_PY[10,2]<- GRE_fin_table_PY1$K
Fin_tab_PY[,2]<- round(Fin_tab_PY[,2], 3)
#save(Fin_tab_PY, file = paste0(folderpath,"Fin_tabPY_r5.Rdata"))
#library("xlsx")
# Write the first data set in a new workbook
#write.xlsx(Fin_all, file = "Final_table.xlsx")
#########################################################################################