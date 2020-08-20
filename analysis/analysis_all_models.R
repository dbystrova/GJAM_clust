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
library(xlsx)
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
K_prior=16
r_reduct = 5

folderpath="PCA_analysis/r5/"
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

##### Models description

##Calibration table
##### ALpha posterior

load("PCA_analysis/alpha/aDP16.Rdata")
load("PCA_analysis/alpha/aDP1_8.Rdata")
load("PCA_analysis/alpha/aDP1_56.Rdata")
#"#3E4A89FF" 
DF_alpha <- rbind(aDP1_8,aDP16,aDP1_56)
DF_alpha$K_pr = factor(DF_alpha$K_pr, levels=c('K = 8','K = 16','K = 56'))
#pdf(file = "Plots/DP1_alpha_posterior.pdf", width=6, height =3)
 aDP1<-   DF_alpha%>%
   ggplot(aes(x=Probability, fill=Distribution)) +scale_fill_manual(values=c("#26828EFF", "#453781FF"))+
   geom_density(adjust = 2, alpha=0.5)+
   #ggtitle(c) +
   xlab(TeX(sprintf('$\\alpha$')))+
   ylab("Probability")+
   theme_bw() +facet_grid(.~K_pr, scale="free")+
   theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "none", plot.title = element_text(hjust = 0.5))
# aDP1
 
 #plot(aDP1)
 #dev.off()

######################################## Convergence

### Effective Size
### beta coefficients
df_beta_DP1<- load_object("PCA_analysis/r5/DP1_analysis/Conv_beta_DP1.Rdata")
df_beta_PY<- load_object("PCA_analysis/r5/PY_analysis/Conv_beta_PY1.Rdata")
df_beta_DP<- load_object("PCA_analysis/r5/DP_analysis/Conv_beta_DP.Rdata")
### sigma coefficients
df_sigma_DP1<- load_object("PCA_analysis/r5/DP1_analysis/Conv_sigma_DP1.Rdata")
df_sigma_PY<- load_object("PCA_analysis/r5/PY_analysis/Conv_sigma_PY1.Rdata")
df_sigma_DP<- load_object("PCA_analysis/r5/DP_analysis/Conv_sigma_DP.Rdata")
beta_ess<- rbind(df_beta_DP1 %>% mutate(Model = "DP"),
                   df_beta_PY %>%  mutate(Model = "DP1"), 
                   df_beta_DP %>% mutate(Model = "PY")
)
beta_ess$parameter <- "Regression coefficients"
sigma_ess<- rbind(df_sigma_DP %>% mutate(Model = "DP"),
                    df_sigma_DP1 %>%  mutate(Model = "DP1"), 
                    df_sigma_PY %>% mutate(Model = "PY"))
sigma_ess$parameter<- "Correlation coefficients"
## full dataframe
ess_df<- rbind(beta_ess,sigma_ess)
#pdf
#pdf(file = "Plots/ESS_all.pdf", width= 6.7, height =3.5)
ess_plot<- ggplot(ess_df,aes(ES,color=Model, fill = Model))+ geom_histogram( position="identity",alpha=0.3, bins=50) +
  xlab("Effective sample size")+
  ylab("Counts") + facet_wrap(.~parameter,scales = "free")+
  scale_color_manual(guide="legend",values=c("#B63679FF","#31688EFF" ,"#35B779FF"), labels=c("DP ", expression(DP[c]), expression(PY[c])))+
  scale_fill_manual(guide="legend",values=c("#B63679FF","#31688EFF" ,"#35B779FF"), labels=c("DP ", expression(DP[c]), expression(PY[c])))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5))
 ess_plot  
plot(ess_plot)
#dev.off()

# 
# #pdf
# pdf(file = "Plots/Sigma_ESS.pdf", width= 3.5, height = 2.5)
# p_sigma_ess<- rbind(df_sigma_DP %>% mutate(Model = "DP"),
#       df_sigma_DP1 %>%  mutate(Model = "DP1"), 
#       df_sigma_PY %>% mutate(Model = "PY")
# ) %>% 
#   ggplot(aes(ES, color = Model, fill = Model))+ geom_histogram( position="identity",alpha=0.3, bins=50) +
#   xlab("Effective size")+
#   ylab("Counts") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "none", plot.title = element_text(hjust = 0.5))
# plot(p_sigma_ess)
# dev.off()
# 
# #pdf
# pdf(file = "Plots/Beta_ESS.pdf", width= 3.5, height =2.5)
# p_beta_ess<- rbind(df_beta_DP1 %>% mutate(Model = "DP"),
#       df_beta_PY %>%  mutate(Model = "DP1"), 
#       df_beta_DP %>% mutate(Model = "PY")
# ) %>% 
#   ggplot(aes(ES, color = Model, fill = Model))+ geom_histogram( position="identity",alpha=0.3, bins=50) +
#   ylab("")+
#   xlab("Effective size") + theme_bw()+
#   theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "right", plot.title = element_text(hjust = 0.5))
# plot(p_beta_ess)
# dev.off()
### Gelman-Rubin diagnostics
## sigma coefficients
GR_sigma_DP<- load_object("PCA_analysis/r5/DP_analysis/GR_value_sigma_DP.Rdata")
GR_sigma_DP1<- load_object("PCA_analysis/r5/DP1_analysis/GR_value_sigma_DP1.Rdata")
GR_sigma_PY<- load_object("PCA_analysis/r5/PY_analysis/GR_value_sigma_PY1.Rdata")

#beta coefficients
GR_beta_DP<- load_object("PCA_analysis/r5/DP_analysis/GR_value_beta_DP.Rdata")
GR_beta_DP1<- load_object("PCA_analysis/r5/DP1_analysis/GR_value_beta_DP1.Rdata")
GR_beta_PY<- load_object("PCA_analysis/r5/PY_analysis/GR_value_beta_PY1.Rdata")

GR_sigma<- tibble( DP = GR_sigma_DP$psrf[,1],
                   DP1 = GR_sigma_DP1$psrf[,1],
                   PY = GR_sigma_PY$psrf[,1]) %>% gather(Model, R_hat, DP:PY)
GR_sigma$parameter <-  "Correlation coefficients"


GR_beta<- tibble( DP = GR_beta_DP$psrf[,1],
        DP1 = GR_beta_DP1$psrf[,1],
        PY = GR_beta_PY$psrf[,1]) %>% gather(Model, R_hat, DP:PY)
GR_beta$parameter <-  "Regression coefficients"


GR<- rbind(GR_beta,GR_sigma)
GR<-GR[GR$R_hat<=1.1,]
GR$Model <- as.factor(GR$Model)
levels(GR$Model) <- c("DP",  expression(DP[c]), expression(PY[c]))
#GR$Model <- factor(GR$Model ,levels =c("DP ", expression(DP[c]), expression(PY[c])) )
#pdf(file = "Plots/Rhat_all.pdf", width= 6.4, height =5)
gr_plot<- ggplot(GR,aes(R_hat, color = Model, fill = Model))+
  geom_histogram(position="identity",alpha=0.3, bins=50)+
  xlab("Potential scale reduction factor")+
  scale_color_manual(guide="legend",values=c("#B63679FF","#31688EFF" ,"#35B779FF"))+
  scale_fill_manual(guide="legend",values=c("#B63679FF","#31688EFF" ,"#35B779FF"))+
  ylab("Counts") + facet_grid(Model~parameter,scales = "free",space="free", labeller = labeller(Model = label_parsed))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "none", plot.title = element_text(hjust = 0.5))
gr_plot  
plot(gr_plot)
#dev.off()

#pdf(file = "Plots/Rhat_sigma.pdf", width= 8.27, height = 9.69)
# tibble( DP = GR_sigma_DP$psrf[,1],
#         DP1 = GR_sigma_DP1$psrf[,1],
#         PY = GR_sigma_PY$psrf[,1]) %>% gather(Model, R_hat, DP:PY)%>%
#   ggplot(aes(R_hat, color = Model, fill = Model))+
#   geom_histogram( position="identity",alpha=0.3, bins=50)+ facet_wrap(~ Model, ncol=1)+
#   xlim(min(hm),1.1)
# #dev.off()
# 
# 
# GR_beta_DP<- load_object("PCA_analysis/r5/DP_analysis/GR_value_beta_DP.Rdata")
# GR_beta_DP1<- load_object("PCA_analysis/r5/DP1_analysis/GR_value_beta_DP1.Rdata")
# GR_beta_PY<- load_object("PCA_analysis/r5/PY_analysis/GR_value_beta_PY1.Rdata")
# 
# 
# #pdf(file = "Plots/Rhat_beta.pdf", width= 8.27, height = 9.69)
# tibble( DP = GR_beta_DP$psrf[,1],
#         DP1 = GR_beta_DP1$psrf[,1],
#         PY = GR_beta_PY$psrf[,1]) %>% gather(Model, R_hat, DP:PY)%>%
#   ggplot(aes(R_hat, color = Model, fill = Model))+ geom_histogram( position="identity",alpha=0.3, bins=50)+ facet_wrap(~ Model, ncol=1)
# #dev.off()


### Results



### Prediction
DP_pred<- load_object("PCA_analysis/r5/DP_analysis/DP_prediction.Rdata")
DP1_pred<- load_object("PCA_analysis/r5/DP1_analysis/DP1_prediction.Rdata")
PY_pred<- load_object("PCA_analysis/r5/PY_analysis/PY1_prediction.Rdata")

AUC_data<- tibble(DP=DP_pred$AUC_out, DP1=DP1_pred$AUC_out, PY= PY_pred$AUC_out)
AUC_data_in<-  tibble(DP=DP_pred$AUC_in, DP1=DP1_pred$AUC_in, PY= PY_pred$AUC_in)
AUC_data$species<- colnames(Ydata_test)[1:ncol(Ydata_test)]
AUC_fin<- melt(AUC_data)

p2<-ggplot(data=AUC_fin)+geom_boxplot(aes(y=as.numeric(value),x=as.factor(variable),fill=as.factor(variable)))+
  scale_y_continuous(name="AUC")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","PY1"))+xlab("Models")+ theme_bw() 
p2

AUC_fin_in_table<- as.data.frame(t(apply(AUC_data_in,2,mean)))
AUC_fin_table<- as.data.frame(t(apply(AUC_data[,1:3],2,mean)))
WAUC_fin_table<- as.data.frame(t(apply(AUC_data[,1:3],2,function(x) sum(x*p_w))))
names(AUC_fin_table)<- c("DP","DP1","PY")
kable(cbind(data.frame( Measure = rbind("AUC")), rbind(AUC_fin_table)), format="pandoc", caption= "Prediction")
Pred_tab <-as.data.frame( rbind(cbind(data.frame( Measure = rbind("AUC_out")), rbind(AUC_fin_table)),
            cbind(data.frame( Measure = rbind("AUC_in")), rbind(AUC_fin_in_table))))
 
##############################################################################################################
### Clustering: K, distances
### Trace chains 
DP1_k_trace<- load_object("PCA_analysis/r5/DP1_analysis/DP1_k_chains.Rdata")
DP_k_trace<- load_object("PCA_analysis/r5/DP_analysis/DP_k_chains.Rdata")
PY_k_trace<- load_object("PCA_analysis/r5/PY_analysis/PY1_k_chains.Rdata")

df<-cbind(DP_k_trace,DP1_k_trace[,2:3],PY_k_trace[,2:3] )
df%>% gather(Model, trace, DP_1:PY1_2)%>%
  ggplot(aes(x=it,y=trace,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))

##  Clustering distance
DP_clust_tab_SW <- load_object("PCA_analysis/r5/DP_analysis/SW_tab_DP.Rdata")
DP1_clust_tab_SW <- load_object("PCA_analysis/r5/DP1_analysis/SW_tab_DP1.Rdata")
PY_clust_tab_SW <- load_object("PCA_analysis/r5/PY_analysis/SW_tab_PY1.Rdata")

SW_tab<- rbind(DP_clust_tab_SW,DP1_clust_tab_SW,PY_clust_tab_SW)
#write.csv(SW_tab, file = "SW_tab.csv")


DP1_clust_tab_GRE <- load_object("PCA_analysis/r5/DP1_analysis/GRE_tab_DP1.Rdata")
DP1_clust_tab_GRE$Model<- "DP1"
DP_clust_tab_GRE <- load_object("PCA_analysis/r5/DP_analysis/GRE_tab_DP.Rdata")
PY_clust_tab_GRE <- load_object("PCA_analysis/r5/PY_analysis/GRE_tab_PY1.Rdata")

GRE_tab<- rbind(DP1_clust_tab_GRE,DP_clust_tab_GRE,PY_clust_tab_GRE)
#write.xlsx(GRE_tab, file = "GRE_tab.xlsx")

############################################ Clustring  distances


############################################ Clusters
load("PCA_analysis/r5/DP1_analysis/Cluster_DP1_1.Rdata")
load("PCA_analysis/r5/DP1_analysis/Cluster_DP1_2.Rdata")
load("PCA_analysis/r5/PY_analysis/Cluster_PY1_1.Rdata")
load("PCA_analysis/r5/PY_analysis/Cluster_PY1_2.Rdata")
load("PCA_analysis/r5/DP_analysis/Cluster_DP_2.Rdata")
load("PCA_analysis/r5/DP_analysis/Cluster_DP_1.Rdata")

Clusters_all_1<- merge(Cluster_DP_1,Cluster_DP1_1, by = c("CODE_CBNA"))
Clusters_all_1<- merge(Clusters_all_1, Cluster_PY1_1,  by = c("CODE_CBNA"))
Clusters_all_1<- merge(Clusters_all_1, True_clust[, c("CODE_CBNA", "PFG", "K_n")], by = c("CODE_CBNA"))
#save(Clusters_all_1, file =  paste0("PCA_analysis/r5/Clusters/Clusters_all_1.Rdata"))


Clusters_all_2<- merge(Cluster_DP_2,Cluster_DP1_2, by = c("CODE_CBNA"))
Clusters_all_2<- merge(Clusters_all_2, Cluster_PY1_2,  by = c("CODE_CBNA"))
Clusters_all_2<- merge(Clusters_all_2, True_clust[, c("CODE_CBNA", "PFG", "K_n")], by = c("CODE_CBNA"))
#save(Clusters_all_2, file =  paste0("PCA_analysis/r5/Clusters/Clusters_all_2.Rdata"))


#Clusters_all_1_old <- load_object("PCA_analysis/r5/Clusters/Clusters_all_1_old.Rdata")
#arandi(Clusters_all_1_old$ClustPY1, Clusters_all_1$ClustPY1)

#Clusters_all_2_old <- load_object("PCA_analysis/r5/Clusters/Clusters_all_2_old.Rdata")
#arandi(Clusters_all_2_old$ClustDP, Clusters_all_2$ClustDP)

#arandi(Clusters_all_1$ClustPY1, Clusters_all_1$ClustDP1)


#Clusters_all_2$DP <- as.numeric(Clusters_all_2$CODE_CBNA)
#Clusters_all_1$CBN <- as.numeric(Clusters_all_1$CODE_CBNA)

#M= Clusters_all_1[order(Clusters_all_1$CBN),]
#B= Clusters_all_2[order(Clusters_all_2$CBN),]

#arandi(PY1_clust2_full$VI_est[[3]], B$ClustPY1)

#### Pairwise distances

########### Cluster distances with random 
#load("PCA_analysis/r5/Cluster_DP1_1.Rdata")
#load("~/Documents/GitHub/GJAM_clust/PCA_analysis/r5/DP1_analysis/Cluster_DP1_1.Rdata")
Clusters_all_2$CBN <- as.numeric(Clusters_all_2$CODE_CBNA)

Clusters_all_2_sorted= Clusters_all_2[order(Clusters_all_2$CBN),]
Random_cluster <-  sample(1:16, size = 111, replace = TRUE)
Random_cluster_w <-  sample(1:16, size = 111, replace = TRUE, prob=table(Clusters_all_2_sorted$K_n)/111)

#arandi(Cluster_PY1_2$ClustPY1,Clusters_all_2_sorted$ClustPY1 )
## Two heatmap plots
M_all<- data.frame()
Clusters_by_Rand<-cbind(Clusters_all_2_sorted, tibble(RU=Random_cluster), tibble(RW= Random_cluster_w))
Clusters_by_Rand_mat<- Clusters_by_Rand[,c("ClustDP","ClustDP1", "ClustPY1", "PFG", "RU", "RW")]
Mat_dist<- matrix(NA, nrow=6, ncol=6)
for(j in 1:6){
  for (k in 1:6){
    Mat_dist[j,k]= arandi(Clusters_by_Rand_mat[,j],Clusters_by_Rand_mat[,k])
  }
}

Mat_dist2<- matrix(NA, nrow=3, ncol=3)
for(j in 1:3){
  for (k in 1:3){
    Mat_dist2[j,k]= arandi(Clusters_by_Rand_mat[,j],Clusters_by_Rand_mat[,k+3])
  }
}

M_print_1<-  Dist_matrix[1:4,1:4]
M_print_1$Model <- colnames(Clusters_by_Rand_mat)[1:4]
M_rand_1<- melt(M_print_1, id=c("Model"))
M_rand_1$Model<- factor(M_rand_1$Model,levels = colnames(M_print_1))
Dist_clust.heatmap <- ggplot(data = M_rand_1, mapping = aes(x = variable,
                                                         y = Model,
                                                         fill = value, color='white')) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2)), size=5)+
  scale_fill_gradient(low = 'white', high = 'black', space = 'Lab') +
  xlab("")+
  ylab("")+
  theme_bw() +
  scale_y_discrete(breaks=c("ClustDP", "ClustDP1", "ClustPY1","PFG"),labels=c("DP", expression(DP[c]), expression(PY[c]),expression(PFG)))+
  scale_x_discrete(breaks=c("ClustDP", "ClustDP1", "ClustPY1","PFG"),labels=c("DP", expression(DP[c]),expression(PY[c]),expression(PFG)))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5))+
  guides(color=FALSE)

Dist_clust.heatmap


Dist_matrix_random<- as.data.frame(round(Mat_dist2,3))
colnames(Dist_matrix_random)<- colnames(Clusters_by_Rand_mat[4:6])
rownames(Dist_matrix_random)<- colnames(Clusters_by_Rand_mat[,1:3])
Dist_matrix_random$Clustering <- colnames(Clusters_by_Rand_mat[,1:3])
Dist_matrix_random.long <- Dist_matrix_random%>% gather(Reference, Distance, PFG:RW)


Dist_clust_random.heatmap <- ggplot(data = Dist_matrix_random.long, mapping = aes(x = Clustering,
                                                            y = Reference,
                                                            fill = Distance, color='white')) +
  geom_tile(aes(fill = Distance)) +
  geom_text(aes(label = round(Distance, 2)), size=5)+
  scale_fill_gradient(low = 'white', high = 'black', space = 'Lab', limits=c(0,1)) +
  xlab("")+
  ylab("")+
  theme_bw() +
  scale_x_discrete(breaks=c("ClustDP", "ClustDP1", "ClustPY1"),labels=c("DP", expression(DP[c]),expression(PY[c])))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5))+
  guides(color=FALSE)

Dist_clust_random.heatmap

pdf(file = "Plots/Final_cluster_distances_random.pdf", width=5, height =3)
print(Dist_clust_random.heatmap)
dev.off()



# M_print<- as.data.frame(round(Mat_dist,3))
# colnames(M_print)<- colnames(Clusters_by_Rand_mat)
# M_print_1<-  M_print[1:4,1:4]
# M_print_1$Model <- colnames(Clusters_by_Rand_mat)[1:4]
# M_rand_1<- melt(M_print_1, id=c("Model"))
# M_rand_1$Model<- factor(M_rand_1$Model,levels = colnames(M_print_1))
# Dist_clust.heatmap <- ggplot(data = M_rand_1, mapping = aes(x = variable,
#                                                          y = Model,
#                                                          fill = value, color='white')) +
#   geom_tile(aes(fill = value)) +
#   geom_text(aes(label = round(value, 2)), size=5)+
#   scale_fill_gradient(low = 'white', high = 'black', space = 'Lab') +
#   xlab("")+
#   ylab("")+
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5))+
#   guides(color=FALSE)
#
# Dist_clust.heatmap

pdf(file = "Plots/Final_cluster_distances_between_models.pdf", width=5, height =3)
print(Dist_clust.heatmap)
dev.off()



M_print_2<-  M_print[3:6,3:6]
M_print_1$Model <- colnames(Clusters_by_Rand_mat)[1:4]
M_rand_1<- melt(M_print_1, id=c("Model"))
M_rand_1$Model<- factor(M_rand_1$Model,levels = colnames(M_print_1))
Dist_clust.heatmap <- ggplot(data = M_rand_1, mapping = aes(x = variable,
                                                            y = Model,
                                                            fill = value, color='white')) + 
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 2)), size=5)+
  scale_fill_gradient(low = 'white', high = 'black', space = 'Lab')+
  xlab(label = "Distance between the optimal clusters and PFG partition ")+
  ylab("")+ 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5))

Dist_clust.heatmap

pdf(file = "Final_cluster_distances_AR.pdf", width=4, height =3)
print(Dist_clust.heatmap)
dev.off()
########################################################################################

arandi(Clusters_all_2_sorted$ClustDP1, Clusters_all_2_sorted$ClustPY1)
arandi(Clusters_all_1_sorted$ClustPY1, Clusters_all_2_sorted$ClustPY1)


arandi(Clusters_all_2_sorted$ClustDP1, DP1_clust2_wp$VI_est)
arandi(Clusters_all_2_sorted$ClustPY1, PY1_clust2_wp$VI_est)
arandi(Cluster_models_2$ClustPY1, PY1_clust2_wp$VI_est[[1]])


load("~/Documents/GitHub/GJAM_clust/PCA_analysis/r5/Clusters_all_1.Rdata")
Clusters_all_1$CBN <- as.numeric(Clusters_all_1$CODE_CBNA)
Clusters_all_1_sorted= Clusters_all_1[order(Clusters_all_1$CBN),]

arandi(Clusters_all_1_sorted$ClustDP1, DP1_clust2_wp$VI_est)
arandi(Clusters_all_1_sorted$ClustPY1, PY1_clust2_wp$VI_est[[1]])


load("~/Documents/GitHub/GJAM_clust/PCA_analysis/r5/Clusters_all_2.Rdata")
Clusters_all_2$CBN <- as.numeric(Clusters_all_2$CODE_CBNA)
Clusters_all_2_sorted= Clusters_all_2[order(Clusters_all_2$CBN),]



########################################################################################
########################################################################################
########################################################################################

#### Prior/Posterior distrbution - clustering  
# 
# DP_post_dist<- load_object("DP_post_clust.Rdata")
# DP1_post_dist16<- load_object("DP1_post_clust16.Rdata")
# PY_post_dist16<- load_object("PY_post_clust16.Rdata")
# 
# Dist_K_post_16<- tibble( k = 1:112,
#                          type= rep("posterior",112),
#                      DP = DP_post_dist,
#                      DP1 = DP1_post_dist16,
#                      PY = PY_post_dist16)
# 
# 
# x_vec<- 1:112
# load("IJulia_part/Cnk_mat_112_H05.Rdata")
# pks_05<- sapply(x_vec, PY_prior,H=112, n=112, alpha=0.47, sigma=0.5,Cnk_mat=Cnk_112_112_H05)
# plot(x_vec, pks_05)
# 
# PY_prior_16_05<- pks_05
# 
# load('IJulia_part/DPM_original_prior.Rdata')
# load('IJulia_part/DPM_prior_16.Rdata')
# 
# Dist_K_prior_16<- tibble( k = 1:112, 
#                           type= rep("prio",112),
#                          DP = DPM_original_prior,
#                          DP1 = DPM_prior_16,
#                          PY = PY_prior_16_05)
# Dist_16<- rbind(Dist_K_post_16, Dist_K_prior_16)
# pdf(file = "Plots/Prior_all_16.pdf", width= 8.27, height = 9.69)
# a <- Dist_16%>% gather(Model, p_k, DP:PY)%>%
#   ggplot(aes(x=k, y= p_k, color= Model, linetype=type)) + geom_line() +scale_x_continuous(breaks = seq(6, 112, by = 5))+
#   ggtitle("Prior and posterior distribution for E[K_n]=16 for different models")+ theme_bw() 
# a
# plot(a)
# dev.off()


############### Prior specification 16  ############################################################
load('IJulia_part/DPM_original_prior.Rdata')

load('IJulia_part/DPM_prior_16.Rdata')


x_vec<- 1:112
load("IJulia_part/Cnk_mat_112_H05.Rdata")
pks_05<- sapply(x_vec, PY_prior,H=112, n=112, alpha=0.47, sigma=0.5,Cnk_mat=Cnk_112_112_H05)
plot(x_vec, pks_05)
load("IJulia_part/Cnk_mat_112_H025.Rdata")
pks_025<- sapply(x_vec, PY_prior,H=112, n=112, alpha=2.574, sigma=0.25,Cnk_mat=Cnk_112_112_H025)
plot(x_vec, pks_025)
load("IJulia_part/Cnk_mat_112_H08.Rdata")
pks_08<- sapply(x_vec, PY_prior,H=112, n=112, alpha=2.574, sigma=0.8,Cnk_mat=Cnk_112_112_H08)
plot(x_vec, pks_08)
exp08<-sum(x_vec*pks_08)
exp08



exp025<-sum(x_vec*pks_025)
PY_prior_16_025<- pks_025
PY_prior_16_05<- pks_05

prior_DP_16_fixed <- load_object("IJulia_part/DPM_prior16_fixed.Rdata")

Prior_spec_16<- tibble( k = 1:112,
                          DP = DPM_original_prior,
                          DP_1 = prior_DP_16_fixed,
                          DP1 = DPM_prior_16,
                          PY_2 = PY_prior_16_025,
                          PY_1= PY_prior_16_05)
pdf(file = "Plots/Prior_16_specification.pdf", width=6.5, height = 4)
a <- Prior_spec_16%>% gather(Model, p_k, DP:PY_1)%>%
  ggplot(aes(x=k, y= p_k, color= Model)) + geom_line(size=0.9) +
  scale_color_manual(values=c("#B63679FF", "#FB8861FF","#31688EFF" ,"#35B779FF","#7AD151FF"), labels = c("DP(1)", "DP(2)",expression(DP[c]), expression(PY[c](1)), expression(PY[c](2)) )) + 
  scale_x_continuous(limits=c(0,70))+
  ylab("Probability")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5),legend.text.align = 0)
a
plot(a)
dev.off()




#### Prior/Posterior distrbution - clustering  

DP_post_dist<- load_object("DP_post_clust.Rdata")
DP1_post_dist16<- load_object("DP1_post_clust16.Rdata")
PY_post_dist16<- load_object("PY_post_clust16.Rdata")

Dist_K_post_16<- tibble( k = 1:112,
                         type= rep("posterior",112),
                         DP = DP_post_dist,
                         DP1 = DP1_post_dist16,
                         PY = PY_post_dist16)


x_vec<- 1:112
load("IJulia_part/Cnk_mat_112_H05.Rdata")
pks<- sapply(x_vec, PY_prior,H=112, n=112, alpha=0.47, sigma=0.5,Cnk_mat=Cnk_112_112_H05)
plot(x_vec, pks)

exp<-sum(x_vec*pks)
exp
PY_prior_16<- pks

load('IJulia_part/DPM_original_prior.Rdata')
load('IJulia_part/DPM_prior_16.Rdata')

breaks= c(0,10,20,30,40,50,60)
Dist_K_prior_16<- tibble( k = 1:112, 
                          type= rep("prior",112),
                          DP = DPM_original_prior,
                          DP1 = DPM_prior_16,
                          PY = PY_prior_16)
Dist_16<- rbind(Dist_K_post_16, Dist_K_prior_16)
pdf(file = "Plots/Prior_all_16.pdf", width= 6.5, height = 4)
a <- Dist_16%>% gather(Model, p_k, DP:PY)%>%
  ggplot(aes(x=k, y= p_k, color= Model, linetype=type)) + geom_line(size=0.9) +
  scale_color_manual(values=c("#B63679FF","#31688EFF" ,"#35B779FF"),labels= c("DP", expression(DP[c]), expression(PY[c]))) + 
  scale_x_continuous(breaks = breaks,minor_breaks = c(16), limits = c(0,70))+
  #scale_x_continuous(breaks = seq(6, 112, by = 5))+
  ylab("Probability")+
  #ggtitle("Prior and posterior distribution forE[K_n]=16 for different models")+ 
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5),legend.text.align = 0)+
  guides(color = guide_legend(order = 1),linetype = guide_legend(order = 2))
a
plot(a)
dev.off()
#################################################################################### Wrong prior
#### Prior/Posterior distrbution - clustering  

DP_post_dist<- load_object("DP_post_clust.Rdata")
DP1_post_dist56<- load_object("DP1_post_clust_56.Rdata")
PY_post_dist56<- load_object("prob_clust_PY_56.Rdata")

Dist_K_post_56<- tibble( k = 1:112,
                         type= rep("posterior",112),
                         DP = DP_post_dist,
                         DP1 = DP1_post_dist56,
                         PY = PY_post_dist56)

x_vec<- 1:112
load("IJulia_part/Cnk_mat_112_H08.Rdata")
pks<- sapply(x_vec, PY_prior,H=112, n=112, alpha=7.7, sigma=0.8,Cnk_mat=Cnk_112_112_H08)
plot(x_vec, pks)

exp<-sum(x_vec*pks)
exp
PY_prior_56<- pks

load('IJulia_part/DPM_original_prior.Rdata')
DPM_prior_56_ <-   load_object('IJulia_part/DPM_prior56.Rdata')

Dist_K_prior_56<- tibble( k = 1:112, 
                          type= rep("prior",112),
                          DP = DPM_original_prior,
                          DP1 = DPM_prior_56_,
                          PY = PY_prior_56)
Dist_56<- rbind(Dist_K_post_56, Dist_K_prior_56)

Df_56 <- Dist_56%>% gather(Model, p_k, DP:PY)
#Df_56$prior <- TeX(sprintf('Prior and posterior distribution for mis-calibrated DP1 $\\alpha$'))
Df_56$prior <- "Mis-callibrated (K=56)"


#pdf(file = "Plots/Prior_all_56.pdf", width= 8.27, height = 9.69)
a <- Dist_56%>% gather(Model, p_k, DP:PY)%>%
  ggplot(aes(x=k, y= p_k, color= Model, linetype=type)) + geom_line() +scale_x_continuous(breaks = seq(6, 112, by = 5))+
  ggtitle("Prior and posterior distribution forE[K_n]=56 for different models")+ theme_bw() 
a
plot(a)
dev.off()





#### Prior/Posterior distrbution - clustering  

DP_post_dist<- load_object("DP_post_clust.Rdata")
DP1_post_dist8<- load_object("/Users/dariabystrova/Documents/GitHub/GJAM_clust/PCA_analysis/r_wp/wp_8/DP1/DP1_post_clust_8.Rdata")
PY_post_dist8<- load_object("/Users/dariabystrova/Documents/GitHub/GJAM_clust/PCA_analysis/r_wp/wp_8/PY/PY_post_clust_8.Rdata")

Dist_K_post_8<- tibble( k = 1:112,
                         type= rep("posterior",112),
                         DP = DP_post_dist,
                         DP1 = DP1_post_dist8,
                         PY = PY_post_dist8)


x_vec<- 1:112
load("IJulia_part/Cnk_mat_112_H025.Rdata")
pks<- sapply(x_vec, PY_prior,H=112, n=112, alpha=0.64, sigma=0.25,Cnk_mat=Cnk_112_112_H025)
plot(x_vec, pks)

exp<-sum(x_vec*pks)
exp
PY_prior_8<- pks

load('IJulia_part/DPM_original_prior.Rdata')
load('IJulia_part/DPM_prior56.Rdata')

DPM_prior_8 <-load_object("/Users/dariabystrova/Documents/GitHub/GJAM_clust/PCA_analysis/r_wp/wp_8/DP1/DPM_prior_8.Rdata")

Dist_K_prior_8<- tibble( k = 1:112, 
                          type= rep("prior",112),
                          DP = DPM_original_prior,
                          DP1 = DPM_prior_8,
                          PY = PY_prior_8)
Dist_8<- rbind(Dist_K_post_8, Dist_K_prior_8)
#pdf(file = "Plots/Prior_all_8.pdf", width= 8.27, height = 9.69)
a <- Dist_8%>% gather(Model, p_k, DP:PY)%>%
  ggplot(aes(x=k, y= p_k, color= Model, linetype=type)) + geom_line() +scale_x_continuous(breaks = seq(6, 112, by = 5))+
  ggtitle("Prior and posterior distribution forE[K_n]=8 for different models")+ theme_bw() 
plot(a)
#dev.off()

Df_8 <- Dist_8%>% gather(Model, p_k, DP:PY)
Df_8$prior <-"Mis-callibrated (K=8)"

DF_wrong_prior<- rbind(Df_8,Df_56)

pdf(file = "Plots/Wrong_prior.pdf", width= 6.5, height =6)
df_WP_plot <- DF_wrong_prior%>%
    ggplot(aes(x=k, y= p_k, color= Model, linetype=type)) + 
    geom_line(size=0.9) +
    ylab("Probability")+ 
    scale_color_manual(values=c("#B63679FF","#31688EFF" ,"#7AD151FF"),labels= c("DP", expression(DP[c]), expression(PY[c]))) + 
    scale_x_continuous(breaks = breaks,minor_breaks = c(16), limits = c(0,70))+
    facet_wrap(.~prior,scales = "free", nrow=2)+
     theme_bw()  +
     theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5),legend.text.align = 0)+
     guides(color = guide_legend(order = 1),linetype = guide_legend(order = 2))
df_WP_plot
dev.off()

##########################################################################################################




Clust2<- load_object("~/Documents/GitHub/GJAM_clust/PCA_analysis/r5/Clusters/Clusters_all_2.Rdata")

Clust2$CBN <- as.numeric(Clust2$CODE_CBNA)
Clust2_sorted= Clust2[order(Clust2$CBN),]



arandi(DP1_clust_full_56$VI_est[[3]],Clust2_sorted$ClustDP1)
arandi(PY_clust_full_56$VI_est[[3]], Clust2_sorted$ClustPY1)


WP_56_DP1<-tibble(Model="DP1",K=length(unique(DP1_clust_full_56$VI_est[[3]])),Dist_PFG = arandi(DP1_clust_full_56$VI_est[[3]], True_clust$K_n), Dist_CM = arandi(DP1_clust_full_56$VI_est[[3]], Clust2_sorted$ClustDP1))
save(WP_56_DP1, file =  paste0(folderpath,"WP_56_DP1.Rdata"))
WP_56_PY<- tibble(Model="PY",K=length(unique(PY_clust_full_56$VI_est[[3]])),Dist_PFG = arandi(PY_clust_full_56$VI_est[[3]], True_clust$K_n), Dist_CM = arandi(PY_clust_full_56$VI_est[[3]], Clust2_sorted$ClustPY1))
save(WP_56_PY, file =  paste0(folderpath,"WP_56_PY.Rdata"))


############################################################################################################
## DIC

fit_gjamDP<- load_object("PCA_analysis/r5_3/fit_gjam.Rdata")
fit_gjamDP1<- load_object("PCA_analysis/r5_4/fit_gjamDP1.Rdata")
fit_gjamPY1<- load_object("PCA_analysis/test/fit_gjamPY1.Rdata")

df_DIC_5<- tibble(r= c(5), DP= fit_gjam$fit$DIC, DP1= fit_gjamDP1$fit$DIC, PY1= fit_gjamPY1$fit$DIC)
df_DIC_10<- load_object("PCA_analysis/different_r/df_DIC_10.Rdata")
df_DIC_10$r=10
df_DIC_15<- load_object("PCA_analysis/different_r/df_DIC_15.Rdata")
df_DIC_15$r = 15
df_DIC_20<- load_object("PCA_analysis/different_r/df_DIC_20.Rdata")
df_DIC_20$r = 20

DIC_df<- rbind(df_DIC_5,df_DIC_10, df_DIC_15, df_DIC_20)


DIC_df<- rbind(df_DIC_5,df_DIC_10,df_DIC_15, df_DIC_20)
DIC_df_scaled <- DIC_df %>%gather(Models, DIC, c(DP,DP1,PY1))
DIC_df_scaled$DIC<- DIC_df_scaled$DIC/10000
DIC_plot<- DIC_df_scaled%>%ggplot(aes(x=r,y=DIC,col=Models))+geom_line(alpha=0.7)+
  xlab("number of latent factors (r)")+ylab(TeX(sprintf('DIC value $\\times 10^4$'))) +theme_bw()+
  scale_color_manual(values=c("#B63679FF","#31688EFF" ,"#35B779FF"),labels= c("DP", expression(DP[c]), expression(PY[c]))) + 
  
  
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = "right", plot.title = element_text(hjust = 0.5),legend.text.align = 0)
DIC_plot
pdf(file = "Plots/DIC.pdf", width=6, height =3)
plot(DIC_plot)
dev.off()
