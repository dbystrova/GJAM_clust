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


######################################## Convergence

### Effective Size

df_beta_DP1<- load_object("PCA_analysis/r5/Conv_beta_DP1.Rdata")
df_beta_PY<- load_object("PCA_analysis/r5/")
df_beta_DP<- load_object("PCA_analysis/r5/")

df_sigma_DP1<- load_object("PCA_analysis/r5/Conv_sigma_DP1.Rdata")
df_sigma_PY<- load_object("PCA_analysis/r5/Conv_sigma_PY1.Rdata")
df_sigma_DP<- load_object("PCA_analysis/r5/Conv_sigma_DP.Rdata")

#pdf
rbind(df_sigma_DP %>% mutate(Model = "DP"),
      df_sigma_DP1 %>%  mutate(Model = "DP1"), 
      df_sigma_PY %>% mutate(Model = "PY")
) %>% 
  ggplot(aes(ES, color = Model, fill = Model))+ geom_histogram( position="identity",alpha=0.3, bins=50) 
#dev.off()

### Gelman-Rubin diagnostics
GR_sigma_DP<- load_object("PCA_analysis/r5/GR_value_sigma_DP.Rdata")
GR_sigma_DP1<- load_object("PCA_analysis/r5/GR_value_sigma_DP1.Rdata")


GR_beta_DP<- load_object("PCA_analysis/r5/GR_value_beta_DP.Rdata")
GR_beta_DP1<- load_object("PCA_analysis/r5/GR_value_beta_DP1.Rdata")

### Results


### Clustering: K, distances
### Trace chains 
DP1_k_trace<- load_object("PCA_analysis/r5/DP1_k_chains.Rdata")
DP_k_trace<- load_object("PCA_analysis/r5/DP_k_chains.Rdata")
PY_k_trace<- load_object("PCA_analysis/r5/PY1_k_chains.Rdata")

df<-cbind(DP_k_trace,DP1_k_trace[,2:3],PY_k_trace[,2:3] )
df%>% gather(Model, trace, DP_1:PY1_2)%>%
  ggplot(aes(x=it,y=trace,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))

### Prediction
DP_pred<- load_object("PCA_analysis/r5/DP_prediction.Rdata")
DP1_pred<- load_object("PCA_analysis/r5/DP1_prediction.Rdata")
PY_pred<- load_object("PCA_analysis/r5/PY1_prediction.Rdata")

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
##  Clustering distance
DP_clust_tab_SW <- load_object("PCA_analysis/r5/SW_tab_DP.Rdata")
DP1_clust_tab_SW <- load_object("PCA_analysis/r5/SW_tab_DP1.Rdata")
PY_clust_tab_SW <- load_object("PCA_analysis/r5/SW_tab_PY1.Rdata")

SW_tab<- rbind(DP_clust_tab_SW,DP1_clust_tab_SW,PY_clust_tab_SW)
#write.xlsx(SW_tab, file = "SW_tab.xlsx")


DP1_clust_tab_GRE <- load_object("PCA_analysis/r5/GRE_tab_DP1.Rdata")
DP1_clust_tab_GRE$Model<- "DP1"
DP_clust_tab_GRE <- load_object("PCA_analysis/r5/GRE_tab_DP.Rdata")
PY_clust_tab_GRE <- load_object("PCA_analysis/r5/GRE_tab_PY1.Rdata")

GRE_tab<- rbind(DP1_clust_tab_GRE,DP_clust_tab_GRE,PY_clust_tab_GRE)
#write.xlsx(GRE_tab, file = "GRE_tab.xlsx")

############################################ Clustring  distances


############################################ Clusters

load("PCA_analysis/r5/Cluster_DP1_1.Rdata")
load("PCA_analysis/r5/Cluster_DP1_2.Rdata")
load("PCA_analysis/r5/Cluster_PY1_1.Rdata")
load("PCA_analysis/r5/Cluster_PY1_2.Rdata")
load("PCA_analysis/r5/Cluster_DP_2.Rdata")
load("PCA_analysis/r5/Cluster_DP_1.Rdata")

Clusters_all_1<- merge(Cluster_DP1_1, Cluster_DP_1, by = c("CODE_CBNA"))
Clusters_all_1<- merge(Clusters_all_1, Cluster_PY1_1,  by = c("CODE_CBNA"))
Clusters_all_1<- merge(Clusters_all_1, True_clust[, c("CODE_CBNA", "PFG", "K_n")], by = c("CODE_CBNA"))
#save(Clusters_all_1, file =  paste0(folderpath,"Clusters_all_1.Rdata"))


Clusters_all_2<- merge(Cluster_DP1_2, Cluster_DP_2, by = c("CODE_CBNA"))
Clusters_all_2<- merge(Clusters_all_2, Cluster_PY1_2,  by = c("CODE_CBNA"))
Clusters_all_2<- merge(Clusters_all_2, True_clust[, c("CODE_CBNA", "PFG", "K_n")], by = c("CODE_CBNA"))
#save(Clusters_all_2, file =  paste0(folderpath,"Clusters_all_2.Rdata"))


load("PCA_analysis/r5/Cluster_DP1_1.Rdata")
Clusters_all_2$CBN <- as.numeric(Clusters_all_2$CODE_CBNA)
c1= Clusters_all_2[order(Clusters_all_2$CBN),3 ]
c2= Cluster_DP_2$ClustDP
############################################################################################################
#### Plot prior for PY1
x_vec<- 1:112
load("IJulia_part/Cnk_mat_112_H025.Rdata")
pks<- sapply(x_vec, PY_prior,H=112, n=112, alpha=42.48226, sigma=0.5,Cnk_mat=Cnk_112_112_H05)
plot(x_vec, pks)

exp<-sum(x_vec*pks)
exp
