## Script for analysis plots

### General include
setwd("~/Documents/GitHub/GJAMF/gjam 4/")
#rm(list=ls())
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
library(corrplot)
library("FactoMineR")
library("ggsci")
li
Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions.R")
source("R/gjam.R")

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

###DATA#############################################
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
iterations=5000
burn_period=1500
K_prior=16
r_reduct = 5



##conditional prediction
columns<-1:ncol(Ydata_train)
ycs<- sample(columns, 10)
y_c_p <-columns[ !columns %in% ycs]


#Preparing species functional groups 
Species_names_groups<- read.csv("Bauges/PFG_Bauges_Description_2017.csv", sep="\t")
# K=16 functional groups
true_names<- as.data.frame(names(table(Species_names_groups$PFG)))
names(true_names)<- c("PFG")
true_names$K_n<- 1:16
Species_names_groups_num<- merge(Species_names_groups,true_names, by="PFG" )



####################################################################################
## Plots for the DIC and choice of r
Tabr5_25<-load_object("PCA_analysis/r5/Fin_tab_r5_25.Rdata")
Tabr5_15<-load_object("PCA_analysis/r5/Fin_tab_r5_L.Rdata")
Tabr5<-load_object("PCA_analysis/r5/Fin_tab_r5.Rdata")


Tabr10<-load_object("PCA_analysis/r10/Fin_tab_r10.Rdata")
Tabr20<-load_object("PCA_analysis/r20/Fin_tab_r20.Rdata")
Full_tab<- rbind(Tabr5[,1:7],Tabr10[,1:7],Tabr20[,1:7])

## Plots for different parameters of "r"

#Tabr20<-load_object("PCA_analysis/r20/fin_tab.Rdata")
df_to_plot<-Full_tab[Full_tab$Parameter=="DIC",] 
df_to_plot2 <- melt(Full_tab, id=c("Parameter", "r"))


p1<- ggplot(df_to_plot2, aes(x=r, y=value, color=variable)) +  geom_point() +
  geom_line()  +  facet_wrap(~ Parameter,scale="free")


#pdf(file = "Final_table_r.pdf")
p1
#dev.off()


## Plots for different number of MCMC samples
Full_tab_MCMC<- rbind(Tabr5[,c(1:6,8)],Tabr5_15[,c(1:6,8)],Tabr5_25[,c(1:6,8)])
df_to_plot_MCMC <- melt(Full_tab_MCMC, id=c("Parameter", "iter"))
#pdf(file = "Final_table_MCMC.pdf")
p2<- ggplot(df_to_plot_MCMC, aes(x=iter, y=value, color=variable)) +  geom_point() +
  geom_line()  +  facet_wrap(~ Parameter,scale="free")
p2
#dev.off()

## So we choose r=5
## Cluster analysis#########################################################

#Species_names_groups_num
Species<- colnames(Ydata)[1: (ncol(Ydata))]
True_clustering<- as.data.frame(Species)
names(True_clustering)<- c("CODE_CBNA")
True_clust<- merge(Species_names_groups_num,True_clustering, by="CODE_CBNA")
#True_clust<- True_clust[order(True_clust$CODE_CBNA),]

##Clusters for different number of iterations
Clust5_L<-load_object("PCA_analysis/r5/Clusters_r5_L.Rdata")
Clust5_25<-load_object("PCA_analysis/r5/Clusters_r5_25.Rdata")
Clust5_25<-load_object("PCA_analysis/r5/Clusters_r5_252.Rdata")
Clust5<-load_object("PCA_analysis/r5/Clusters_r5.Rdata")

##Clusters for different r
Clust10<-load_object("PCA_analysis/r10/Clusters_r10.Rdata")
Clust20<-load_object("PCA_analysis/r20/Clusters_r20.Rdata")


## Two random clusters
Random_cluster <-  sample(1:16, size = 111, replace = TRUE)
Random_cluster_w <-  sample(1:16, size = 111, replace = TRUE, prob=p_w)

Clusters_by_rDP<-data.frame(PFG=True_clust$K_n,R5=Clust5$DP,R10=Clust10$DP,R20=Clust20$DP, Random=Random_cluster )
Clusters_by_rDP1<-data.frame(PFG=True_clust$K_n,R5=Clust5$DP1,R10=Clust10$DP1,R20=Clust20$DP1, Random=Random_cluster )
Clusters_by_rDP2<-data.frame(PFG=True_clust$K_n,R5=Clust5$DP2,R10=Clust10$DP2,R20=Clust20$DP2, Random=Random_cluster )
Clusters_by_rPY1<-data.frame(PFG=True_clust$K_n,R5=Clust5$PY1,R10=Clust10$PY1,R20=Clust20$PY1, Random=Random_cluster )
Clusters_by_rPY2<-data.frame(PFG=True_clust$K_n,R5=Clust5$PY2,R10=Clust10$PY2,R20=Clust20$PY2, Random=Random_cluster )

Mat_dist<- matrix(NA, nrow=5, ncol=5)
for(j in 1:5){
  for (k in 1:5){
    Mat_dist[j,k]= vi.dist(Clusters_by_r[,j],Clusters_by_r[,k])
  }
}
kable(as.data.frame(round(Mat_dist,3)), col.names = colnames(Clusters_by_r))

M_all<- data.frame()
for(l in 1:4){
  Clusters_by_r<-data.frame(PFG=True_clust$K_n,R5=Clust5[,l+1],R10=Clust10[,l+1],R20=Clust20[,l+1])
  Mat_dist<- matrix(NA, nrow=4, ncol=4)
  for(j in 1:4){
    for (k in 1:4){
      Mat_dist[j,k]= arandi(Clusters_by_r[,j],Clusters_by_r[,k])
    }
  }
  M_print<- as.data.frame(round(Mat_dist,3))
  rownames(M_print)<- c("PFG","R5","R10","R20")
  M_all<- rbind(M_all,cbind(M_print,rep(colnames(Clust5)[l+1],4), c("PFG","R5","R10","R20")))
}

colnames(M_all)<- c("PFG","R5","R10","R20","Model","R")
M_f<- melt(M_all, id=c("Model","R"))
Dist.heatmap <- ggplot(data = M_f, mapping = aes(x = variable,
                                                 y = R,
                                                 fill = value, group=)) +  facet_wrap(.~ Model,scale="free")+
  geom_tile() +
  xlab(label = "AR distance")
Dist.heatmap
#pdf(file = "Final_cluster_distances_AR.pdf")
Dist.heatmap
grid.table(M_all[,1:5])
#dev.off()



M_all<- data.frame()
for(l in 1:4){
  Clusters_by_it<-data.frame(PFG=True_clust$K_n,It10=Clust5[,l+1],It15=Clust5_L[,l+1],It25=Clust5_25[,l+1])
  Mat_dist<- matrix(NA, nrow=4, ncol=4)
  for(j in 1:4){
    for (k in 1:4){
      Mat_dist[j,k]= arandi(Clusters_by_it[,j],Clusters_by_it[,k])
    }
  }
  M_print<- as.data.frame(round(Mat_dist,3))
  M_all<- rbind(M_all,cbind(M_print,rep(colnames(Clust5)[l+1],4), c("PFG","It10","It15","It25")))
}

colnames(M_all)<- c("PFG","It10","It15","It25")
M_f<- melt(M_all, id=c("Model","R"))
Dist.heatmap <- ggplot(data = M_f, mapping = aes(x = variable,
                                                 y = R,
                                                 fill = value, group=)) +  facet_wrap(.~ Model,scale="free")+
  geom_tile() +
  xlab(label = "AR distance")
Dist.heatmap
#pdf(file = "Final_cluster_distances_AR_it.pdf")
Dist.heatmap
grid.table(M_all[,1:5])
#dev.off()



########### Cluster distances with random 


M_all<- data.frame()
Clusters_by_Rand<-data.frame(PFG=True_clust$K_n,DP=Clust5_25$DP,DP1=Clust5_25$DP1,
                             DP2=Clust5_25$DP2, PY1=Clust5_25$PY1, PY2=Clust5_25$PY2, RU=Random_cluster, RW= Random_cluster_w)
Mat_dist<- matrix(NA, nrow=8, ncol=8)
  for(j in 1:8){
    for (k in 1:8){
      Mat_dist[j,k]= vi.dist(Clusters_by_Rand[,j],Clusters_by_Rand[,k])
    }
}
M_print<- as.data.frame(round(Mat_dist,3))
colnames(M_print)<- colnames(Clusters_by_Rand)
M_print$Model <- colnames(Clusters_by_Rand)
M_rand<- melt(M_print, id=c("Model"))
M_rand$Model<- factor(M_rand$Model,levels = colnames(M_print)[1:8])
Dist_rand.heatmap <- ggplot(data = M_rand, mapping = aes(x = variable,
                                                 y = Model,
                                                 fill = value, color='white')) + 
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 2)), size=2)+
  scale_fill_gradient(low = 'white', high = 'black', space = 'Lab', na.value = 'white')+
  xlab(label = "VI distance for final cluster with Random clusters ")
Dist_rand.heatmap

pdf(file = "Final_cluster_distances_VI_rand.pdf")
print(Dist_rand.heatmap)
#grid.table(M_all[,1:5])
dev.off()


#### Anova





do.call(rbind, lapply(summary.aov(aov(formulaov, data=dataCluster)), as.data.frame))




##############################################################################################################

##Clusters explained by the traits

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

expandSigma_rmd <- function(sigma, S){
  
  ss <- diag(S)
  ss[lower.tri(ss,diag=T)] <- sigma
  ss[upper.tri(ss)] <- t(ss)[upper.tri(ss)]
  ss
}




convert_to_m<-function(ar, S){
  C <- matrix(0,d,d)
  i.lwr <- which(lower.tri(C, diag = TRUE), arr.ind=TRUE)
  C[i.lwr] <- ar
  C<-makeSymm(C)
  return(t(C))
}

postH<- apply(fit_gjam$chains$sgibbs,2, quantile,0.95)
postL<-apply(fit_gjam$chains$sgibbs,2, quantile,0.05)
post_mean<-apply(fit_gjam$chains$sgibbs,2, mean)


pH<-convert_to_m(postH, 112)
pL<-convert_to_m(postL,113)
post_mean_s<-convert_to_m(post_mean)
S_mean<-cov2cor(post_mean_s)
#R_sign<-cov2cor(mod_gjam1$parameters$sigMu)*(!(pH>0 & pL<0))
R_sign<-cov2cor(post_mean_s)*(!(pH>0 & pL<0))

sgibbs<-abind(gj_mod$m1$chains$sgibbs[-(1:burn),],gj_mod$m2$chains$sgibbs[-(1:burn),],along=1)

tau<-array(NA,dim=c(data$J,data$J,dim(sgibbs)[1]))
for(j in 1:dim(sgibbs)[1]){
  ss <- expandSigma_rmd(sgibbs[j,], S = data$J)
  si <- solve(ss)
  tau[,,j] <- -cov2cor(si)
}

tau_mean<-apply(tau,c(1,2), mean)
tau_HI<-apply(tau,c(1,2),quantile,0.95)
tau_LO<-apply(tau,c(1,2),quantile,0.05)

Tau_sign<-tau_mean*(!(tau_HI>0 & tau_LO<0))

x<- data$X
#p1<-gjamPredict(mod_gjam1)
mu<-array(NA,dim=c(data$n,data$J,it))
for(k in 1:it){
  for(j in 1:data$J){
    
    mu[,j,k] <- pnorm(x%*%gj_mod$m1$chains$bgibbs[k,(3*(j-1)+1):(3*j)])
    
  }
}

par(mfrow=c(2,3),oma = c(1, 1, 1, 1))
corrplot(cor(data$Y), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Correlation cor(Y)")
corrplot(S_mean, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("R")
corrplot(gj_mod$m1$parameters$ematrix, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("E matrix")
corrplot(Tau_sign, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Tau signif")
corrplot(R_sign, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("R signif")
corrplot(interact, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions")

sgibbs<-fit_gjam$chains$sgibbs[burn_period:iterations,]
sigErrGibbs<-fit_gjam$chains$sigErrGibbs[burn_period:iterations]
kgibbs<-fit_gjam$chains$kgibbs[burn_period:iterations,]
sigma<-invsigma<-array(NA,dim=c(S,S,iterations-burn_period))


#sgibbs<-mod_gjam_red$chains$sgibbs
#sigErrGibbs<-mod_gjam_red$chains$sigErrGibbs
#kgibbs<-mod_gjam_red$chains$kgibbs
N<-fit_gjam$modelList$reductList$N
r<-fit_gjam$modelList$reductList$r
N_dim<-iterations-burn_period
for(j in 1:N_dim){
  Z  <- matrix(sgibbs[j,],N,r)
  sigma[,,j] <- .expandSigma(sigErrGibbs[j], S, Z = Z, kgibbs[j,], REDUCT = T) #sigma
  invsigma[,,j] <- invWbyRcpp(sigErrGibbs[j], Z[kgibbs[j,],]) #inverse sigma
} 

sigma_mean<-apply(sigma,c(1,2),mean) 
sigma_q05<-apply(sigma,c(1,2),quantile,0.05) 
sigma_q95<-apply(sigma,c(1,2),quantile,0.95) 
Sigma_sign<--cov2cor(sigma_mean*(!(sigma_q95>0 & sigma_q05<0)))

invsigma_mean<-apply(invsigma,c(1,2),mean) 
invsigma_q05<-apply(invsigma,c(1,2),quantile,0.05) 
invsigma_q95<-apply(invsigma,c(1,2),quantile,0.95) 
INVSigma_sign<--cov2cor(sigma_mean*(!(invsigma_q95>0 & invsigma_q05<0)))

cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))

gcols = colorRampPalette(c( "White", "White", "Black"))


corrplot(Sigma_sign, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Correlation from dim. reduction")
corrplot(INVSigma_sign, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")


# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix( INVSigma_sign, weighted=T, mode="undirected", diag=F)

# Basic chart
plot(network)



Cluster_DP<- data.frame( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustDP=Clust5_25$DP )
Compare_clust<- merge(Cluster_DP, True_clust,by = "CODE_CBNA")


Cluster_DP1<- data.frame( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustDP1=Clust5_25$DP1 )
Compare_clust1<- merge(Cluster_DP1, True_clust,by = "CODE_CBNA")
kable(table(Compare_clust1$ClustDP1,Compare_clust1$PFG), caption="DP1")



Cluster_DP2<- data.frame( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustDP2=Clust5_25$DP2 )
Compare_clust2<- merge(Cluster_DP2, True_clust,by = "CODE_CBNA")
kable(table(Compare_clust2$ClustDP2,Compare_clust2$PFG), caption = "DP2")


Cluster_PY1<- data.frame(CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustPY1=Clust5_25$PY1 )
Compare_clust3<- merge(Cluster_PY1, True_clust,by = "CODE_CBNA")
kable(table(Compare_clust3$ClustPY1,Compare_clust3$PFG), caption = "PY1")


Cluster_PY2<- data.frame( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustPY2=Clust5_25$PY2 )
Compare_clust4<- merge(Cluster_PY2, True_clust,by = "CODE_CBNA")
kable(table(Compare_clust4$ClustPY2,Compare_clust4$PFG), caption = "PY2")

All_clusters<-data.frame(PFG=True_clust$K_n,DP=Clust5_25$DP,DP1=Clust5_25$DP1,DP2=Clust5_25$DP2,PY1=Clust5_25$PY1,PY2=Clust5_25$PY2)



########## Traits
### Traits
T_PFGdata<- load_object("Bauges/PFG.mat.traits.pfg.RData")
T_data<- load_object("Bauges/PFG.mat.traits.sp.RData")
New_traits_data<- load_object("Bauges/traitsBauges.RData")

### Visualize in LF

T_PFGdata<- load_object("Bauges/PFG.mat.traits.pfg.RData")
T_data<- load_object("Bauges/PFG.mat.traits.sp.RData")

Compare_clust1$species<- sapply(Compare_clust1$CODE_CBNA, function(x) paste0("X",x))
Compare_clust$species<- sapply(Compare_clust$CODE_CBNA, function(x) paste0("X",x))
Compare_clust2$species<- sapply(Compare_clust2$CODE_CBNA, function(x) paste0("X",x))
Compare_clust3$species<- sapply(Compare_clust3$CODE_CBNA, function(x) paste0("X",x))
Compare_clust4$species<- sapply(Compare_clust4$CODE_CBNA, function(x) paste0("X",x))

traits<-merge(T_data,Compare_clust1[,c("ClustDP1","species")], by="species" )
traits<-merge(traits,Compare_clust[,c("ClustDP","species")], by="species" )
traits<-merge(traits,Compare_clust2[,c("ClustDP2","species")], by="species" )
traits<-merge(traits,Compare_clust3[,c("ClustPY1","species")], by="species" )
traits<-merge(traits,Compare_clust4[,c("ClustPY2","species")], by="species" )

row.names(traits)<- traits$species
#traits_nm<-  na.omit(traits[,c("dispersal","light","soil_contrib","soil_tolerance","height")])
#traits_nm<-  na.omit(traits)
L<- rowSums(is.na(traits[,c("dispersal","light","soil_contrib","soil_tolerance","height")])) > 0
traits<- traits[!L,]

tdf<- traits[,-c(5:7)]
col_names <- c("PFG" ,"type","dispersal", "light","soil_tolerance","soil_contrib", "ClustDP1", "ClustDP" ,"ClustDP2" ,  "ClustPY1"      , "ClustPY2")

# do it for some names in a vector named 'col_names'
tdf[col_names] <- lapply(tdf[col_names] , factor)
res.famd <- FAMD(tdf,sup.var=c(1,2,9:13), graph = FALSE)


fviz_famd_ind(res.famd, 
              habillage =  "PFG", # color by groups 
              addEllipses = TRUE, ellipse.type = "confidence", 
              repel = TRUE, # Avoid text overlapping,
              mean.point = TRUE,
              label="none",
              legend.title = "Groups") +
  theme_minimal() +
  theme(legend.position = "bottom")


fviz_famd_ind(res.famd, 
              habillage =  "ClustDP2", # color by groups 
              # addEllipses = TRUE, ellipse.type = "confidence", 
              repel = TRUE, # Avoid text overlapping,
              mean.point = TRUE,
              label= "none") +
  theme_minimal() +
  theme(legend.position = "bottom")



fviz_famd_ind(res.famd, 
              habillage =  "ClustPY1", # color by groups 
              # addEllipses = TRUE, ellipse.type = "confidence", 
              repel = TRUE, # Avoid text overlapping,
              mean.point = TRUE,
              label= "none") +
  theme_minimal() +
  theme(legend.position = "bottom")


New_traits_ds<- load_object("Bauges/traitsBauges.RData")
SP<- Species_names_groups[,c(1,3)]
colnames(SP)<- c("code_cbna","species_name")
New_traits_data<- merge(New_traits_ds[,c(1:4,7,13)],SP,by="code_cbna")
Clusters<- Reduce(function(x, y) merge(x, y, by="species",all=T), list(Compare_clust1[,c("ClustDP1","species")], Compare_clust[,c("ClustDP","species")], Compare_clust2[,c("ClustDP2","species")],Compare_clust3[,c("ClustPY1","species")],Compare_clust4[,c("ClustPY2","species")] ))
New_traits_data<- merge(New_traits_data,as.data.frame(Clusters),by="species")

row.names(New_traits_data)<- New_traits_data$species_name
L<- rowSums(is.na(New_traits_data[,c("LNC","LCCp","LDMC","SLA")])) > 0
New_traits_data<- New_traits_data[!L,]

td<- New_traits_data[,c(3:6,8)]

res.pca <- prcomp( New_traits_data[,c(3:6)],  scale = TRUE)
fviz_pca_ind(res.pca, col.ind=new_traits_nm$PFG, repel=TRUE)


p <- fviz_pca_ind(res.pca, label="none", geom="point",habillage=new_traits_nm$PFG)
print(p)


p <- fviz_pca_ind(res.pca, label="none", habillage=New_traits_data$ClustDP)
print(p)


p <- fviz_pca_biplot(res.pca, label="var", habillage=New_traits_data$ClustDP1)
print(p)



p <- fviz_pca_ind(res.pca, label="none", habillage=new_traits_nm$ClustDP, addEllipses = TRUE) 
print(p)

library(viridis)

tdf[col_names] <- lapply(tdf[col_names] , factor)
res.famd <- FAMD(tdf,sup.var=c(1,2,9:13), graph = FALSE)
fviz_famd_var(res.famd, repel=TRUE)+ scale_color_brewer(palette="Dark2") +
  fviz_famd_ind(res.famd, col.ind= "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=TRUE)


fviz_famd_ind(res.famd, 
              habillage =  "PFG", # color by groups 
              addEllipses = TRUE, ellipse.type = "confidence", 
              repel = TRUE, # Avoid text overlapping,
              mean.point = TRUE,
              label="none",
              legend.title = "Groups") +scale_color_viridis(discrete = TRUE, option = "D",alpha=rep(0.3,nrow(tdf)))+
  theme_minimal() +
  theme(legend.position = "bottom")


fviz_famd_ind(res.famd, 
              habillage =  "ClustDP", # color by groups 
              # addEllipses = TRUE, ellipse.type = "confidence", 
              repel = TRUE, # Avoid text overlapping,
              mean.point = TRUE,
              label= "none") +scale_color_viridis(discrete = TRUE, option = "D")+
  theme_minimal() +
  theme(legend.position = "bottom")


library(afex)

fviz_famd_ind(res.famd, 
              habillage =  "ClustDP2", # color by groups 
              # addEllipses = TRUE, ellipse.type = "confidence", 
              repel = TRUE, # Avoid text overlapping,
              mean.point = TRUE,
              label= "none") +scale_color_viridis(discrete = TRUE, option = "D")+
  theme_minimal() +
  theme(legend.position = "bottom")



############################## Function of traits  #############################################
library(nnet)
new_df<- data.frame(FAMDC1=res.famd$ind$coord[,1], FAMDC2=res.famd$ind$coord[,2], ClustDP=traits$ClustDP)
new_df1<- data.frame(FAMDC1=res.famd$ind$coord[,1], FAMDC2=res.famd$ind$coord[,2], ClustDP1=traits$ClustDP1)
new_df2<- data.frame(FAMDC1=res.famd$ind$coord[,1], FAMDC2=res.famd$ind$coord[,2], ClustDP2=traits$ClustDP2)
new_df3<- data.frame(FAMDC1=res.famd$ind$coord[,1], FAMDC2=res.famd$ind$coord[,2], ClustPY1=traits$ClustPY1)
new_df4<- data.frame(FAMDC1=res.famd$ind$coord[,1], FAMDC2=res.famd$ind$coord[,2], ClustPY2=traits$ClustPY2)
new_df5<- data.frame(FAMDC1=res.famd$ind$coord[,1], FAMDC2=res.famd$ind$coord[,2], ClustPFG=traits$PFG)


#model <- multinom(ClustDP ~ FAMDC1 +FAMDC2, data=new_df)
# Run the model
model0 <- multinom(PFG ~ type+ height+palatability+dispersal+light+soil_contrib+soil_tolerance, data=traits)
summary(model0)
#model <- multinom(ClustDP ~ FAMDC1 +FAMDC2, data=new_df)
#model2 <- multinom(ClustDP1 ~ FAMDC1 +FAMDC2, data=new_df2)

model <- multinom(ClustDP ~ type+ height+palatability+dispersal+light+soil_contrib+soil_tolerance, data=traits)
summary(model)
model2 <- multinom(ClustDP1 ~ type+ height+palatability+dispersal+light+soil_contrib+soil_tolerance, data=traits)
summary(model2)
model3 <- multinom(ClustDP2 ~ type+ height+palatability+dispersal+light+soil_contrib+soil_tolerance, data=traits)
summary(model3)
model4 <- multinom(ClustPY1 ~ type+ height+palatability+dispersal+light+soil_contrib+soil_tolerance, data=traits)
summary(model4)
model5 <- multinom(ClustPY2 ~ type+ height+palatability+dispersal+light+soil_contrib+soil_tolerance, data=traits)
summary(model5)
model6 <- multinom(PFG ~ type+ height+palatability+dispersal+light+soil_contrib+soil_tolerance, data=traits)
summary(model6)


kable(data.frame(Metric="AIC",DP=model$AIC, DP1=model2$AIC,DP2=model3$AIC,PY1=model4$AIC,PY2=model5$AIC, PFG= model6$AIC))
kable(data.frame(Metric="Residual_error",DP=model$deviance, DP1=model2$deviance,DP2=model3$deviance,PY1=model4$deviance,PY2=model5$deviance, PFG= model6$deviance))



##########################################################################################



# MANOVA test
res.man <- manova(cbind(height, longevity,maturity, palatability,dispersal,light,soil_contrib,soil_tolerance) ~ PFG, data = traits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(height,longevity) ~ ClustDP, data = traits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(height,soil_tolerance) ~ ClustDP2, data = traits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(height,soil_tolerance) ~ ClustPY1, data = traits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(height, longevity,maturity, palatability,dispersal,light,soil_contrib,soil_tolerance) ~ ClustDP1, data = traits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(height, longevity,maturity, palatability,dispersal,light,soil_contrib,soil_tolerance) ~ ClustDP2, data = traits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(height, longevity,maturity, palatability,dispersal,light,soil_contrib,soil_tolerance) ~ ClustPY1, data = traits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(height, longevity,maturity, palatability,dispersal,light,soil_contrib,soil_tolerance) ~ ClustPY2, data = traits)
summary(res.man)
summary.aov(res.man)


#### Additional traits








###########################Conditional prediction########################################################


### Function of another traits


New_traits_data<- load_object("Bauges/traitsBauges.RData")

Ntraits<- New_traits_data[,c(2:7,11,13,14)]
#traits<-merge(T_data,Compare_clust1[,c("ClustDP1","species")], by="species" )
#traits<-T_data
Ntraits<-merge(Ntraits,Compare_clust[,c("ClustDP","species")], by="species" )
Ntraits<-merge(Ntraits,Compare_clust1[,c("ClustDP1","species")], by="species" )

Ntraits<-merge(Ntraits,Compare_clust2[,c("ClustDP2","species")], by="species" )
Ntraits<-merge(Ntraits,Compare_clust3[,c("ClustPY1","species")], by="species" )
Ntraits<-merge(Ntraits,Compare_clust4[,c("ClustPY2","species")], by="species" )

new_traits_nm<-  na.omit(Ntraits[, c(2:4,7,9,10,11)])


res.pca <- prcomp(new_traits_nm[,-c(5:7)],  scale = TRUE)
fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, col.ind=new_traits_nm$PFG, repel=TRUE)


p <- fviz_pca_ind(res.pca, label="none", habillage=new_traits_nm$PFG)
print(p)


p <- fviz_pca_ind(res.pca, label="none", habillage=new_traits_nm$ClustDP)
print(p)


p <- fviz_pca_ind(res.pca, label="none", habillage=new_traits_nm$ClustPY1) 
print(p)


tdf[col_names] <- lapply(tdf[col_names] , factor)
res.famd <- FAMD(tdf,sup.var=c(1,2,9:13), graph = FALSE)
fviz_famd_var(res.famd, repel=TRUE)+ scale_color_brewer(palette="Dark2") +
  fviz_famd_ind(res.famd, col.ind= "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel=TRUE)



#################################### Run the model ##############################
model0 <- multinom(PFG ~  LNC+ LCCp+LDMC + LKC+LPC+ SLA + ELLNIT_INDK , data=Ntraits)
summary(model0) 
#model <- multinom(ClustDP ~ FAMDC1 +FAMDC2, data=new_df)
#model2 <- multinom(ClustDP1 ~ FAMDC1 +FAMDC2, data=new_df2)

model <- multinom(ClustDP ~  LNC+ LCCp+LDMC + LKC+LPC+ SLA + ELLNIT_INDK , data=Ntraits)
summary(model)
model2 <- multinom(ClustDP1 ~ LNC+ LCCp+LDMC + LKC+LPC+ SLA + ELLNIT_INDK , data=Ntraits)
summary(model2)
model3 <- multinom(ClustDP2 ~  LNC+ LCCp+LDMC + LKC+LPC+ SLA + ELLNIT_INDK , data=Ntraits)
summary(model3)
model4 <- multinom(ClustPY1 ~  LNC+ LCCp+LDMC + LKC+LPC+ SLA + ELLNIT_INDK , data=Ntraits)
summary(model4)
model5 <- multinom(ClustPY2 ~  LNC+ LCCp+LDMC + LKC+LPC+ SLA + ELLNIT_INDK , data=Ntraits)
summary(model5)
model6 <- multinom(PFG ~  LNC+ LCCp+LDMC + LKC+LPC+ SLA + ELLNIT_INDK , data=Ntraits)
summary(model6)
kable(data.frame(Metric="AIC",DP=model$AIC, DP1=model2$AIC,DP2=model3$AIC,PY1=model4$AIC,PY2=model5$AIC, PFG= model6$AIC))
kable(data.frame(Metric="Residual_error",DP=model$deviance, DP1=model2$deviance,DP2=model3$deviance,PY1=model4$deviance,PY2=model5$deviance, PFG= model6$deviance))
#################################

# MANOVA test
res.man <- manova(cbind(LNC,LCCp,LDMC,LPC,SLA ,ELLNIT_INDK) ~ PFG, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC,LPC,SLA ,ELLNIT_INDK)  ~ ClustDP, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC,LPC,SLA ,ELLNIT_INDK)  ~ ClustDP2, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC,LPC,SLA ,ELLNIT_INDK)  ~ ClustDP1, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC,LPC,SLA ,ELLNIT_INDK)  ~ ClustPY1, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC,LPC,SLA ,ELLNIT_INDK)  ~ ClustPY2, data = Ntraits)
summary(res.man)
summary.aov(res.man)


# MANOVA test
res.man <- manova(cbind(LNC,LCCp,LDMC ,SLA ,ELLNIT_INDK) ~ PFG, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC ,SLA ,ELLNIT_INDK)  ~ ClustDP, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC, SLA ,ELLNIT_INDK)  ~ ClustDP2, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC ,SLA ,ELLNIT_INDK)  ~ ClustDP1, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC ,SLA ,ELLNIT_INDK)  ~ ClustPY1, data = Ntraits)
summary(res.man)
summary.aov(res.man)

res.man <- manova(cbind(LNC,LCCp,LDMC ,SLA ,ELLNIT_INDK)  ~ ClustPY2, data = Ntraits)
summary(res.man)
summary.aov(res.man)




T_M<- melt(traits, id=c("species","type","height","longevity", "maturity", "palatability", "dispersal", "light", "soil_contrib", "soil_tolerance"))
T_M$value <- as.factor(T_M$value)


pdf(file = "Box_plot_height.pdf")

p<- ggplot(T_M, aes(x=value, y=height)) + 
  geom_boxplot() + facet_wrap(~ variable,scale="free")
print(p)

dev.off()
pdf(file = "Box_plot_longevity.pdf")

p<- ggplot(T_M, aes(x=value, y=longevity)) + 
  geom_boxplot() + facet_wrap(~ variable,scale="free")
p

print(p)

dev.off()







res.man <- manova(cbind(LNC,LCCp,LDMC ,SLA ,ELLNIT_INDK)  ~ ClustPY2, data = Ntraits)
summary(res.man)
summary.aov(res.man)




########## Final table to PDF

library(gridExtra)


pdf(file = "Final_table_prediction.pdf")


grid.table(Tabr5_25[1:5,1:6])
dev.off()


########## Final table to PDF clustering


pdf(file = "Final_table_prediction_K_dist.pdf")


grid.table(Tabr5_25[6:9,1:6])
dev.off()



pdf(file = "Final_table_prediction_K.pdf")


grid.table(Tabr5_25[10:12,1:6])
dev.off()

## Species condi
## Species conditioning on 

SP <- colnames(Ydata)[ycs]
names_conditon_species <- True_clust[True_clust$CODE_CBNA%in%as.numeric(SP),3]


#######


T_M<- melt(traits, id=c("species","type","height","longevity", "maturity", "palatability", "dispersal", "light", "soil_contrib", "soil_tolerance"))




do.call(rbind, lapply(summary.aov(aov(formulaov, data=dataCluster)), as.data.frame))



