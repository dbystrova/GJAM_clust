#PCA Plot


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


Rcpp::sourceCpp('src/cppFns.cpp')
source("R/gjamHfunctions.R")
source("R/gjam.R")

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

###Load coordinates 
B_coords_xy<- load_object("Bauges/DB.XY.RData")
#Load abundance/PA data from folder
PA_data<-load_object("Bauges/DOM.mat.sites.species.PA.RData")
AB_data<-load_object("Bauges/DOM.mat.sites.species.abund.RData")


#Presence/Absence
PA_data_df<- as.data.frame(PA_data)
PA_data_df$cite<- rownames(PA_data)
#Abundance
AB_data_df<- as.data.frame(AB_data)
AB_data_df$cite<- rownames(AB_data)
#Intersection


###### Environmental covarites
folder.name="Bauges"
zone.name="ENV_VARIABLES"
zone.env.folder="EOBS_1970_2005"
zone.env.variables=c("bio_1_0","bio_12_0","bio_19_0","bio_8_0","slope")

##function from the package FATEHD
getSDM_env = function(zone.name, zone.env.folder, zone.env.variables, maskSimul)
{
  env.files = list.files(path = paste0(folder.name, "/",zone.name, "/", zone.env.folder)
                         , pattern = paste0(paste0(zone.env.variables, ".img", collapse = "|")
                                            , "|"
                                            , paste0(zone.env.variables, ".tif", collapse = "|"))
                         , full.names = TRUE)
  zone.env.stk.CALIB = raster::stack(env.files)
  maskSimul=raster("Bauges/MASK_100m.tif")
  origin(maskSimul) = origin(zone.env.stk.CALIB)
  zone.env.stk.PROJ = stack(zone.env.stk.CALIB * maskSimul)
  names(zone.env.stk.PROJ) = names(zone.env.stk.CALIB)
  
  return(list(env.CALIB = zone.env.stk.CALIB
              , env.PROJ = zone.env.stk.PROJ))
}

new_B_env<-getSDM_env(zone.name, zone.env.folder, zone.env.variables, maskSimul=raster("Bauges/MASK_100m.tif")) 
B_env_raster<- new_B_env$env.PROJ
B_env<-as.data.frame(extract(B_env_raster, B_coords_xy))
B_env$cite<- rownames(B_coords_xy)
##Delete 0 values and NA's for the environmental variables
NAs_values<- is.na(B_env$bio_1_0)&is.na(B_env$bio_12_0)&is.na(B_env$bio_19_0)&is.na(B_env$bio_8_0)&is.na(B_env$slope)
B_env_1<- B_env[!NAs_values,]
zeros_values<- (B_env_1$bio_1_0==0)&(B_env_1$bio_12_0==0)&(B_env_1$bio_19_0==0)&(B_env_1$bio_8_0==0)&(B_env_1$slope==0)
B.envm<- B_env_1[!zeros_values,]

PA_env <- merge(B.envm,PA_data_df,by="cite")
#The same dataset [coincides in 1353 cites]
#AB_env<- merge(B.envm,AB_data_df,by="cite")

#summary(PA_env)
num_of_NA_in_PA<- rowSums(is.na(PA_env[,7:131]))/(ncol(PA_env)-6)
num_of_NA_in_PA_sp<- colSums(is.na(PA_env[,7:131]))/nrow(PA_env)
p1<- ggplot(as.data.frame(cbind(1:nrow(PA_env),num_of_NA_in_PA)), aes(x=V1, y=num_of_NA_in_PA)) + 
  geom_point(col="steelblue", size=2) +
  geom_hline(yintercept=0.6,linetype="dashed", color = "red")+
  labs(title="Proportion of missing entries in rows in PA dataset", subtitle=" number of species by cite", y="Species richness", x="Cite", caption="")
#p1

p2<- ggplot(as.data.frame(cbind(1:(ncol(PA_env)-6),num_of_NA_in_PA_sp)), aes(x=V1, y=num_of_NA_in_PA_sp)) + 
  geom_point(col="steelblue", size=2) +
  labs(title="Number of missing entries in columns in PA dataset", subtitle=" number of species by cite", y="Species richness", x="Cite", caption="")
#p2

grid.arrange(p1, p2, ncol=2)
#####################################################################################################################
#we consider dataset without missings
PA_envf<- PA_env[which(num_of_NA_in_PA<0.6),]
env_data_P<- PA_envf[,1:6]


#####################################################################################################################
## Check for duplicates
#No dup on cites
#Duplictes on env covariates
duplicated_env<- env_data_P[,2:6] %>% duplicated()
env_data_dup<- env_data_P[duplicated_env,2:6] #226 dup covariates
env_data_nodup<- env_data_P[!duplicated_env,2:6] #>200 dup covariates



env.pca <- princomp(as.data.frame(apply(env_data_nodup,MARGIN=2, scale)))
fviz_eig(env.pca)
pca_plot<- fviz_pca_var(env.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


pdf(file = "Plots/PCA_plot.pdf", width=6, height =4)
plot(pca_plot)
dev.off()

#####################################################################################################################
###centering PA data

PA_envfin<- PA_envf[!duplicated_env,]
PA_env.pca<- env.pca$scores[,1:2]
PA_envfin$PC1<- PA_env.pca[,1]
PA_envfin$PC2<- PA_env.pca[,2]

