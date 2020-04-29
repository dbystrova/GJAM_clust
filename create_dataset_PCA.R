setwd("~/Documents/GitHub/GJAMF/gjam 4/")
#setwd("~/Bureau/GJAM_paper/GJAMF/gjam 4")
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
library(AUC)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
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

create_DS<- function(rare, data, file_main, file_train, file_test,seed, percent=0.70){
  Sp_cite_counts <- colSums(data[,7:(ncol(data)-2)], na.rm = TRUE, dims = 1)
  Sp_rare<- c(rep(FALSE,6),Sp_cite_counts< rare)
  DS_rare<- data[,!Sp_rare]
  save(DS_rare, file=file_main)
  set.seed(seed)
  smp_size <- floor(percent * nrow(DS_rare))
  train_ind <- sample(seq_len(nrow(DS_rare)), size = smp_size)
  train_rare <- DS_rare[train_ind,]
  save(train_rare, file=file_train)
  test_rare <- DS_rare[-train_ind,]
  save(test_rare, file=file_test)
}
#Test
#lab<- names(data)%in% names(DS_rare)
#names(data)[!lab]

########### Environmental data

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

########## Covariates
PA_env <- merge(B.envm,PA_data_df,by="cite")
#The same dataset [coincides in 1353 cites]
#AB_env<- merge(B.envm,AB_data_df,by="cite")
  
#summary(PA_env)
num_of_NA_in_PA<- rowSums(is.na(PA_env[,7:ncol(PA_env)]))/(ncol(PA_env)-6)
num_of_NA_in_PA_sp<- colSums(is.na(PA_env[,7:ncol(PA_env)]))/nrow(PA_env)
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


###centering PA data

PA_envfin<- PA_envf[!duplicated_env,]
PA_env.pca<- env.pca$scores[,1:2]
PA_envfin$PC1<- PA_env.pca[,1]
PA_envfin$PC2<- PA_env.pca[,2]


############ To delete the species that never detected
Sites_per_species <- colSums(PA_envfin[,7:131], na.rm = TRUE, dims = 1)
sp_list<- names(PA_envfin[,7:131])[Sites_per_species==0]
### To delete the species that never detected
#remove from ydata types never present:14607, 14797
PA_data<- subset(PA_envfin, select = -c(`14607`,`14797`))
#save(PA_data, file = "Bauges_dataset/PA_data_clean_PCA.RData")
PA_data<- load_object("Bauges_dataset/PA_data_clean_PCA.RData")


############ Rare species
set.seed(123)
data<-PA_data
smp_size <- floor(0.70 * nrow(data))
train_ind <- sample(seq_len(nrow(data)), size = smp_size)
train <- data[train_ind,]
test <- data[-train_ind,]
##Create DS  deleting rare species, less than 10, 20 , 30, 50
Sp_cite_counts <- colSums(PA_data[,7:(ncol(PA_data)-2)], na.rm = TRUE, dims = 1)


create_DS(rare=10, data= PA_data,file_main="Bauges_dataset/TrainTestPCA10/PA_pcadata_10.Rdata", file_train="Bauges_dataset/TrainTestPCA10/PA_pcatrain_10.Rdata", file_test="Bauges_dataset/TrainTestPCA10/PA_pcatest_10.Rdata",seed=16, percent=0.70 )
create_DS(rare=20, data= PA_data,file_main="Bauges_dataset/TrainTestPCA20/PA_pcadata_20.Rdata", file_train="Bauges_dataset/TrainTestPCA20/PA_pcatrain_20.Rdata", file_test="Bauges_dataset/TrainTestPCA20/PA_pcatest_20.Rdata",seed=16, percent=0.70 )

###Test
pca10<- load_object("Bauges_dataset/TrainTestPCA10/PA_pcadata_10.Rdata")
S_prev <- sort(colSums(pca10[,7:(ncol(pca10)-2)], na.rm = TRUE, dims = 1), decreasing = TRUE)


pca20<- load_object("Bauges_dataset/TrainTestPCA20/PA_pcadata_20.Rdata")
S_prev <- sort(colSums(pca20[,7:(ncol(pca20)-2)], na.rm = TRUE, dims = 1), decreasing = TRUE)
