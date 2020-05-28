### This is the version of analysis where use two chains
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
fit_gjam<- load_object(paste0(folderpath,"fit_gjam.Rdata"))
fit_gjamDP1<- load_object(paste0(folderpath,"fit_gjamDP1.Rdata"))
fit_gjamPY1<- load_object(paste0(folderpath3,"fit_gjamPY1.Rdata"))
fit_gjamPY2<- load_object(paste0(folderpath,"fit_gjamPY2.Rdata"))
fit_gjamPY12<- load_object(paste0(folderpath,"fit_gjamPY1_05.Rdata"))


## Load models 2nd run
fit_gjam_2<- load_object(paste0(folderpath2,"fit_gjam.Rdata"))
fit_gjamDP1_2<- load_object(paste0(folderpath2,"fit_gjamDP1.Rdata"))
fit_gjamPY1_2<- load_object(paste0(folderpath4,"fit_gjamPY1.Rdata"))
fit_gjamPY12_2<- load_object(paste0(folderpath2,"fit_gjamPY1_05.Rdata"))

#fit_gjamPY2<- load_object(paste0(folderpath,"fit_gjamPY2.Rdata"))

########################################Prediction########################################################
#Prediction out-of sample on xtest
new_out <- list(xdata =xdata_test,  nsim = 1000) # effort unchanged 
#Prediction in-sample
new_in <- list(xdata =xdata_train,  nsim = 1000) # effort unchanged 
#Conditional prediction 
new_cond <- list(ydataCond = Ydata_test[,ycs], xdata =xdata_test, nsim=1000)   # cond on obs CA data
# WAUC 


model_prediction_summary<- function(model, list_out_s,list_in_s , list_cond, Ytest=Ydata_test, Ytrain=Ydata_train, pw=p_w){
  pred_out_s <- gjamPredict(output = model, newdata = list_out_s)
  #save(predict_gjam, file = paste0(folderpath,"predict_gjam.Rdata"))
  AUC_ous<-vector()
  for(i in 1:ncol(Ytest)){ 
    label<- is.na(Ytest[,i])
    predict<- pred_out_s$sdList$yMu[!label,i]
    test_value<- Ytest[!label,i]
    if(sum(test_value)>0){
      AUC_ous<-c(AUC_ous,auc(AUC::roc(predict,factor(test_value))))
    }
  }
  pred_in_s <- gjamPredict(output = model, newdata = list_in_s)
  AUC_GJAMDin<-vector()
  for(i in 1:ncol(Ytrain)){ 
    label<- is.na(Ytrain[,i])
    predict<- pred_in_s$sdList$yMu[!label,i]
    test_value<- Ytrain[!label,i]
    if(sum(test_value)>0){
      AUC_GJAMDin<-c(AUC_GJAMDin,auc(roc(predict,factor(test_value))))
    }
  }
  
  predict_gjam_cond  <- gjamPredict(output = model, newdata = list_cond)
  AUC_GJAM_cond<-vector()
  for(i in y_c_p){ 
    label<- is.na(Ytest[,i])
    predict<- predict_gjam_cond$sdList$yMu[!label,i]
    test_value<- Ytest[!label,i]
    if(sum(test_value)>0){
      AUC_GJAM_cond<-c(AUC_GJAM_cond,auc(roc(predict,factor(test_value))))
    }
  }
  return(list(AUC_out=AUC_ous,AUC_in=AUC_GJAMDin, AUC_cond= AUC_GJAM_cond, WAUC= AUC_ous*pw))
}

####
gjamDP<- model_prediction_summary(model=fit_gjam, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamDP1<- model_prediction_summary(model=fit_gjamDP1, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamPY1<- model_prediction_summary(model=fit_gjamPY1, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamPY2<- model_prediction_summary(model=fit_gjamPY2, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)

gjamDP_2<- model_prediction_summary(model=fit_gjam_2, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamDP1_2<- model_prediction_summary(model=fit_gjamDP1_2, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamPY1_2<- model_prediction_summary(model=fit_gjamPY1_2, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)

##Correlation list
Cor_DP<- list()
Cor_DP1<- list()
Cor_PY1<- list()

#Convergence of parameters
Cor_DP<- list.append(Cor_DP, cor(gjamDP$AUC_out, gjamDP_2$AUC_out), cor(gjamDP$AUC_in, gjamDP_2$AUC_in))
Cor_DP1<- list.append(Cor_DP1, cor(gjamDP1$AUC_out, gjamDP1_2$AUC_out), cor(gjamDP1$AUC_in, gjamDP1_2$AUC_in))
Cor_PY1<- list.append(Cor_PY1, cor(gjamPY1$AUC_out, gjamPY1_2$AUC_out), cor(gjamPY1$AUC_in, gjamPY1_2$AUC_in))

####Prediction for all models 
AUC_data<- tibble(GJAM=gjamDP$AUC_out, GJAM1=gjamDP1$AUC_out, PY1= gjamPY1$AUC_out , PY2= gjamPY2$AUC_out)
AUC_data_in<-  tibble(GJAM=gjamDP$AUC_in, GJAM1=gjamDP1$AUC_in, PY1= gjamPY1$AUC_in , PY2= gjamPY2$AUC_in)
AUC_data_cond<- tibble(GJAM=gjamDP$AUC_cond, GJAM1=gjamDP1$AUC_cond, PY1= gjamPY1$AUC_cond, PY2= gjamPY2$AUC_cond )
AUC_data$species<- colnames(Ydata_test)[1:ncol(Ydata_test)]
AUC_fin<- melt(AUC_data)

p2<-ggplot(data=AUC_fin)+geom_boxplot(aes(y=as.numeric(value),x=as.factor(variable),fill=as.factor(variable)))+
  scale_y_continuous(name="AUC")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM1","PY1","PY2"))+xlab("Models")+ theme_bw() 
p2

AUC_fin_cond_table<- as.data.frame(t(apply(AUC_data_cond,2,mean)))
AUC_fin_in_table<- as.data.frame(t(apply(AUC_data_in,2,mean)))
AUC_fin_table<- as.data.frame(t(apply(AUC_data[,1:4],2,mean)))
WAUC_fin_table<- as.data.frame(t(apply(AUC_data[,1:4],2,function(x) sum(x*p_w))))
names(AUC_fin_table)<- c("GJAM","GJAM1","PY1","PY2")
kable(cbind(data.frame( Measure = rbind("AUC")), rbind(AUC_fin_table)), format="pandoc", caption= "Prediction")

################################################################################################################

#Convergence of Sigma
#Sigma
df1<- tibble(ES= effectiveSize(mcmc(fit_gjam$chains$sgibbs[fit_gjam$modelList$burnin:fit_gjam$modelList$ng,])))
df2<- tibble( ES=effectiveSize(mcmc(fit_gjamDP1$chains$sgibbs[fit_gjamDP1$modelList$burnin:fit_gjamDP1$modelList$ng,]) ))
df3<- tibble(ES=effectiveSize(mcmc(fit_gjamPY1$chains$sgibbs[fit_gjamPY1$modelList$burnin:fit_gjamPY1$modelList$ng,])))
df4<- tibble(ES =effectiveSize(mcmc(fit_gjamPY2$chains$sgibbs[fit_gjamPY2$modelList$burnin:fit_gjamPY2$modelList$ng,])))
#pdf(file = "Plots/Effective_size_sigma.pdf", width= 8.27, height = 9.69)
rbind(df1 %>% mutate(var = "DP"),
      df2 %>%  mutate(var = "DP1"), 
      df3 %>% mutate(var = "PY1")
      #,df4 %>%  mutate(var = "PY2")
      ) %>% 
  ggplot(aes(ES, color = var, fill = var, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 
#dev.off()

df1<- tibble(ES= effectiveSize(mcmc(fit_gjam$chains$sgibbs[thin_5,])))
df2<- tibble( ES=effectiveSize(mcmc(fit_gjamDP1$chains$sgibbs[thin_5,]) ))
df3<- tibble(ES=effectiveSize(mcmc(fit_gjamPY1$chains$sgibbs[thin_5,])))
df4<- tibble(ES =effectiveSize(mcmc(fit_gjamPY2$chains$sgibbs[thin_5,])))
#pdf(file = "Plots/Effective_size_sigma.pdf", width= 8.27, height = 9.69)
rbind(df1 %>% mutate(var = "DP"),
      df2 %>%  mutate(var = "DP1"), 
      df3 %>% mutate(var = "PY1")
      #,df4 %>%  mutate(var = "PY2")
) %>% 
  ggplot(aes(ES, color = var, fill = var, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 
#dev.off()



thin_5<- seq(20001, 70000, by =5)
s1<- mcmc(fit_gjam$chains$sgibbs[thin_5,])
effectiveSize(s1)
##Convergence of the sigma_matrix_coefficients


#pdf(file = "Plots/Effective_size_beta.pdf", width= 8.27, height = 9.69)

dfb1<- tibble(ES= effectiveSize(mcmc(fit_gjam$chains$bgibbs)))
dfb2<- tibble( ES=effectiveSize(mcmc(fit_gjamDP1$chains$bgibbs) ))
dfb3<- tibble(ES=effectiveSize(mcmc(fit_gjamPY1$chains$bgibbs)))
dfb4<- tibble(ES =effectiveSize(mcmc(fit_gjamPY2$chains$bgibbs)))


rbind(dfb1 %>% mutate(var = "DP"),
      dfb2 %>%  mutate(var = "DP1"), 
      dfb3 %>% mutate(var = "PY1")
    #  , dfb4 %>%  mutate(var = "PY2")
      ) %>% 
  ggplot(aes(ES, color = var, fill = var, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 
#dev.off()

############################################################################################################

### Prior/Posterior and convergence for hyperparameters 
p_DP1 =tibble(it= 1: length(fit_gjamDP1$chains$alpha.DP_g),
       DP1= fit_gjamDP1$chains$alpha.DP_g) %>%
  ggplot(aes(x=it,y=DP1))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title=TeX(sprintf('Trace plot for DP1 parameter $\\alpha$')))+xlab("iterations")+ylab("Concentration parameter DP1") +theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15)) +theme(legend.text=element_text(size=15))



### convergence to alpha DP1

alpha_ch_1<- fit_gjamDP1$chains$alpha.DP_g[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng]
alpha_ch_2<- fit_gjamDP1_2$chains$alpha.DP_g[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng]

plot(1:50000,alpha_ch_2, type="l", col="red")
lines(1:50000,alpha_ch_1)
alpha_chains  = mcmc.list(mcmc(alpha_ch_1), mcmc(alpha_ch_2))
gelman.diag(alpha_chains)
gelman.plot(alpha_chains)
lattice::xyplot(alpha_chains, autoburinin=FALSE)
effectiveSize(mcmc(alpha_ch_1))

non_c_chains_1 <- fit_gjamDP1$chains$sgibbs[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng,]
non_c_chains_2 <- fit_gjamDP1_2$chains$sgibbs[(fit_gjamDP1_2$modelList$burnin+1):fit_gjamDP1_2$modelList$ng,]
non_c_chain  = mcmc.list(mcmc(non_c_chains_1), mcmc(non_c_chains_2))
x =gelman.diag(non_c_chain)
summary(x$psrf[,1])

non_c_chains_1 <- fit_gjamPY1$chains$sgibbs[(fit_gjamPY1$modelList$burnin+1):fit_gjamPY1$modelList$ng,]
non_c_chains_2 <- fit_gjamPY1_2$chains$sgibbs[(fit_gjamPY1_2$modelList$burnin+1):fit_gjamPY1_2$modelList$ng,]
non_c_chain  = mcmc.list(mcmc(non_c_chains_1), mcmc(non_c_chains_2))
x =gelman.diag(non_c_chain)
summary(x$psrf[,1])
hist(x$psrf[,1])
x$psrf[x$psrf[,1]>1.1,1]
names<- names(x$psrf[x$psrf[,1]>1.1,1])

##### Ploting chains with large psf
plots<- list()
ESS <- list()
for(key in names) {
  non_c_chains_1 <- fit_gjamPY1$chains$sgibbs[(fit_gjamPY1$modelList$burnin +1):fit_gjamPY1$modelList$ng,c(key)]
  non_c_chains_2 <- fit_gjamPY1_2$chains$sgibbs[(fit_gjamPY1_2$modelList$burnin +1):fit_gjamPY1_2$modelList$ng,c(key)]
  non_c_chain  = mcmc.list(mcmc(non_c_chains_1), mcmc(non_c_chains_2))
  x =gelman.diag(non_c_chain)
  x$psrf[1,]
  gelman.plot(non_c_chain)
  p = tibble(iterations = 1:(fit_gjamPY1_2$modelList$ng - fit_gjamPY1_2$modelList$burnin),
         chain_1 =non_c_chains_1,
         chain_2 = non_c_chains_2 ) %>%
    gather(Chains, trace, chain_1:chain_2)%>%
    ggplot(aes(x=iterations,y=trace,col=Chains))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
    labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()
  plots<- list.append(plots,p)
  ESS<- list.append(ESS, effectiveSize(mcmc(c(non_c_chains_2,non_c_chains_1))))
}

prow <- plot_grid(
     plots[[1]] + theme(legend.position='none'),
     plots[[2]]+ theme(legend.position='none'),
     plots[[3]]+ theme(legend.position='none'),
     plots[[4]]+ theme(legend.position='none'),
     plots[[5]]+ theme(legend.position='none'),
     plots[[6]]+ theme(legend.position='none'),
     plots[[7]]+ theme(legend.position='none'),
     plots[[8]]+ theme(legend.position='none'),
    nrow = 3, ncol=3)

legend_b <- get_legend(plots[[8]]+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))

pdf
plot(p)
#pdf(file = "Plots/PY1_Rhat_large.pdf", width= 8.27, height = 9.69)
#plot(p)
#dev.off()


plots2<- list()
ESS2 <- list()
for(key in names) {
  non_c_chains_1 <- fit_gjamPY1$chains$sgibbs[(fit_gjamPY1$modelList$burnin +1):fit_gjamPY1$modelList$ng,c(key)]
  non_c_chains_2 <- fit_gjamPY1_2$chains$sgibbs[(fit_gjamPY1_2$modelList$burnin +1):fit_gjamPY1_2$modelList$ng,c(key)]
  non_c_chain  = mcmc.list(mcmc(non_c_chains_1), mcmc(non_c_chains_2))
  x =gelman.diag(non_c_chain)
  x$psrf[1,]
  gelman.plot(non_c_chain)
  plots2<- list.append(plots2, lattice::xyplot(non_c_chain, autoburinin=FALSE))
  ESS2<- list.append(ESS2, effectiveSize(mcmc(c(non_c_chains_2,non_c_chains_1))))
}


prow <- plot_grid(
  plots[[1]] + theme(legend.position='none'),
  plots[[2]]+ theme(legend.position='none'),
  plots[[3]]+ theme(legend.position='none'),
  plots[[4]]+ theme(legend.position='none'),
  plots[[5]]+ theme(legend.position='none'),
  plots[[6]]+ theme(legend.position='none'),
  nrow = 3, ncol=2)

legend_b <- get_legend(plots[[6]]+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))




lattice::xyplot(non_c_chain, autoburinin=FALSE)
#r-2_N-64  r-2_N-97 r-2_N-108  r-3_N-64   r-4_N-1  r-4_N-64  r-4_N-97 r-4_N-102  r-5_N-64  r-5_N-95 
#r-1_N-64  r-2_N-47  r-2_N-64  r-2_N-68  r-2_N-97 r-2_N-108  r-4_N-97 r-4_N-102 
non_c_chains_1 <- fit_gjamPY1$chains$sgibbs[(fit_gjamPY1$modelList$burnin +1):fit_gjamPY1$modelList$ng,c("r-1_N-642")]
non_c_chains_2 <- fit_gjamPY1_2$chains$sgibbs[(fit_gjamPY1_2$modelList$burnin +1):fit_gjamPY1_2$modelList$ng,c("r-1_N-64")]
non_c_chain  = mcmc.list(mcmc(non_c_chains_1), mcmc(non_c_chains_2))
x =gelman.diag(non_c_chain)
x$psrf[1,]
gelman.plot(non_c_chain)
#gelman.plot(alpha_chains)
lattice::xyplot(non_c_chain, autoburinin=FALSE)
effectiveSize(mcmc(c(fit_gjam$chains$sgibbs)))
summary(non_c_chains_1)
summary(x$psrf[,1])

gelman.plot(non_c_chain)
lattice::xyplot(non_c_chain, autoburinin=FALSE)


sigma_chain_1 <- mcmc(fit_gjam$chains$sgibbs[(fit_gjamPY1$modelList$burnin+10000):fit_gjamPY1$modelList$ng,])
sigma_chain_2 <- mcmc(fit_gjam_2$chains$sgibbs[(fit_gjamPY1_2$modelList$burnin+10000):fit_gjamPY1_2$modelList$ng,])
sigma_chains  = mcmc.list(sigma_chain_1, sigma_chain_2)
gelman_sigma <- gelman.diag(sigma_chains)
summary(gelman_sigma$psrf[,1])
hist(gelman_sigma$psrf[,1])


sigma_chain_1 <- mcmc(fit_gjamDP1$chains$sgibbs[(fit_gjamDP1$modelList$burnin):fit_gjamPY1$modelList$ng,])
sigma_chain_2 <- mcmc(fit_gjamDP1_2$chains$sgibbs[(fit_gjamDP1_2$modelList$burnin):fit_gjamPY1_2$modelList$ng,])
sigma_chains  = mcmc.list(sigma_chain_1, sigma_chain_2)
gelman_sigma <- gelman.diag(sigma_chains)
summary(gelman_sigma$psrf[,1])
hist(gelman_sigma$psrf[,1])


sigma_chain_1 <- mcmc(fit_gjamPY12$chains$sgibbs[(fit_gjamPY1$modelList$burnin):fit_gjamPY1$modelList$ng,])
sigma_chain_2 <- mcmc(fit_gjamPY12_2$chains$sgibbs[(fit_gjamPY1_2$modelList$burnin):fit_gjamPY1_2$modelList$ng,])
sigma_chains  = mcmc.list(sigma_chain_1, sigma_chain_2)
gelman_sigma <- gelman.diag(sigma_chains)
summary(gelman_sigma$psrf[,1])
hist(gelman_sigma$psrf[,1])


sigma_chain_1 <- mcmc(fit_gjamDP1$chains$sgibbs[,c("r-1_N-2")])
sigma_chain_2 <- mcmc(fit_gjamDP1_2$chains$sgibbs[,c("r-1_N-2")])
sigma_chains  = mcmc.list(sigma_chain_1, sigma_chain_2)
gelman_sigma <- gelman.diag(sigma_chains)
summary(gelman_sigma$psrf[,1])
gelman.plot(sigma_chain_1)
lattice::xyplot(sigma_chain_1, autoburinin=FALSE)
effectiveSize(mcmc(gelman_sigma))
geweke.diag(sigma_chain_1)
geweke.plot(sigma_chain_1)
gelman_sigma$psrf[which(gelman_sigma$psrf[,1]>1.1),1]
lattice::xyplot(sigma_chains, autoburinin=FALSE)
#tibble(ES =effectiveSize(mcmc(fit_gjamDP2$chains$alpha.DP_g)))

#pdf(file = "Plots/Posterior_for_DP_parameters.pdf", width= 8.27, height = 9.69)
plot(p_DP1)
#dev.off()
# 
# prow <- plot_grid(
#   p_DP1 + theme(legend.position='none'),
#   p_PY2 + theme(legend.position='none'),
#   p_PY2_d + theme(legend.position='none'),
#   nrow = 2
# )
#legend_b <- get_legend(p_DP2+theme(legend.position ='top'))
#p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))
#plot(p)
#dev.off()

#df1<- tibble(ES= effectiveSize(mcmc(fit_gjamDP1$chains$alpha.DP_g)))
#df2<- tibble( ES=effectiveSize(mcmc(fit_gjamPY2$chains$alpha.PY_g) ))
#df3<- tibble(ES=effectiveSize(mcmc(fit_gjamPY2$chains$discount.PY_g)))
#rbind(df1 %>% mutate(var = "DP2"),
#      df2 %>%  mutate(var = "PY2_a"), 
#      df3 %>% mutate(var = "PY2_d")) %>% 
#  ggplot(aes(ES, color = var, fill = var, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 




#### Add prior/posterior plot
prior_nu1 <- fit_gjamDP1$modelList$reductList$otherpar$shape
prior_nu2 <- fit_gjamDP1$modelList$reductList$otherpar$rate
alpha_vec<- rgamma(50000, prior_nu1,prior_nu2)
x<- sapply(alpha_vec, functionDPM,n=112,N=112)
mean(x)

###Posterior for alpha parameter DP2

aDP1 =tibble(Prior =alpha_vec,
       Posterior = fit_gjamDP1$chains$alpha.DP_g[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng])%>%
  gather(Distribution, density, Prior:Posterior)%>%
 ggplot(aes(x=density, fill=Distribution)) +
  geom_density(adjust = 1, alpha=0.7)+ggtitle(TeX(sprintf('Prior and posterior distribution for D1 $\\alpha$'))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))


#prior_nu1 <- fit_gjamPY2$modelList$reductList$otherpar$shape
#prior_nu2 <- fit_gjamPY2$modelList$reductList$otherpar$rate
#alpha_vecPY<- rgamma(18000, prior_nu1,prior_nu2)
#x<- sapply(alpha_vecPY, functionPY,n=112, sigma_py=0.25)
#mean(x)


#pdf(file = "Plots/Posterior_for_DP1_parameters.pdf", width= 8.27, height = 9.69)
#plot(aDP1)
#dev.off()


##### Posterior distribution for the number of clusters


############################Trace plot##########################################################################
#pdf(file = "Plots/Trace_plot_partitions.pdf", width= 8.27, height = 9.69)

tibble(it= 1: length(apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))),
              DP= apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x))),
              DP1 =apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x))),
              PY1=apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x))),
             PY2=apply(fit_gjamPY2$chains$kgibbs,1,function(x) length(unique(x))) 
       ) %>%
gather(Model, trace, DP:PY2)%>%
 ggplot(aes(x=it,y=trace,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
 labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
 theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
 theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
       axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
       plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))

#dev.off()

##### Truncation weights ############################################################

last_pk<- round(mean(fit_gjamDP1$chains$pk_g[-c(1:fit_gjamDP1$modelList$burnin),fit_gjamDP1$modelList$reductList$N]),5)
df_weights <- tibble(pw =apply(fit_gjamDP1$chains$pk_g[-c(1:fit_gjamDP1$modelList$burnin),],2,mean),
                     tr= 1:ncol(fit_gjamDP1$chains$pk_g))
pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
  geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
  geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Mean weights for DP1, p_last= ",last_pk), 
                                                                                                    caption=paste0("Number of iterations: ",fit_gjamDP1$modelList$ng," burnin: ",fit_gjamDP1$modelList$burnin))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
pl_weigths


last_pk<- round(mean(fit_gjamPY1$chains$pk_g[-c(1:fit_gjamPY1$modelList$burnin),fit_gjamPY1$modelList$reductList$N]),5)
df_weights <- tibble(pw =apply(fit_gjamPY1$chains$pk_g[-c(1:fit_gjamPY1$modelList$burnin),],2,mean),
                     tr= 1:ncol(fit_gjamPY1$chains$pk_g))
pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
  geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
  geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Mean weights for PY1, p_last= ",last_pk), 
                                                                                                    caption=paste0("Number of iterations: ",fit_gjamPY1$modelList$ng," burnin: ",fit_gjamPY1$modelList$burnin))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
pl_weigths



last_pk<- round(mean(fit_gjamPY2$chains$pk_g[-c(1:fit_gjamPY2$modelList$burnin),fit_gjamPY2$modelList$reductList$N]),5)
df_weights <- tibble(pw =apply(fit_gjamPY2$chains$pk_g[-c(1:fit_gjamPY2$modelList$burnin),],2,mean),
                     tr= 1:ncol(fit_gjamPY2$chains$pk_g))
pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
  geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
  geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Mean weights for PY1, p_last= ",last_pk), 
                                                                                                    caption=paste0("Number of iterations: ",fit_gjamPY2$modelList$ng," burnin: ",fit_gjamPY2$modelList$burnin))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
pl_weigths



##################################Clustering #####################################


relabel_clust<- function(Mat ){
  Mat_new<- matrix(NA, nrow=nrow(Mat),ncol=ncol(Mat))
  for(i in 1:nrow(Mat) ){
    label<-unique(Mat[i,])
    change<- label[which(label>ncol(Mat))]
    all_n<- 1:ncol(Mat)
    admis<- all_n[which(!(all_n %in% label))]
    Mat_new[i,]<-Mat[i,]
    if(length(change)>0){
      for(k in 1:length(change)){
        old_row<- Mat_new[i,]
        rep <-old_row==change[k]
        Mat_new[i,] <- replace(Mat_new[i,], rep, admis[k])
      }
    }
  }
  return(Mat_new)
}




gjamClust<- function(model){
  if("other" %in% colnames(model$inputs$y)){
    sp_num<- ncol(model$inputs$y)-1
  } 
   MatDP<- relabel_clust(model$chains$kgibbs[(model$modelList$burnin+1):model$modelList$ng,1:sp_num])
  CM_DP<-  comp.psm(MatDP)
  mbind_DP.ext <- minbinder.ext(CM_DP,MatDP, method="all",include.greedy = TRUE)
  vi_DP.ext <- minVI(CM_DP,MatDP, method="all", include.greedy = TRUE)
  #start.cl =  True_clust$K_n
  return(list(Bind_est=mbind_DP.ext$cl[1,], VI_est = vi_DP.ext$cl[1,]))
}



gjamClust_full_test<- function(model, K=K_prior, true_clust =True_clust$K_n ){
  if("other" %in% colnames(model$inputs$y)){
    sp_num<- ncol(model$inputs$y)-1
  } 
  LF_value<- c()
  VI_list<- list()
  MatDP<- relabel_clust(model$chains$kgibbs[(model$modelList$burnin+1):model$modelList$ng,1:sp_num])
  CM_DP<-  comp.psm(MatDP)
  start_list<- list(1:(dim(CM_DP)[1]), sample(1:K,dim(CM_DP)[1],replace = TRUE), true_clust)
  #mbind_DP.ext <- minbinder.ext(CM_DP,MatDP, method="all",include.greedy = TRUE)
  for (i in 1:length(start_list)){
    vi_DP.ext <- minVI(CM_DP,MatDP, method="all", include.greedy = TRUE, start.cl = start_list[[i]])
    VI_list<- list(VI_list, vi_DP.ext$cl[1,])
    LF_value<- c(LF_value, vi_DP.ext$value[1])
  }
  #start.cl =  True_clust$K_n
  return(list(VI_loss_value = LF_value, VI_est = VI_list))
}



#DP_clust_full<- gjamClust_full_test(model= fit_gjam)
#DP1_clust_full<- gjamClust_full_test(model= fit_gjamDP1)
#PY1_clust_full<- gjamClust_full_test(model= fit_gjamPY1)
#PY2_clust_full<- gjamClust_full_test(model= fit_gjamPY2)

if("other" %in% colnames(model$inputs$y)){
  sp_num<- ncol(model$inputs$y)-1
} 
MatDP<- relabel_clust(chains_2[,1:sp_num])
CM_DP<-  comp.psm(MatDP)
mbind_DP.ext <- minbinder.ext(CM_DP,MatDP, method="all",include.greedy = TRUE)
vi_DP.ext <- minVI(CM_DP,MatDP, method="all", include.greedy = TRUE)

non_c_chains_1 <- fit_gjamPY1$chains$kgibbs[(fit_gjamPY1$modelList$burnin +1):fit_gjamPY1$modelList$ng,1:111]
non_c_chains_2 <- fit_gjamPY1_2$chains$kgibbs[(fit_gjamPY12_2$modelList$burnin +1):fit_gjamPY12_2$modelList$ng,1:111]


DP_grEPL<- MinimiseEPL(MatDP, pars =list("VI"))

load("PCA_analysis/r5/Clusters_modells_1.Rdata")
DP_grEPL$decision
arandi(DP_grEPL$decision, Cluster_models_1$ClustPY1)

DP_clust<- gjamClust(model= fit_gjam)
DP1_clust<- gjamClust(model= fit_gjamDP1)
PY1_clust<- gjamClust(model= fit_gjamPY1)
PY2_clust<- gjamClust(model= fit_gjamPY2)



DP_clust_2<- gjamClust(model= fit_gjam_2)
DP1_clust_2<- gjamClust(model= fit_gjamDP1_2)
PY1_clust_2<- gjamClust(model= fit_gjamPY1_2)

arandi(DP_clust$VI_est,DP_clust_2$VI_est)
arandi(DP1_clust$VI_est,DP1_clust_2$VI_est)
arandi(PY1_clust$VI_est,PY1_clust_2$VI_est)

#### Compare the obtained estimates with the PFG clusters

BI_fin_table_K<- tibble(GJAM = length(unique(DP_clust$Bind_est)),GJAM1=length(unique(DP1_clust$Bind_est)),PY1=length(unique(PY1_clust$Bind_est)),PY2=length(unique(PY2_clust$Bind_est)))
BI_fin_table_VIdist<- tibble(GJAM = vi.dist(DP_clust$Bind_est, True_clust$K_n),GJAM1 = vi.dist(DP1_clust$Bind_est, True_clust$K_n),PY1= vi.dist(PY1_clust$Bind_est, True_clust$K_n),PY2= vi.dist(PY2_clust$Bind_est, True_clust$K_n))
BI_fin_table_ARdist<- tibble(GJAM = arandi(DP_clust$Bind_est, True_clust$K_n),GJAM1 = arandi(DP1_clust$Bind_est, True_clust$K_n),PY1 = arandi(PY1_clust$Bind_est, True_clust$K_n),PY2 = arandi(PY2_clust$Bind_est, True_clust$K_n))

VI_fin_table_K<-tibble(GJAM = length(unique(DP_clust$VI_est)),GJAM1=length(unique(DP1_clust$VI_est)),PY1=length(unique(PY1_clust$VI_est)),PY2=length(unique(PY2_clust$VI_est)))
VI_fin_table_VIdist<-  tibble(GJAM = vi.dist(DP_clust$VI_est, True_clust$K_n),GJAM1 = vi.dist(DP1_clust$VI_est, True_clust$K_n),PY1= vi.dist(PY1_clust$VI_est, True_clust$K_n),PY2= vi.dist(PY2_clust$VI_est, True_clust$K_n))
VI_fin_table_ARdist<- tibble(GJAM = arandi(DP_clust$VI_est, True_clust$K_n),GJAM1 = arandi(DP1_clust$VI_est, True_clust$K_n),PY1 = arandi(PY1_clust$VI_est, True_clust$K_n),PY2 = arandi(PY2_clust$VI_est, True_clust$K_n))





#MatDP<- fit_gjam$chains$kgibbs[(fit_gjam$modelList$burnin+1):fit_gjam$modelList$ng,1:112]
#DP_grEPL<- MinimiseEPL(MatDP, pars = list(Kup = 1, loss_type ="VI"))

#
#DP_clust2<- gjamClust2(model= fit_gjam, pars =list(decision_init=True_clust$K_n, "VI"))



gjamClust2_full<- function(model, K = K_prior,  true_clust =True_clust$K_n){
  if("other" %in% colnames(model$inputs$y)){
    sp_num<- ncol(model$inputs$y)-1
  } 
  LF_value<- c()
  VI_list<- list()
  MatDP<- relabel_clust(model$chains$kgibbs[(model$modelList$burnin+1):model$modelList$ng,1:sp_num])
  start_list<- list(1:(dim(MatDP)[2]), sample(1:K,dim(MatDP)[2],replace = TRUE), true_clust)
  for (i in 1:length(start_list)){
    DP_grEPL<- MinimiseEPL(MatDP, pars = list(decision_init=start_list[[i]], loss_type = "VI"))
    VI_list<- list.append(VI_list, DP_grEPL$decision)
    LF_value<- c(LF_value, DP_grEPL$EPL)
  }
  #start.cl =  True_clust$K_n
   return(list( VI_est = VI_list, EPL_value = LF_value))
}


gjamClust2<- function(model, pars =list("VI")){
  if("other" %in% colnames(model$inputs$y)){
    sp_num<- ncol(model$inputs$y)-1
  } 
  MatDP<- relabel_clust(model$chains$kgibbs[(model$modelList$burnin+1):model$modelList$ng,1:sp_num])
  DP_grEPL<- MinimiseEPL(MatDP, pars =pars)
  return(list( VI_est = DP_grEPL$decision, EPL_value =  DP_grEPL$EPL))
}

#MatDP<- fit_gjam$chains$kgibbs[(fit_gjam$modelList$burnin+1):fit_gjam$modelList$ng,1:112]
#DP_grEPL<- MinimiseEPL(MatDP, pars = list(Kup = 1, loss_type ="VI"))

#
#DP_clust2<- gjamClust2(model= fit_gjam, pars =list(decision_init=True_clust$K_n, "VI"))

if("other" %in% colnames(model$inputs$y)){
  sp_num<- ncol(model$inputs$y)-1
} 
MatDP<- relabel_clust(chains_2[1:sp_num])
DP_grEPL<- MinimiseEPL(MatDP, pars =pars)
return(list( VI_est = DP_grEPL$decision, EPL_value =  DP_grEPL$EPL))
chains_2

DP_clust2_full<- gjamClust2_full(model= fit_gjam)
DP1_clust2_full<- gjamClust2_full(model= fit_gjamDP1)
PY1_clust2_full<- gjamClust2_full(model= fit_gjamPY1)
PY2_clust2_full<- gjamClust2_full(model= fit_gjamPY2)


DP_clust2<- gjamClust2(model= fit_gjam)
DP1_clust2<- gjamClust2(model= fit_gjamDP1)
PY1_clust2<- gjamClust2(model= fit_gjamPY1)
PY2_clust2<- gjamClust2(model= fit_gjamPY2)
PY12_clust2<- gjamClust2(model= fit_gjamPY12)


DP_clust2_2<- gjamClust2(model= fit_gjam_2)
DP1_clust2_2<- gjamClust2(model= fit_gjamDP1_2)
PY1_clust2_2<- gjamClust2(model= fit_gjamPY1_2)
PY12_clust2_2<- gjamClust2(model= fit_gjamPY12_2)

arandi(DP_clust2$VI_est,DP_clust2_2$VI_est)
arandi(DP1_clust2$VI_est,DP1_clust2_2$VI_est)
arandi(PY1_clust2$VI_est,PY1_clust2_2$VI_est)
arandi(PY12_clust2$VI_est,PY12_clust2_2$VI_est)
arandi(PY1_clust2$VI_est, PY12_clust2$VI_est)
arandi(PY1_clust2_2$VI_est, PY12_clust2_2$VI_est)

VI_fin_table_K2<-tibble(GJAM = length(unique(DP_clust2$VI_est)),GJAM1=length(unique(DP1_clust2$VI_est)),PY1=length(unique(PY1_clust2$VI_est)),PY2=length(unique(PY2_clust2$VI_est)))
VI_fin_table_VIdist2<-  tibble(GJAM = vi.dist(DP_clust2$VI_est, True_clust$K_n),GJAM1 = vi.dist(DP1_clust2$VI_est, True_clust$K_n),PY1= vi.dist(PY1_clust2$VI_est, True_clust$K_n),PY2= vi.dist(PY2_clust2$VI_est, True_clust$K_n))
VI_fin_table_ARdist2<- tibble(GJAM = arandi(DP_clust2$VI_est, True_clust$K_n),GJAM1 = arandi(DP1_clust2$VI_est, True_clust$K_n),PY1 = arandi(PY1_clust2$VI_est, True_clust$K_n),PY2 = arandi(PY2_clust2$VI_est, True_clust$K_n))




###############################################################################################################
# sp_num<- ncol(Ydata)-1
# MatDP<- relabel_clust(fit_gjam$chains$kgibbs[(fit_gjam$modelList$burnin+1):fit_gjam$modelList$ng,1:sp_num])
# MatDP1<- relabel_clust(fit_gjamDP1$chains$kgibbs[(fit_gjamDP1$modelList$burnin+1):fit_gjamDP1$modelList$ng,1:sp_num])
# MatPY1<- relabel_clust(fit_gjamPY1$chains$kgibbs[(fit_gjamPY1$modelList$burnin+1):fit_gjamPY1$modelList$ng,1:sp_num])
# MatPY2<- relabel_clust(fit_gjamPY2$chains$kgibbs[(fit_gjamPY2$modelList$burnin+1):fit_gjamPY2$modelList$ng,1:sp_num])

############ greedy EPL ##### another algorithm 
# DP_grEPL<- MinimiseEPL(MatDP, pars = list())
# length(unique(DP_grEPL$decision))
# arandi(DP_grEPL$decision,DP_clust$VI_est )
# DP1_grEPL<- MinimiseEPL(MatDP1, pars = list(loss_type="VI"))
# length(unique(DP1_grEPL$decision))
# 
# PY1_grEPL<- MinimiseEPL(MatPY1, pars = list(Kup=5, loss_type="VI"))
# length(unique(PY1_grEPL$decision))
# 
# PY2_grEPL<- MinimiseEPL(MatPY2, pars = list(loss_type="VI"))
# length(unique(PY2_grEPL$decision))
# 
# arandi(PY1_grEPL$decision,PY1_clust$VI_est)
# arandi(DP1_grEPL$decision,DP1_clust$VI_est)
# arandi(DP_grEPL$decision,DP1_grEPL$decision)
# arandi(DP_grEPL$decision,DP1_grEPL$decision)

# arandi(DP1_clust$VI_est,PY1_clust$VI_est)
# 
# arandi(PY1_grEPL$decision,DP1_grEPL$decision)
# arandi(PY1_grEPL$decision,DP1_grEPL$decision)
# 
##### confidence intervals  ############################################################

#DP_cb = credibleball(DP_clust$VI_est, MatDP, c.dist = c("VI","Binder"), alpha = 0.05)

############################################################################################################
# 
# A=load_object("PCA_analysis/r5/Clusters_modells_1.Rdata")
# c1<-PY1_clust$VI_est
# c2<-DP1_clust$VI_est
# 
# arandi(A$ClustPY1,c1 )
# arandi(A$ClustDP1,c2)
# arandi(c1,c2)
# 
# arandi(A$ClustPY1,PY1_clust2$VI_est )
# arandi(A$ClustDP2,DP1_clust2$VI_est )
# 
 # arandi(PY1_clust2$VI_est,PY1_clust$VI_est)
 # arandi(PY1_clust2$VI_est,DP1_clust2$VI_est )
 # arandi(PY1_clust$VI_est,DP1_clust$VI_est )
 # arandi(DP1_clust2$VI_est,DP1_clust$VI_est )
 # 
# 
#summary(DP_cb)
#The credible ball characterizes the uncertainty in the clustering esitmate.
# 
# DP_cb = credibleball(DP_clust$VI_est, MatDP2, c.dist = c("VI"), alpha = 0.05)
# DP1_cb = credibleball(DP1_clust$VI_est, MatDP1, c.dist = c("VI"), alpha = 0.05)
# PY1_cb = credibleball(PY1_clust$VI_est, MatPY1, c.dist = c("VI"), alpha = 0.05)
# PY2_cb = credibleball(PY2_clust$VI_est, MatPY2, c.dist = c("VI"), alpha = 0.05)
# 

#### Cluster names
# Cluster_models1<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustDP=DP_clust$VI_est,
#                          ClustDP1=DP1_clust$VI_est,ClustPY1=PY1_clust$VI_est,ClustPY2=PY2_clust$VI_est, PFG=True_clust$K_n)
Cluster_models_1<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustDP=DP_clust$VI_est,
                          ClustDP1=DP1_clust$VI_est,ClustPY1=PY1_clust$VI_est,PFG=True_clust$K_n)

save(Cluster_models_1, file =  paste0(folderpath,"Clusters_modells_1.Rdata"))
Cluster_models_2<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustDP=DP_clust2$VI_est,
                           ClustDP1=DP1_clust2$VI_est,ClustPY1=PY1_clust2$VI_est,PFG=True_clust$K_n)

save(Cluster_models_2, file =  paste0(folderpath,"Clusters_modells_2.Rdata"))

##################################Covariance matrix######################################################################
#### Add cluster labels 
Colnames_Y_clust<- merge(Colnames_Y, Cluster_models_2, by ="CODE_CBNA")

### Covariance matrix for the mean 
#pdf(file = "Plots/Correlation_matrix_DP.pdf", width= 8.27, height = 9.69)
MDP= fit_gjam$parameters$corMu
Colnames_Y_clust$Sp_name_DP <- paste(Colnames_Y_clust$ClustDP, Colnames_Y_clust$species,sep="_")
rownames(MDP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP"],"other")
colnames(MDP)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP"],"other")
colors_vir=viridis(length(unique(Colnames_Y_clust$ClustDP))+1, option = "magma")
LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_DP),"ClustDP"],length(unique(Colnames_Y_clust$ClustDP))+1), function(x) colors_vir[x])
cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))

corrplot(MDP, diag = FALSE, order = "hclust", tl.cex = 0.45,tl.srt=45, method = "color", tl.col=LabelCol,col=cols(200),
         type = "full", title= "Correlation for the DP model (original)", mar=c(0,0,1,0))

corrplot(MDP, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45,  tl.col=LabelCol,
         method = "color",col=cols(200), type = "full",title= "Correlation for the DP model (DP groups)", mar=c(0,0,1,0))
#dev.off()

#pdf(file = "Plots/Correlation_matrix_DP2.pdf", width= 8.27, height = 9.69)
MDP1= fit_gjamDP1$parameters$corMu
Colnames_Y_clust$Sp_name_DP1 <- paste(Colnames_Y_clust$ClustDP1, Colnames_Y_clust$species,sep="_")
rownames(MDP1)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP1"],"other")
colnames(MDP1)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP1"],"other")
colors_vir=viridis(length(unique(Colnames_Y_clust$ClustDP1))+1, option = "magma")
LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_DP1),"ClustDP1"],15), function(x) colors_vir[x])

corrplot(MDP1, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.45,tl.srt=45, method = "color",col=cols(200),
         type = "lower", title= "Correlation for the DP1 model (original)", mar=c(0,0,1,0))

corrplot(MDP1, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45, tl.col=LabelCol,method = "color",col=cols(200),
         type = "full", title= "Correlation for the DP1 model (DP2 groups)", mar=c(0,0,1,0))

#dev.off()

#pdf(file = "Plots/Correlation_matrix_PY1.pdf", width= 8.27, height = 9.69)
MPY1= fit_gjamPY1$parameters$corMu
Colnames_Y_clust$Sp_name_PY1 <- paste(Colnames_Y_clust$ClustPY1, Colnames_Y_clust$species,sep="_")
rownames(MPY1)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_PY1"],"other")
colnames(MPY1)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_PY1"],"other")
colors_vir=viridis(length(unique(Colnames_Y_clust$ClustPY1))+1, option = "magma")
LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_PY1),"ClustPY1"],19), function(x) colors_vir[x])

corrplot(MPY1, diag = FALSE, order = "hclust", tl.cex = 0.45,tl.srt=45, method = "color",col=cols(200),
         type = "full", title= "Correlation for the PY1 model (hclust)", mar=c(0,0,1,0))

corrplot(MPY1, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45, tl.col=LabelCol,method = "color",col=cols(200),
         type = "full", title= "Correlation for the PY1 model (PY1 groups)", mar=c(0,0,1,0))

#dev.off()

#pdf(file = "Plots/Correlation_matrix_PY2.pdf", width= 8.27, height = 9.69)
# MPY2= fit_gjamPY2$parameters$corMu
# Colnames_Y_clust$Sp_name_PY2 <- paste(Colnames_Y_clust$ClustPY1, Colnames_Y_clust$species,sep="_")
# rownames(MPY2)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_PY2"],"other")
# colnames(MPY2)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_PY2"],"other")
# colors_vir=viridis(length(unique(Colnames_Y_clust$ClustPY2))+1, option = "magma")
# LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_PY2),"ClustPY2"],length(unique(Colnames_Y_clust$ClustPY2))+1), function(x) colors_vir[x])
# 
# corrplot(MPY1, diag = FALSE, order = "hclust", tl.cex = 0.45,tl.srt=45, method = "color",col=cols(200),
#          type = "full", title= "Correlation for the PY2 model (original)", mar=c(0,0,1,0))
# 
# corrplot(MPY1, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45, tl.col=LabelCol,method = "color",col=cols(200),
#          type = "full", title= "Correlation for the PY2 model (PY2 groups)", mar=c(0,0,1,0))

#dev.off()
###### Look at the graph  

##################################Covariance matrix using the  MCMC samples#####################################################################

##Covariance matrix
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


cov_matrix<- function(fit, burn_period,iterations){
  sgibbs<-fit$chains$sgibbs[burn_period:iterations,]
  sigErrGibbs<-fit$chains$sigErrGibbs[burn_period:iterations]
  kgibbs<-fit$chains$kgibbs[burn_period:iterations,]
  sigma<-invsigma<-array(NA,dim=c(S,S,iterations-burn_period))
  N<-fit$modelList$reductList$N
  r<-fit$modelList$reductList$r
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
  ## Inverse
  invsigma_mean<-apply(invsigma,c(1,2),mean) 
  invsigma_q05<-apply(invsigma,c(1,2),quantile,0.05) 
  invsigma_q95<-apply(invsigma,c(1,2),quantile,0.95) 
  INVSigma_sign<--cov2cor(sigma_mean*(!(invsigma_q95>0 & invsigma_q05<0)))
  return(list(Sigma =Sigma_sign, InvS=  INVSigma_sign, S_mean=sigma_mean, IS_mean=invsigma_mean ))
}


A= cov_matrix(fit=fit_gjam, burn_period=fit_gjam$modelList$burnin,iterations=fit_gjam$modelList$ng)


cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))

gcols = colorRampPalette(c( "White", "White", "Black"))
corrplot(A$InvS, diag = FALSE, order = "hclust",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
corrplot(A$IS_mean, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")

A= cov_matrix(fit=fit_gjamDP1, burn_period=fit_gjamDP2$modelList$burnin +1,iterations=fit_gjamDP1$modelList$ng)


Mat<- cov2cor(A$S_mean)

Mat[abs(Mat)<0.2] <- 0
Graph_pcor <- qgraph(Mat,graph="pcor", layout="spring")
Graph_pcor
# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix( Mat, weighted=T, mode="undirected", diag=F)
plot(network)
library(corrr)
network_plot(Mat)

library(RColorBrewer)
coul <- brewer.pal(nlevels(as.factor(mtcars$cyl)), "Set2")

# Map the color to cylinders
my_color <- coul[as.numeric(as.factor(mtcars$cyl))]

# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)
plot(network, 
     vertex.size=12,
     vertex.color=my_color, 
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent"
)

# title and legend
text(0,0,"mtcars network",col="white", cex=1.5)
legend(x=-0.2, y=-0.12, 
       legend=paste( levels(as.factor(mtcars$cyl)), " cylinders", sep=""), 
       col = coul , 
       bty = "n", pch=20 , pt.cex = 2, cex = 1,
       text.col="white" , horiz = F)library(RColorBrewer)
coul <- brewer.pal(nlevels(as.factor(mtcars$cyl)), "Set2")

# Map the color to cylinders
my_color <- coul[as.numeric(as.factor(mtcars$cyl))]

# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)
plot(network, 
     vertex.size=12,
     vertex.label.cex=0.7,
     vertex.label.color="white",
     vertex.frame.color="transparent"
)

# title and legend
text(0,0,"mtcars network",col="white", cex=1.5)
legend(x=-0.2, y=-0.12, 
       legend=paste( levels(as.factor(mtcars$cyl)), " cylinders", sep=""), 
       col = coul , 
       bty = "n", pch=20 , pt.cex = 2, cex = 1,
       text.col="white" , horiz = F)

########################################################################################

#### Final Table
form<-c(formula)
Fin_all<-as.data.frame(matrix(NA,nrow=9,ncol=8))
names(Fin_all)<- c("Parameter","GJAM","GJAM1","PY1","r", "iter", "burn","formula")
Fin_all$iter<- fit_gjam$modelList$ng
Fin_all$burn<- fit_gjam$modelList$burnin
Fin_all$r<-fit_gjam$modelList$reductList$r
Fin_all$formula<-as.character(form)
Fin_all[1,1]<- "DIC"
Fin_all[1,2:4]<- c(fit_gjam$fit$DIC,fit_gjamDP1$fit$DIC,fit_gjamPY1$fit$DIC)/10000
Fin_all[2,1]<- "mean AUC"
Fin_all[2,2:4]<- AUC_fin_table[,1:3]
Fin_all[3,1]<- "mean WAUC"
Fin_all[3,2:4]<- WAUC_fin_table[,1:3]
Fin_all[4,1]<- "AUC in"
Fin_all[4,2:4]<- AUC_fin_in_table[,1:3]
Fin_all[5,1]<- "AR dist VI loss 2"
Fin_all[5,2:4]<- VI_fin_table_ARdist2[,1:3]
Fin_all[6,1]<- "AR dist VI loss"
Fin_all[6,2:4]<- VI_fin_table_ARdist[,1:3]
Fin_all[7,1]<- "mean K"
Fin_all[7,2:4]<- c(mean(apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjam$modelList$burnin:fit_gjam$modelList$ng]),
                    mean(apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamDP1$modelList$burnin:fit_gjamDP1$modelList$ng]),
                    mean(apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamPY1$modelList$burnin:fit_gjamPY1$modelList$ng]))
Fin_all[8,1]<- "K VI2"
Fin_all[8,2:4]<- VI_fin_table_K2[,1:3]
Fin_all[9,1]<- "K VI"
Fin_all[9,2:4]<- VI_fin_table_K[,1:3]
Fin_all[,2:4]<- round(Fin_all[,2:5], 3)
#save(Fin_all, file = paste0(folderpath,"Fin_tab_r5.Rdata"))
#library("xlsx")
# Write the first data set in a new workbook
#write.xlsx(Fin_all, file = "Final_table.xlsx")
#########################################################################################


#### Final Table
form<-c(formula)
Fin_all<-as.data.frame(matrix(NA,nrow=9,ncol=8))
names(Fin_all)<- c("Parameter","GJAM","GJAM1","PY1","PY2","r", "iter", "burn","formula")
Fin_all$iter<- fit_gjam$modelList$ng
Fin_all$burn<- fit_gjam$modelList$burnin
Fin_all$r<-fit_gjam$modelList$reductList$r
Fin_all$formula<-as.character(form)
Fin_all[1,1]<- "DIC"
Fin_all[1,2:5]<- c(fit_gjam$fit$DIC,fit_gjamDP1$fit$DIC,fit_gjamPY1$fit$DIC,fit_gjamPY2$fit$DIC)/10000
Fin_all[2,1]<- "mean AUC"
Fin_all[2,2:5]<- AUC_fin_table[,1:4]
Fin_all[3,1]<- "mean WAUC"
Fin_all[3,2:5]<- WAUC_fin_table[,1:4]
Fin_all[4,1]<- "AUC in"
Fin_all[4,2:5]<- AUC_fin_in_table[,1:4]
Fin_all[5,1]<- "AR dist VI loss 2"
Fin_all[5,2:5]<- VI_fin_table_ARdist2[,1:4]
Fin_all[6,1]<- "AR dist VI loss"
Fin_all[6,2:5]<- VI_fin_table_ARdist[,1:4]
Fin_all[7,1]<- "mean K"
Fin_all[7,2:5]<- c(mean(apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjam$modelList$burnin:fit_gjam$modelList$ng]),
                   mean(apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamDP1$modelList$burnin:fit_gjamDP1$modelList$ng]),
                   mean(apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamPY1$modelList$burnin:fit_gjamPY1$modelList$ng]),
                   mean(apply(fit_gjamPY2$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamPY2$modelList$burnin:fit_gjamPY2$modelList$ng]))
Fin_all[8,1]<- "K VI2"
Fin_all[8,2:5]<- VI_fin_table_K2[,1:4]
Fin_all[9,1]<- "K VI"
Fin_all[9,2:5]<- VI_fin_table_K[,1:4]
Fin_all[,2:5]<- round(Fin_all[,2:5], 3)

