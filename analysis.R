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
library(cowplot)
library(rootSolve)
library(FactoMineR)
library(ggsci)
library(viridis)
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
#K_prior=16
#r_reduct = 5

folderpath="PCA_analysis/r5_25/"

##conditional prediction
columns<-1:ncol(Ydata_train)
ycs<- sample(columns, 10)
y_c_p <-columns[ !columns %in% ycs]


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





## Load models
fit_gjam<- load_object(paste0(folderpath,"fit_gjam.Rdata"))
fit_gjamDP2<- load_object(paste0(folderpath,"fit_gjamDP2.Rdata"))
fit_gjamPY1<- load_object(paste0(folderpath,"fit_gjamPY1.Rdata"))
fit_gjamPY2<- load_object(paste0(folderpath,"fit_gjamPY2.Rdata"))

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
gjamDP2<- model_prediction_summary(model=fit_gjamDP2, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamPY1<- model_prediction_summary(model=fit_gjamPY1, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)
gjamPY2<- model_prediction_summary(model=fit_gjamPY2, list_out_s=new_out,list_in_s=new_in , list_cond= new_cond, Ytest=Ydata_test, Ytrain=Ydata_train)

####Prediction for all models 
AUC_data<- tibble(GJAM=gjamDP$AUC_out, GJAM2=gjamDP2$AUC_out, PY1= gjamPY1$AUC_out, PY2=gjamPY2$AUC_out )
AUC_data_in<-  tibble(GJAM=gjamDP$AUC_in, GJAM2=gjamDP2$AUC_in, PY1= gjamPY1$AUC_in, PY2=gjamPY2$AUC_in )
AUC_data_cond<- tibble(GJAM=gjamDP$AUC_cond, GJAM2=gjamDP2$AUC_cond, PY1= gjamPY1$AUC_cond, PY2=gjamPY2$AUC_cond )
AUC_data$species<- colnames(Ydata_test)[1:ncol(Ydata_test)]
AUC_fin<- melt(AUC_data)

p2<-ggplot(data=AUC_fin)+geom_boxplot(aes(y=as.numeric(value),x=as.factor(variable),fill=as.factor(variable)))+
  scale_y_continuous(name="AUC")+
  scale_fill_discrete(name = "Models", labels = c("GJAM","GJAM2","PY1","PY2"))+xlab("Models")+ theme_bw() 
p2

AUC_fin_cond_table<- as.data.frame(t(apply(AUC_data_cond,2,mean)))
AUC_fin_in_table<- as.data.frame(t(apply(AUC_data_in,2,mean)))
AUC_fin_table<- as.data.frame(t(apply(AUC_data[,1:4],2,mean)))
WAUC_fin_table<- as.data.frame(t(apply(AUC_data[,1:4],2,function(x) sum(x*p_w))))
names(AUC_fin_table)<- c("GJAM","GJAM2","PY1","PY2")
kable(cbind(data.frame( Measure = rbind("AUC")), rbind(AUC_fin_table)), format="pandoc", caption= "Prediction")

################################################################################################################

#Convergence of Sigma
#Sigma
df1<- tibble(ES= effectiveSize(mcmc(fit_gjam$chains$sgibbs[fit_gjam$modelList$burnin:fit_gjam$modelList$ng,])))
df2<- tibble( ES=effectiveSize(mcmc(fit_gjamDP2$chains$sgibbs[fit_gjam$modelList$burnin:fit_gjam$modelList$ng,]) ))
df3<- tibble(ES=effectiveSize(mcmc(fit_gjamPY1$chains$sgibbs[fit_gjam$modelList$burnin:fit_gjam$modelList$ng,])))
df4<- tibble(ES =effectiveSize(mcmc(fit_gjamPY2$chains$sgibbs[fit_gjam$modelList$burnin:fit_gjam$modelList$ng,])))

#pdf(file = "Plots/Effective_size_sigma.pdf", width= 8.27, height = 9.69)
rbind(df1 %>% mutate(var = "DP"),
      df2 %>%  mutate(var = "DP2"), 
      df3 %>% mutate(var = "PY1")
      #df4 %>%  mutate(var = "PY2")
      ) %>% 
  ggplot(aes(ES, color = var, fill = var, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 
#dev.off()

#pdf(file = "Plots/Effective_size_beta.pdf", width= 8.27, height = 9.69)

df1<- tibble(ES= effectiveSize(mcmc(fit_gjam$chains$bgibbs)))
df2<- tibble( ES=effectiveSize(mcmc(fit_gjamDP2$chains$bgibbs) ))
df3<- tibble(ES=effectiveSize(mcmc(fit_gjamPY1$chains$bgibbs)))
df4<- tibble(ES =effectiveSize(mcmc(fit_gjamPY2$chains$bgibbs)))


rbind(df1 %>% mutate(var = "DP"),
      df2 %>%  mutate(var = "DP2"), 
      df3 %>% mutate(var = "PY1"),
      df4 %>%  mutate(var = "PY2")) %>% 
  ggplot(aes(ES, color = var, fill = var, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 
#dev.off()

############################################################################################################

### Prior/Posterior and convergence for hyperparameters 
p_DP2 =tibble(it= 1: length(fit_gjamDP2$chains$alpha.DP_g),
       DP2= fit_gjamDP2$chains$alpha.DP_g) %>%
  ggplot(aes(x=it,y=DP2))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title=TeX(sprintf('Trace plot for DP2 parameter $\\alpha$')))+xlab("iterations")+ylab("Concentration parameter DP2") +theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15)) +theme(legend.text=element_text(size=15))




### Prior/Posterior and convergence for hyperparameters 
p_PY2 =tibble(it= 1: length(fit_gjamPY2$chains$alpha.PY_g),
       PY2= fit_gjamPY2$chains$alpha.PY_g) %>%
  ggplot(aes(x=it,y=PY2))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title=TeX(sprintf('Trace plot for PY1 parameter $\\alpha$')))+xlab("iterations")+ylab("Concentration parameter PY2") +theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15)) +theme(legend.text=element_text(size=15))



### Prior/Posterior and convergence for hyperparameters 
p_PY2_d =tibble(it= 1: length(fit_gjamPY2$chains$discount.PY_g),
       PY2= fit_gjamPY2$chains$discount.PY_g) %>%
  ggplot(aes(x=it,y=PY2))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  labs(title=TeX(sprintf('Trace plot for parameter $\\sigma$')))+xlab("iterations")+ylab("Concentration parameter PY2") +theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15)) +theme(legend.text=element_text(size=15))

#pdf(file = "Plots/Posterior_for_PY_parameters.pdf", width= 8.27, height = 9.69)
prow <- plot_grid(
  p_DP2 + theme(legend.position='none'),
  p_PY2 + theme(legend.position='none'),
  p_PY2_d + theme(legend.position='none'),
  nrow = 2
)
#legend_b <- get_legend(p_DP2+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))
plot(p)
#focdev.off()


#df1<- tibble(ES= effectiveSize(mcmc(fit_gjamDP2$chains$alpha.DP_g)))
#df2<- tibble( ES=effectiveSize(mcmc(fit_gjamPY2$chains$alpha.PY_g) ))
#df3<- tibble(ES=effectiveSize(mcmc(fit_gjamPY2$chains$discount.PY_g)))
#rbind(df1 %>% mutate(var = "DP2"),
#      df2 %>%  mutate(var = "PY2_a"), 
#      df3 %>% mutate(var = "PY2_d")) %>% 
#  ggplot(aes(ES, color = var, fill = var, alpha = 0.3))+ geom_histogram( position="identity", alpha=0.2) 



alpha<-mcmc(fit_gjamDP2$chains$alpha.DP_g)
plot(alpha,main="alpha DP")
acfplot(alpha)
cumuplot(alpha)

alpha<-mcmc(fit_gjamPY2$chains$discount.PY_g)
plot(alpha,main="alpha DP")
acfplot(alpha)
cumuplot(alpha)

#### Add prior/posterior plot
prior_nu1 <- fit_gjamDP2$modelList$reductList$otherpar$shape
prior_nu2 <- fit_gjamDP2$modelList$reductList$otherpar$rate
alpha_vec<- rgamma(18000, prior_nu1,prior_nu2)
x<- sapply(alpha_vec, functionDPM,n=112,N=112)
mean(x)

###Posterior for alpha parameter DP2

aDP2 =tibble(Prior =alpha_vec,
       Posterior = fit_gjamDP2$chains$alpha.DP_g[(fit_gjamDP2$modelList$burnin+1):fit_gjamDP2$modelList$ng])%>%
  gather(Distribution, density, Prior:Posterior)%>%
 ggplot(aes(x=density, fill=Distribution, alpha=0.7)) +
  geom_density(adjust = 2)+ggtitle(TeX(sprintf('Prior and posterior distribution for DP2 $\\alpha$'))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))


prior_nu1 <- fit_gjamPY2$modelList$reductList$otherpar$shape
prior_nu2 <- fit_gjamPY2$modelList$reductList$otherpar$rate
alpha_vecPY<- rgamma(18000, prior_nu1,prior_nu2)
x<- sapply(alpha_vecPY, functionPY,n=112, sigma_py=0.25)
mean(x)

###Posterior for alpha parameter PY
aPY =tibble(Prior =alpha_vecPY,
       Posterior = fit_gjamPY2$chains$alpha.PY_g[(fit_gjamPY2$modelList$burnin+1):fit_gjamPY2$modelList$ng])%>%
  gather(Distribution, density, Prior:Posterior)%>%
  ggplot(aes(x=density, fill=Distribution)) +
  geom_density(adjust = 2, alpha=0.5)+ggtitle(TeX(sprintf('Prior and posterior distribution for PY2 $\\alpha$'))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))



##Posterior for discount parameter
sPY =tibble(Prior =fit_gjamPY2$chains$discount.PY_g[(fit_gjamPY2$modelList$burnin+1):fit_gjamPY2$modelList$ng],
       Posterior = fit_gjamPY2$chains$discount.PY_g[(fit_gjamPY2$modelList$burnin+1):fit_gjamPY2$modelList$ng])%>%
  gather(Distribution, density, Prior:Posterior)%>%
  ggplot(aes(x=density, fill=Distribution))  +
  geom_density(adjust = 2,alpha=0.5)+ggtitle(TeX(sprintf('Prior and posterior distribution for PY2 $\\sigma$')))  +
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))


#pdf(file = "Plots/Posterior_for_PY_parameters2.pdf", width= 8.27, height = 9.69)
prow <- plot_grid(
  aDP2 + theme(legend.position='none'),
  aPY + theme(legend.position='none'),
  sPY + theme(legend.position='none'),
  nrow = 2
)
legend_b <- get_legend(aDP2+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(9, 1))
plot(p)
#dev.off()


##### Posterior distribution for the number of clusters

x =apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))
hist(x)


############################Trace plot##########################################################################
#pdf(file = "Plots/Trace_plot_partitions.pdf", width= 8.27, height = 9.69)

tibble(it= 1: length(apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))),
              DP= apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x))),
              DP2 =apply(fit_gjamDP2$chains$kgibbs,1,function(x) length(unique(x))),
              PY1=apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x))),
              PY2=apply(fit_gjamPY2$chains$kgibbs,1,function(x) length(unique(x))) ) %>%
gather(Model, trace, DP:PY2)%>%
 ggplot(aes(x=it,y=trace,col=Model))+geom_line(alpha=0.8)+ scale_color_viridis(discrete=TRUE)+
 labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
 theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
 theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
       axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
       plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))

#dev.off()

##### Truncation weights ############################################################

last_pk<- round(mean(fit_gjamDP2$chains$pk_g[-c(1:fit_gjamDP2$modelList$burnin),fit_gjamDP2$modelList$reductList$N]),5)
df_weights <- tibble(pw =apply(fit_gjamDP2$chains$pk_g[-c(1:fit_gjamDP2$modelList$burnin),],2,mean),
                     tr= 1:ncol(fit_gjamDP2$chains$pk_g))
pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
  geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
  geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Mean weights for DP2, p_last= ",last_pk), 
                                                                                                    caption=paste0("Number of iterations: ",fit_gjamDP2$modelList$ng," burnin: ",fit_gjamDP2$modelList$burnin))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
pl_weigths


last_pk<- round(mean(fit_gjamPY1$chains$pk_g[-c(1:fit_gjamPY1$modelList$burnin),fit_gjamPY1$modelList$reductList$N]),5)
df_weights <- tibble(pw =apply(fit_gjamPY1$chains$pk_g[-c(1:fit_gjamPY2$modelList$burnin),],2,mean),
                     tr= 1:ncol(fit_gjamPY1$chains$pk_g))
pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
  geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
  geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Mean weights for PY1, p_last= ",last_pk), 
                                                                                                    caption=paste0("Number of iterations: ",fit_gjamPY2$modelList$ng," burnin: ",fit_gjamPY2$modelList$burnin))+
  theme_bw() + theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))
pl_weigths



last_pk<- round(mean(fit_gjamPY2$chains$pk_g[-c(1:fit_gjamPY2$modelList$burnin),fit_gjamPY2$modelList$reductList$N]),5)
df_weights <- tibble(pw =apply(fit_gjamPY2$chains$pk_g[-c(1:fit_gjamPY2$modelList$burnin),],2,mean),
                     tr= 1:ncol(fit_gjamPY2$chains$pk_g))
pl_weigths<- ggplot(df_weights, aes(x=tr, y=pw)) +
  geom_segment( aes(x=tr,xend=tr,y=0,yend=pw)) +
  geom_point( size=0.5, color="red", fill=alpha("blue", 0.3), alpha=0.4, shape=21, stroke=2)+  labs(title=paste0("Mean weights for PY2, p_last= ",last_pk), 
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



gjamClust<- function(model, true_cl = True_clust$K_n){
  if("other" %in% colnames(model$inputs$y)){
    sp_num<- ncol(model$inputs$y)-1
  } 
  MatDP<- relabel_clust(model$chains$kgibbs[(model$modelList$burnin+1):model$modelList$ng,1:sp_num])
  tr_cl<- true_cl
  CM_DP<-  comp.psm(MatDP)
  mbind_DP.ext <- minbinder.ext(CM_DP,MatDP, method="all",include.greedy = TRUE)
  vi_DP.ext <- minVI(CM_DP,MatDP, method="all", include.greedy = TRUE)
  return(list(Bind_est=mbind_DP.ext$cl[1,], VI_est = vi_DP.ext$cl[1,]))
}

DP_clust<- gjamClust(model= fit_gjam)
DP2_clust<- gjamClust(model= fit_gjamDP2)
PY1_clust<- gjamClust(model= fit_gjamPY1)
PY2_clust<- gjamClust(model= fit_gjamPY2)


#### Compare the obtained estimates with the PFG clusters


BI_fin_table_K<- tibble(GJAM = length(unique(DP_clust$Bind_est)),GJAM2=length(unique(DP2_clust$Bind_est)),PY1=length(unique(PY1_clust$Bind_est)),PY2=length(unique(PY2_clust$Bind_est)))
BI_fin_table_VIdist<- tibble(GJAM = vi.dist(DP_clust$Bind_est, True_clust$K_n),GJAM2 = vi.dist(DP2_clust$Bind_est, True_clust$K_n),PY1= vi.dist(PY1_clust$Bind_est, True_clust$K_n),PY2 = vi.dist(PY2_clust$Bind_est, True_clust$K_n))
BI_fin_table_ARdist<- tibble(GJAM = arandi(DP_clust$Bind_est, True_clust$K_n),GJAM2 = arandi(DP2_clust$Bind_est, True_clust$K_n),PY1 = arandi(PY1_clust$Bind_est, True_clust$K_n),PY2 = arandi(PY2_clust$Bind_est, True_clust$K_n))


VI_fin_table_K<-tibble(GJAM = length(unique(DP_clust$VI_est)),GJAM2=length(unique(DP2_clust$VI_est)),PY1=length(unique(PY1_clust$VI_est)),PY2=length(unique(PY2_clust$VI_est)))
VI_fin_table_VIdist<-  tibble(GJAM = vi.dist(DP_clust$VI_est, True_clust$K_n),GJAM2 = vi.dist(DP2_clust$VI_est, True_clust$K_n),PY1= vi.dist(PY1_clust$VI_est, True_clust$K_n),PY2 = vi.dist(PY2_clust$VI_est, True_clust$K_n))
VI_fin_table_ARdist<- tibble(GJAM = arandi(DP_clust$VI_est, True_clust$K_n),GJAM2 = arandi(DP2_clust$VI_est, True_clust$K_n),PY1 = arandi(PY1_clust$VI_est, True_clust$K_n),PY2 = arandi(PY2_clust$VI_est, True_clust$K_n))





### Old style###############################################################################################################
sp_num<- ncol(Ydata)-1
MatDP<- relabel_clust(fit_gjam$chains$kgibbs[(fit_gjam$modelList$burnin+1):fit_gjam$modelList$ng,1:sp_num])
MatDP2<- relabel_clust(fit_gjamDP2$chains$kgibbs[(fit_gjamDP2$modelList$burnin+1):fit_gjamDP2$modelList$ng,1:sp_num])
MatPY1<- relabel_clust(fit_gjamPY1$chains$kgibbs[(fit_gjamPY1$modelList$burnin+1):fit_gjamPY1$modelList$ng,1:sp_num])
MatPY2<- relabel_clust(fit_gjamPY2$chains$kgibbs[(fit_gjamPY2$modelList$burnin+1):fit_gjamPY2$modelList$ng,1:sp_num])

### Obtain Posterior similarity matrices
tr_cl<-True_clust$K_n
CM_DP<-  comp.psm(MatDP)
CM_DP2<- comp.psm(MatDP2)
CM_PY1<- comp.psm(MatPY1)
CM_PY2<- comp.psm(MatPY2)

mbind_DP.ext <- minbinder.ext(CM_DP,MatDP, method="all",include.greedy = TRUE)
mbind_DP2.ext <- minbinder.ext(CM_DP2,MatDP2, method="all",include.greedy = TRUE)
mbind_PY1.ext <- minbinder.ext(CM_PY1,MatPY1, method="all",include.greedy = TRUE)
mbind_PY2.ext <- minbinder.ext(CM_PY2,MatPY2, method="all",include.greedy = TRUE)


G_BI_cl_fin_table<- tibble(GJAM = t(c(length(unique(mbind_DP.ext$cl[1,])),length(unique(mbind_DP2.ext$cl[1,])),length(unique(mbind_PY1.ext$cl[1,])),length(unique(mbind_PY2.ext$cl[1,])))))
names(G_BI_cl_fin_table)<- c("GJAM","GJAM2","PY1","PY2")
formattable(G_BI_cl_fin_table)


Gbi.dist_DP<- vi.dist(mbind_DP.ext$cl[1,], tr_cl)
Gbi.dist_DP2<- vi.dist(mbind_DP2.ext$cl[1,], tr_cl)
Gbi.dist_PY1<- vi.dist(mbind_PY1.ext$cl[1,], tr_cl)
Gbi.dist_PY2<- vi.dist(mbind_PY2.ext$cl[1,], tr_cl)
BI_VIloss_fin_table<- as.data.frame(t(c(Gbi.dist_DP,Gbi.dist_DP2,Gbi.dist_PY1,Gbi.dist_PY2)))
formattable(BI_VIloss_fin_table, "VI distance, binder loss")


# compare clusterings found by different methods with true grouping for Binder loss
Ar_D_DP<- arandi(mbind_DP.ext$cl[1,], tr_cl)
#Ar_D_DP1<- arandi(mbind_DP1.ext$cl[1,], tr_cl)
Ar_D_DP2<-arandi(mbind_DP2.ext$cl[1,], tr_cl)
Ar_D_PY1<-arandi(mbind_PY1.ext$cl[1,], tr_cl)
Ar_D_PY2<-arandi(mbind_PY2.ext$cl[1,], tr_cl)
Ar_D_fin_table<- as.data.frame(t(c(Ar_D_DP,Ar_D_DP2,Ar_D_PY1,Ar_D_PY2)))
names(Ar_D_fin_table)<- c("GJAM","GJAM2","PY1","PY2")
print(Ar_D_fin_table)
formattable(Ar_D_fin_table, "R index distance, binder loss")



################################### GREEEDY VI loss 
G_vi_DP <- minVI(CM_DP,MatDP, method="all", include.greedy = TRUE)
G_vi_DP2 <- minVI(CM_DP2,MatDP2, method="all", include.greedy = TRUE)
G_vi_PY1 <- minVI(CM_PY1,MatPY1, method="all", include.greedy = TRUE, l=300)
G_vi_PY2 <- minVI(CM_PY2,MatPY2, method="all", include.greedy = TRUE)

G_VI_cl_fin_table<- as.data.frame(t(c(length(unique(G_vi_DP$cl[1,])),length(unique(G_vi_DP1$cl[1,])),length(unique(G_vi_DP2$cl[1,])),length(unique(G_vi_PY1$cl[1,])),length(unique(G_vi_PY2$cl[1,])))))
names(G_VI_cl_fin_table)<- c("GJAM","GJAM1","GJAM2","PY1","PY2")
formattable(G_VI_cl_fin_table)


Ar_VI_DP<- arandi(G_vi_DP$cl[1,], tr_cl)
Ar_VI_DP1<- arandi(G_vi_DP1$cl[1,], tr_cl)
Ar_VI_DP2<-arandi(G_vi_DP2$cl[1,], tr_cl)
Ar_VI_PY1<-arandi(G_vi_PY1$cl[1,], tr_cl)
Ar_VI_PY2<-arandi(G_vi_PY2$cl[1,], tr_cl)
Ar_VI_fin_table<- as.data.frame(t(c(Ar_VI_DP,Ar_VI_DP1,Ar_VI_DP2,Ar_VI_PY1,Ar_VI_PY2)))
names(Ar_VI_fin_table)<- c("GJAM","GJAM1","GJAM2","PY1","PY2")
formattable(Ar_VI_fin_table,"R index distances, VI loss")


Gvi.dist_DP<- vi.dist(G_vi_DP$cl[1,], tr_cl)
Gvi.dist_DP1<- vi.dist(G_vi_DP1$cl[1,], tr_cl)
Gvi.dist_DP2<- vi.dist(G_vi_DP2$cl[1,], tr_cl)
Gvi.dist_PY1<- vi.dist(G_vi_PY1$cl[1,], tr_cl)
Gvi.dist_PY2<- vi.dist(G_vi_PY2$cl[1,], tr_cl)
VI_VIloss_fin_table<- as.data.frame(t(c(Gvi.dist_DP,Gvi.dist_DP1,Gvi.dist_DP2,Gvi.dist_PY1,Gvi.dist_PY2)))
formattable(VI_VIloss_fin_table,"VI distance, VI loss")
###############################################################################################################
############ greedy EPL ##### another algorithm 
DP_grEPL<- MinimiseEPL(MatDP, pars = list())
length(unique(DP_grEPL$decision))
arandi(DP_grEPL$decision,G_vi_DP$cl[1,] )
DP2_grEPL<- MinimiseEPL(MatDP2, pars = list(loss_type="NVI"))
length(unique(DP2_grEPL$decision))
arandi(DP2_grEPL$decision,G_vi_DP2$cl[1,] )
PY1_grEPL<- MinimiseEPL(MatPY1, pars = list())
length(unique(PY1_grEPL$decision))
PY2_grEPL<- MinimiseEPL(MatPY2, pars = list(loss_type="NVI"))
length(unique(PY2_grEPL$decision))
##### confidence intervals  ############################################################

DP_cb = credibleball(DP_clust$VI_est, MatDP, c.dist = c("VI","Binder"), alpha = 0.05)

############################################################################################################

#summary(DP_cb)
#The credible ball characterizes the uncertainty in the clustering esitmate.
#It can be summarized with:
#1. upper vertical bound: partitions in the ball with the fewest clusters that are most distant,
#2. lower vertical bound: partitions in the ball with the most clusters that are most distant,
#3. horizontal bound: partitions in the ball with the greatest distance.
#The upper vertical bound has 89 clusters with a distance of 1.21.
#The lower vertical bound has 111 clusters with a distance of 1.88.
#The horizontal bound has 110 clusters with a distance of 1.88.

DP_cb = credibleball(DP_clust$VI_est, MatDP2, c.dist = c("VI"), alpha = 0.05)
DP2_cb = credibleball(DP2_clust$VI_est, MatDP2, c.dist = c("VI"), alpha = 0.05)
PY1_cb = credibleball(PY1_clust$VI_est, MatPY1, c.dist = c("VI"), alpha = 0.05)
PY2_cb = credibleball(PY2_clust$VI_est, MatPY2, c.dist = c("VI"), alpha = 0.05)


#### Cluster names
Cluster_models<- tibble( CODE_CBNA=colnames(Ydata)[1:(ncol(Ydata)-1)],  ClustDP=DP_clust$VI_est,
                         ClustDP2=DP2_clust$VI_est,ClustPY1=PY1_clust$VI_est,ClustPY2=PY2_clust$VI_est, PFG=True_clust$K_n)

save(Cluster_models, file =  paste0(folderpath,"Clusters_modells.Rdata"))
##################################Covariance matrix######################################################################


#### Add cluster labels 


Colnames_Y_clust<- merge(Colnames_Y, Cluster_models, by ="CODE_CBNA")

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
MDP2= fit_gjamDP2$parameters$corMu
Colnames_Y_clust$Sp_name_DP2 <- paste(Colnames_Y_clust$ClustDP2, Colnames_Y_clust$species,sep="_")
rownames(MDP2)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP2"],"other")
colnames(MDP2)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_DP2"],"other")
colors_vir=viridis(length(unique(Colnames_Y_clust$ClustDP2))+1, option = "magma")
LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_DP2),"ClustDP2"],15), function(x) colors_vir[x])

corrplot(MDP2, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.45,tl.srt=45, method = "color",col=cols(200),
         type = "lower", title= "Correlation for the DP2 model (original)", mar=c(0,0,1,0))

corrplot(MDP2, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45, tl.col=LabelCol,method = "color",col=cols(200),
         type = "full", title= "Correlation for the DP2 model (DP2 groups)", mar=c(0,0,1,0))

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
MPY2= fit_gjamPY2$parameters$corMu
Colnames_Y_clust$Sp_name_PY2 <- paste(Colnames_Y_clust$ClustPY1, Colnames_Y_clust$species,sep="_")
rownames(MPY2)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_PY2"],"other")
colnames(MPY2)=c(Colnames_Y_clust[order(Colnames_Y_clust$CN),"Sp_name_PY2"],"other")
colors_vir=viridis(length(unique(Colnames_Y_clust$ClustPY2))+1, option = "magma")
LabelCol = sapply(c(Colnames_Y_clust[order(Colnames_Y_clust$Sp_name_PY2),"ClustPY2"],length(unique(Colnames_Y_clust$ClustPY2))+1), function(x) colors_vir[x])

corrplot(MPY1, diag = FALSE, order = "hclust", tl.cex = 0.45,tl.srt=45, method = "color",col=cols(200),
         type = "full", title= "Correlation for the PY2 model (original)", mar=c(0,0,1,0))

corrplot(MPY1, diag = FALSE, order = "alphabet", tl.cex = 0.45,tl.srt=45, tl.col=LabelCol,method = "color",col=cols(200),
         type = "full", title= "Correlation for the PY2 model (PY2 groups)", mar=c(0,0,1,0))

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

A= cov_matrix(fit=fit_gjamDP2, burn_period=fit_gjamDP2$modelList$burnin + 8000,iterations=fit_gjamDP2$modelList$ng)


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
Fin_all<-as.data.frame(matrix(NA,nrow=12,ncol=9))
names(Fin_all)<- c("Parameter","GJAM","GJAM2","PY1","PY2","r", "iter", "burn","formula")
Fin_all$iter<- fit_gjam$modelList$ng
Fin_all$burn<- fit_gjam$modelList$burnin
Fin_all$r<-fit_gjam$modelList$reductList$r
Fin_all$formula<-as.character(form)
Fin_all[1,1]<- "DIC"
Fin_all[1,2:5]<- c(fit_gjam$fit$DIC,fit_gjamDP2$fit$DIC,fit_gjamPY1$fit$DIC,fit_gjamPY2$fit$DIC)/10000
Fin_all[2,1]<- "mean AUC"
Fin_all[2,2:5]<- AUC_fin_table
Fin_all[3,1]<- "mean WAUC"
Fin_all[3,2:5]<- WAUC_fin_table
Fin_all[4,1]<- "AUC in"
Fin_all[4,2:5]<- AUC_fin_in_table
Fin_all[5,1]<- "AUC cond"
Fin_all[5,2:5]<- AUC_fin_cond_table
Fin_all[6,1]<- "VI dist Binder loss"
Fin_all[6,2:5]<- BI_fin_table_VIdist
Fin_all[7,1]<- "AR dist Binder loss"
Fin_all[7,2:5]<- BI_fin_table_ARdist
Fin_all[8,1]<- "VI dist VI loss"
Fin_all[8,2:5]<- VI_fin_table_VIdist
Fin_all[9,1]<- "AR dist VI loss"
Fin_all[9,2:5]<- VI_fin_table_ARdist
Fin_all[10,1]<- "mean K"
Fin_all[10,2:5]<- c(mean(apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjam$modelList$burnin:fit_gjam$modelList$ng]),
                    mean(apply(fit_gjamDP2$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamDP2$modelList$burnin:fit_gjamDP2$modelList$ng]),
                    mean(apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamPY1$modelList$burnin:fit_gjamPY1$modelList$ng]),
                    mean(apply(fit_gjamPY2$chains$kgibbs,1,function(x) length(unique(x)))[fit_gjamPY2$modelList$burnin:fit_gjamPY2$modelList$ng]))
Fin_all[11,1]<- "K Bind"
Fin_all[11,2:5]<- BI_fin_table_K
Fin_all[12,1]<- "K VI"
Fin_all[12,2:5]<- VI_fin_table_K
Fin_all[,2:5]<- round(Fin_all[,2:5], 3)
#save(Fin_all, file = paste0(folderpath,"Fin_tab_r5.Rdata"))
#library("xlsx")
# Write the first data set in a new workbook
#write.xlsx(Fin_all, file = "Final_table.xlsx")
#########################################################################################
## Probably re-do the conditional prediction






