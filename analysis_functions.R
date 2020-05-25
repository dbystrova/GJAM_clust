

model_prediction_summary<- function(model, list_out_s,list_in_s , Ytest=Ydata_test, Ytrain=Ydata_train, pw=p_w){
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
  
  
  return(list(AUC_out=AUC_ous,AUC_in=AUC_GJAMDin, WAUC= AUC_ous*pw))
}


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
