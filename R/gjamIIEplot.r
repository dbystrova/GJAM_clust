
gjamIIEplot <- function(fit, response, effectMu, effectSd = NULL, 
                        ylim = NULL, col='black', legLoc = 'topleft',
                        cex=1){
  
  mainEffect <- dirEffect <- intEffect <- indEffectTo <- NULL 
  mainSd <- dirSd <- intSd <- indSdTo <- NULL
  
  for(k in 1:length(fit$IIE))assign( names(fit$IIE)[[k]], fit$IIE[[k]] )
  
  factorList <- fit$modelSummary$factorList
  cnames     <- colnames(mainEffect)
  
  if(length(factorList) > 0){
    for(j in 1:length(factorList)){
      nk <- length(factorList[[j]])
      wk <- match(factorList[[j]],cnames)
      sk <- matrix( unlist(strsplit(cnames[wk],names(factorList)[[j]])),
                    ncol=2,byrow=T)[,2]
      w0 <- which(nchar(sk[1]) == 0)
      if(length(w0) > 0)sk[w0] <- paste(names(factorList)[[j]],c(1:length(w0)), 
                                        sep='_')
      if(length(w0) == 1)sk[w0] <- names(factorList)[[j]]
      cnames[wk] <- sk
    }
  }
  colnames(mainEffect) <- colnames(dirEffect) <- colnames(intEffect) <- 
    colnames(indEffectTo) <- cnames
      
  np <- length(effectMu)
  ns <- length(effectSd)
  
  if(length(col) == 1 & np > 1){
    colF <- colorRampPalette(c('blue','brown','orange'))
    col  <- colF(4)
  }
  
  stackList <- vector(np, mode = 'list')
  names(stackList) <- effectMu
  stackSd <- vector(ns, mode = 'list')
  names(stackSd) <- effectSd
  
  for(j in 1:np){
    if(effectMu[j] == 'main')   stackList[[j]] <- mainEffect[response,] 
    if(effectMu[j] == 'direct') stackList[[j]] <- dirEffect[response,] 
    if(effectMu[j] == 'int')    stackList[[j]] <- intEffect[response,] 
    if(effectMu[j] == 'ind')    stackList[[j]] <- indEffectTo[response,] 
  }
  for(j in 1:ns){
    if(effectSd[j] == 'main')   stackSd[[j]] <- mainSd[response,] 
    if(effectSd[j] == 'direct') stackSd[[j]] <- dirSd[response,] 
    if(effectSd[j] == 'int')    stackSd[[j]] <- intSd[response,] 
    if(effectSd[j] == 'ind')    stackSd[[j]] <- indSdTo[response,] 
  }
  
  .stackedBoxPlot( stackList = stackList, stackSd = stackSd, ylim = ylim,
                              col=col, barnames=names(stackList[[1]]),
                              decreasing=F, cex=cex)
  abline(h=0,col='grey',lwd=2)
  
  leg <- effectMu
  w1  <- grep('indTo',effectMu)
  if(length(w1) > 0)leg[w1] <- paste('indirect on',response)
  
  legend(legLoc,leg, text.col=col, bty='n')
}
  
  