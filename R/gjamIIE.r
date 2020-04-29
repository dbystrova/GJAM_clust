

gjamIIE <- function(output, xvector, MEAN = T, keepNames = NULL,
                    omitY = NULL, sdScaleX = T, sdScaleY = F){
  
  xMu    <- colMeans(output$inputs$x)
  xSd    <- apply(output$inputs$x,2,sd)
  standX <- cbind(xMu,xSd)
  colnames(standX) <- c('xmean','xsd')
  
  
  xii <- which(!names(xvector) %in% colnames(output$inputs$x))
  if(length(xii) > 0)xvector <- xvector[-xii]
  
  xii <- which(!colnames(output$inputs$x) %in% names(xvector))
  if(length(xii) > 0){
    stop('xvector is missing variables in model')
  }
  
  factorList <- output$inputs$factorBeta$factorList
  otherpar   <- output$modelList$reductList$otherpar
  
  IIE <- .directIndirectCoeffs( snames = colnames(output$inputs$y), xvector, 
                                chains = output$chains, MEAN,
                                factorList = factorList,
                                keepNames, omitY, sdScaleY, sdScaleX,
                                standX = standX, 
                                otherpar = otherpar,
                                REDUCT = output$modelList$REDUCT, 
                                ng = output$modelList$ng, 
                                burnin = output$modelList$burnin)
  fit <- append(output, list('IIE' = IIE))
  fit
}
