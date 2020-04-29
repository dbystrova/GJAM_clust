
gjamPriorTemplate <- function(formula, xdata, ydata, lo = NULL, hi = NULL){
  
  # template for prior coefficient matrix
  # lo   - list of lower limits
  # hi   - list of upper limits
  
  tmp <- grep('_',colnames(ydata))
  if(length(tmp) > 0)stop("remove '_' from colnames(ydata)")
  tmp <- grep('_',colnames(xdata))
  if(length(tmp) > 0)stop("remove '_' from colnames(xdata)")
  
  x      <- model.matrix(formula,xdata)
  S      <- ncol(ydata)                    # no. responses
  Q      <- ncol(x)                        # no. predictors
  xnames <- colnames(x)
  xnames[1] <- 'intercept'
  ynames <- colnames(ydata)
  beta   <- matrix(0,Q,S)
  rownames(beta) <- xnames         
  colnames(beta) <- ynames
  
  loBeta <- beta - Inf
  hiBeta <- beta + Inf
  
  if(!is.null(lo)){
    loBeta <- .setLoHi(plist = lo, pmat = loBeta, xnames, ynames)
  }
  if(!is.null(hi)){
    hiBeta <- .setLoHi(plist = hi, pmat = hiBeta, xnames, ynames)
  }
  
  attr(loBeta,'formula') <- attr(hiBeta,'formula') <- formula
  
  list(lo = loBeta, hi = hiBeta)
}


  
