
gjamTrimY <- function(y, minObs = 2, maxCols = NULL, OTHER = TRUE){  
    
    # minObs    - minimum no. of non-zero values in a column of y
    # maxCols   - number of columns to retain, those with highest values
    # OTHER     - logical or names to include in 'other' class
    # if(OTHER) sum of rare are returned in 'other' column
    # if already a column 'other', they are combined
    
    y      <- as.matrix(y)
    nc     <- ncol(y)
    ci     <- 1:nc
    mnames <- colnames(y)
    other  <- numeric(0)
    
    if(is.character(OTHER)){
      other <- y[,OTHER,drop=F]
      ci    <- ci[!mnames %in% OTHER]
      y     <- y[,ci,drop=F]
      OTHER <- T
    }
    
    io <- y
    io[io > 0] <- 1
    
    csum <- colSums(io, na.rm=T)
    ww   <- which(csum >= minObs)
    
    if(!is.null(maxCols)){
      ww <-  ww[ order(csum[ww],decreasing=T) ]
      ww <-  ww[1:maxCols]
    }
    ci     <- ci[ww]
    out    <- y[,ww]
    mnames <- mnames[ww]
    
    if(OTHER){
      other  <- rowSums(cbind(other,y[,-ww]),na.rm=T)
      out    <- cbind(out,other)
      mnames <- c(mnames,'other')
      ww <- which(colnames(out) == 'other')
      if(length(ww) > 1){
        other <- rowSums(out[,ww])
        out <- cbind(out[,-ww],other)
      }
    }
    
    if(!is.matrix(out)){
      out <- matrix(out,ncol=1)
      colnames(out) <- mnames
    }
    csum <- csum[ci]
    
    list(y = out, colIndex = ci, nobs = csum)
  }
  