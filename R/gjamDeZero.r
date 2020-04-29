
gjamDeZero <- function( ymat ){
  
  S <- ncol(ymat)
  n <- nrow(ymat)
  ynames <- colnames(ymat)
  if(is.null(ynames))ynames <- paste('S',1:S,sep='_')
  
  index <- which(ymat > 0)
  
  list(yvec = ymat[index], n = n, S = S, index = index, 
       ynames = ynames) 
}
