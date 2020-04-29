
gjamReZero <- function( yDeZero ){
  
  ymat <- matrix(0, yDeZero$n, yDeZero$S)
  ymat[yDeZero$index] <- yDeZero$yvec
  colnames(ymat) <- yDeZero$ynames
  
  ymat
}
