
gjamSimData <- function( n = 1000, S = 10, Q = 5, x = NULL, nmiss = 0, 
                         typeNames, effort = NULL ){
  
  # nmiss  = number of missing values in x
  # effort = used for DA data, format: 
  # list(columns = 1:S, values = rep(1,n))
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  
  wc <- which(typeNames == 'CC')
  if(length( .between(length(wc),1,2) ) > 0 )
    stop('CC must have at least 3 columns')
  wc <- which(typeNames == 'FC')
  if(length( .between(length(wc),1,2) ) > 0)
    stop('FC must have at least 3 columns')
  
  if( !is.null(x) ){
    n <- nrow(x)
    Q <- ncol(x)
  }
  
  .simData( n, S, Q, x, typeNames, nmiss, effort)
}
