

gjamPoints2Grid <- function(specs, xy, nxy = NULL, dxy = NULL, 
                            predGrid = NULL, effortOnly = TRUE){
  
  dx <- dy <- nx <- ny <- NULL
  
  wna <- which(is.na(xy),arr.ind=T)
  if(length(wna) > 0){
    wna <- unique(wna[,1])
  }
  wna <- c(wna,which(is.na(specs)))
  if(length(wna) > 0){
    specs <- specs[-wna]
    xy    <- xy[-wna,]
  }
  
  if(length(nxy) == 1) nx <- ny <- nxy
  if(length(dxy) == 1) dx <- dy <- dxy
  if(is.null(nxy) & is.null(dxy) & is.null(predGrid)){
    stop('must supply nxy or dxy or predGrid')
  }
  if(length(nxy) == 2) nx <- nxy[1]; ny <- nxy[2]
  if(length(dxy) == 2) dx <- dxy[1]; dy <- dxy[2]
  
  if(length(specs) != nrow(xy))stop('specs must have length = nrow(xy)')
  
  mr <- apply( apply(xy,2,range, na.rm=T), 2, diff )
  
  if(!is.null(dxy)){
    if(mr[1] < (3*dx) | mr[2] < (3*dy))stop('dxy too small')
  }
  .incidence2Grid(specs, lonLat = xy, nx, ny, dx, dy, predGrid, effortOnly)
}
