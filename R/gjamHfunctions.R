.aggData <- function(cnames, data, gather, FUN){
  
  #cnames - column names in data to operate on
  #gather   - list of indices, all of same length
  
  cnames <- cnames[cnames %in% colnames(data)]
  if(length(cnames) == 0)return( list(data = numeric(0)) )
  
  FAC <- F
  df <- data[,cnames]
  
  if(FUN == 'factor'){
    FAC <- T
    tf <- fnames <- character(0)
    df <- numeric(0)
    
    for(j in 1:length(cnames)){
      kf <- as.character(data[,cnames[j]])
      tf <- cbind(tf, kf)
      fnames <- append(fnames,list(sort(unique(kf))))
      df <- cbind(df, match(kf,fnames[[j]]) )
    }
    colnames(df) <- cnames
    FUN <- 'max'
  }
  
  tmp <- aggregate(df, by=gather, FUN=FUN)
  ord <- do.call(order,list(tmp[,names(gather)]))
  colnames(tmp)[-c(1:length(gather))] <- cnames
  
  if(FAC){
    for(j in 1:length(cnames)){
      kf <- fnames[[j]][tmp[,cnames[j]]]
      tmp[,cnames[j]] <- as.factor(kf)
    }
  }
  list(ord = ord, data = tmp[ord,])
}

.setLoHi <- function(plist, pmat, xnames, ynames){
  
  # called by gjamPriorTemplate
  
  pnames <- names(plist)
  
  wx <- which(pnames %in% xnames)
  for(k in wx)pmat[pnames[k],] <- plist[[k]]
  
  wy <- which(pnames %in% ynames)
  for(k in wy)pmat[,pnames[k]] <- plist[[k]]
  
  comb <- grep('_',pnames)       # combination of x, y
  
  for(k in comb){
    ck <- unlist( strsplit(pnames[k],'_') )
    ik <- 1; jk <- 2
    wx <- which(xnames == ck[ik])
    if(length(wx) == 0){
      ik <- 2; jk <- 1
      wx <- which(xnames == ck[ik])
    }
    wy <- which(ynames == ck[jk])
    pmat[wx,wy] <- plist[[comb[k]]]
  }
  pmat
}

.myBoxPlot <- function(mat, tnam, snames, specColor, label, LEG=F){
  
  # tnam is columns of mat, with values of snames used to match specColor
  
  ord <- order(colMeans(mat),decreasing=F)
  mat  <- mat[,ord]
  tnam <- tnam[ord]
  bb   <- specColor[ match(tnam, snames) ]
  ry   <- range(mat)
  ymin <- min(mat) - diff(ry)*.15
  ymax <- max(mat) + diff(ry)*.15
  bx   <- .getColor(bb,.4)
  
  tmp <- .boxplotQuant( mat, xaxt='n',outline=F,ylim=c(ymin,ymax),
                        col=bx, border=bb, xaxt='n',lty=1)
  abline(h=0,lwd=2,col='grey',lty=2)
  
  dy <- .05*diff(par()$yaxp[1:2])
  
  cext <- .fitText2Fig(tnam,fraction=1)
  text((1:length(ord)) - .1,dy + tmp$stats[5,],tnam,srt=70,pos=4,
       col=bb, cex=cext)
  
  pl    <- par('usr')
  xtext <- pl[1]
  ytext <- pl[3] + diff(pl[3:4])*.85
  .plotLabel(label,location='topleft', cex=1.0)
  
  lg <- unique(names(bb))
  if(!is.null(lg) & LEG){
    colb <- bb[lg]
    legend('bottomright',legend=lg,text.col=colb,bty='n')
  }
  
}

.splitNames <- function(nameVec, snames=NULL, split='_'){
  
  vnam <- matrix( unlist( strsplit(nameVec,split) ),ncol=2,byrow=T)
  
  if(is.null(snames))return( list(vnam = vnam, xnam = NULL) )
  
  ix <- 1
  nc <- 2
  y1 <- which(vnam[,1] %in% snames)
  y2 <- which(vnam[,2] %in% snames)
  if(length(y1) > length(y2))ix <- 2
  if(ix == 2)nc <- 1
  
  xnam <- vnam[,ix]
  vnam <- vnam[,nc]
  list(vnam = vnam, xnam = xnam)
}

.updateBetaMet <- function(X, Y, B, lo, hi, loc, REDUCT, sig=NULL,
                           sinv=NULL, sp=.01 ){
  
  # metropolis update with TIME
  
  bnew <- B
  mnow <- X%*%B
  bnew[loc] <- .tnorm(nrow(loc),lo[loc],hi[loc],B[loc],sp)
  mnew <- X%*%bnew
  
  if(REDUCT){
    pnow <-  colSums( dnorm(Y,mnow,sqrt(sig),log=T) ) #cond ind species
    pnew <-  colSums( dnorm(Y,mnew,sqrt(sig),log=T) )
    z <- which( runif(length(pnow),0,1) < exp(pnew - pnow) )
    if(length(z) > 0)B[,z] <- bnew[,z]
  }else{
    pnow <- sum( .dMVN(Y,mnow,sinv=sinv,log=T) )
    pnew <- sum( .dMVN(Y,mnew,sinv=sinv,log=T) )
    z <- runif(1,0,1) < exp(pnew - pnow) 
    if(z)B <- bnew
  }
  B
}

.getPattern <- function(mat, wloc){
  
  mat <- mat*0 + 1
  
  rows <- numeric(0)
  pattern <- numeric(0)
  aa   <- mat*0
  aa[wloc] <- 1             # indicate keepers
  
  U  <- nrow(mat)
  SS <- ncol(mat)
  cc <- 1:U
  
  wk <- which( rowSums(abs(mat+1), na.rm=T) == 0 )  #none to sample
  if(length(wk) > 0){
    aa[cc[wk],] <- NA
    cc <- cc[-wk]
  }
  
  if(length(cc) == 0)return( list(rows = rows, pattern = pattern) )
  
  for(k in 1:U){
    
    if(length(cc) == 0)break
    
    ak <- aa[drop=F,cc,]
    ac <- matrix(ak[1,],nrow(ak),ncol(ak),byrow=T)
    ac[is.na(ac)] <- 0
    am <- ak - ac
    
    w1 <- which( duplicated(am) & rowSums(am,na.rm=T) == 0 )
    w1 <- c(1,w1)
    
    cw <- cc[w1]
    
    rr <- matrix( cw, 1 )
    pp <- matrix( which(aa[rr[1],] != 0), 1) #length-0 means all have no prior
    if(length(pp) == 0)pp <- matrix(1:SS,1)
    
    aa[cw,] <- NA
    cc <- cc[-w1]
    
    if(length(rows) == 0){
      rows <- rr
      pattern <- pp
      next
    } else {
      if(ncol(rr) == ncol(rows))rows <- rbind(rows,rr)
      if(ncol(rr) > ncol(rows)){
        rk <- matrix(NA,nrow(rows),ncol(rr))
        rk[,1:ncol(rows)] <- rows
        rows <- rbind(rk,rr)
      }
      if(ncol(rows) > ncol(rr)){
        rj <- matrix(NA,1,ncol(rows))
        rj[1:length(rr)] <- rr
        rows <- rbind(rows,rj)
      }
      
      if(ncol(pp) == ncol(pattern))pattern <- rbind(pattern,pp)
      if(ncol(pp) > ncol(pattern)){
        rk <- matrix(NA,nrow(pattern),ncol(pp))
        rk[,1:ncol(pattern)] <- pattern
        pattern <- rbind(rk,pp)
      }
      if(ncol(pattern) > ncol(pp)){
        rj <- matrix(NA,1,ncol(pattern))
        rj[1:length(pp)] <- pp
        pattern <- rbind(pattern,rj)
      }
    }
    
    if(length(cc) == 0)break
    
  }
  
  list(rows = rows, pattern = pattern)
}

.alphaPrior <- function(w, tindex, alphaPrior){
  
  S <- ncol(w)
  lo <- alphaPrior$lo
  hi <- alphaPrior$hi
  
  alpha <- (lo + hi)/2
  
  tmp    <- .getURowCol(alpha)
  uindex <- tmp$uindex
  Amat   <- tmp$Amat
  wA     <- tmp$wA
  aindex <- tmp$aindex
  
  loA <- hiA <- Amat
  loA[ wA ] <- lo[is.finite(alpha)]
  hiA[ wA ] <- hi[is.finite(alpha)]
  
  list(Amat = Amat, loAmat = loA, hiAmat = hiA, 
       wA = wA, uindex = uindex, aindex = aindex)
}

.betaPrior <- function(beta, notOther, loBeta, hiBeta){
  
  BPRIOR <- F
  loB <- hiB <- NULL
  
  bg   <- beta
  
  wB      <- which(!is.na(t(loBeta[,notOther])), arr.ind=T)[,c(2,1)]
  colnames(wB) <- c('row','col')
  
  bg <- (loBeta + hiBeta)/2
  
  loB <- loBeta[,notOther]
  hiB <- hiBeta[,notOther]
  list(beta = bg, loB = loB, hiB = hiB, wB = wB, BPRIOR = BPRIOR)
}

.lambdaPrior <- function(lprior, w, x, tindex, xnames, 
                         snames, other, notOther){
  
  loLambda <- lprior$lo
  hiLambda <- lprior$hi
  lambda   <- (loLambda + hiLambda)/2
  
  loLambda[,other] <- hiLambda[,other] <- NA
  
  lkeep    <- which(is.finite(loLambda))
  
  timeZero <- NULL
  
  M  <- nrow(lambda)
  rownames(lambda)[1] <- 'intercept'
  S  <- ncol(lambda)
  SS <- length(notOther)
  n  <- nrow(x)
  wz   <- w
  
  gindex <- kronecker(diag(S),rep(1,M)) 
  gindex <- gindex[lkeep,]
  
  wg     <- which(gindex == 1,arr.ind=T)
  wc     <- matrix(rep(1:M,S*M),S*M,S)[lkeep,]
  rowG   <- wc[wg]
  gindex <- cbind(rowG,wg)
  tmp    <- as.vector( t(outer(colnames(lambda)[notOther],
                               rownames(lambda),paste,sep='_') ) )
  rownames(gindex) <- tmp[lkeep]
  
  colX <- match(rownames(lambda),colnames(x))
  colX <- colX[rowG]
  gindex <- cbind(colX, gindex)
  colnames(gindex)[3:4] <- c('rowL','colW')
  nV <- nrow(gindex)
  
  Vmat <- matrix(0,n,nV)
  wz[wz < 0] <- 0
  Vmat[tindex[,2],] <- wz[tindex[,2], gindex[,'colW']]*x[tindex[,2], gindex[,'colX']]
  Vmat[timeZero,]   <- wz[timeZero, gindex[,'colW']]*x[timeZero, gindex[,'colX']]
  
  Lmat <- matrix(NA,nV,S)
  rownames(Lmat) <- rownames(gindex)
  loLmat <- hiLmat <- Lmat[,notOther]
  
  Lmat[ gindex[,c('rowL','colW')] ] <- lambda[ gindex[,c('rowG','colW')] ]
  
  
  lo <- hi <- Lmat*0
  lo[ gindex[,c('rowL','colW')] ] <- loLambda[ gindex[,c('rowG','colW')] ]
  hi[ gindex[,c('rowL','colW')] ] <- hiLambda[ gindex[,c('rowG','colW')] ]
  
  wL <- which(!is.na(Lmat[,notOther]),arr.ind=T)
  lo[is.na(lo)] <- 0
  hi[is.na(hi)] <- 0
  
  list(Lmat = Lmat, loLmat = lo[,notOther], hiLmat = hi[,notOther], wL = wL, 
       gindex = gindex, Vmat = Vmat)
}

.plotXbyY <- function(xdata, yy, ee, vname, xlim=range(xdata[,vname],na.rm=T)){
  
  plotj <- as.character(xdata[,'plot'])
  sitej <-  .splitNames(plotj)$vnam[,1]
  sitea <- sort(unique(sitej))
  years <- sort(unique(xdata[,'year']))
  nyr   <- length(years)
  labs  <- rep(1:12,nyr)
  week  <- (1:(nyr*12))*30
  midYr <- (1:(nyr*2))*365/2
  S     <- ncol(yy)
  ynames <- colnames(yy)
  
  mfrow <- .getPlotLayout(S)
  xaxt  <- 's'
  if(vname == 'JD'){
    xaxt <- 'n'
    xlim <- c(130,365*nyr)
  }
  
  par(mfrow=mfrow,bty='n',mar=c(4,4,1,1))
  
  for(s in 1:S){
    ys <- yy[,s]/ee[,s]
    plot(NULL,xlim=xlim,ylim=range(ys,na.rm=T),xlab=' ',ylab='count/14 day/m2',
         xaxt=xaxt)
    if(vname == 'JD'){
      abline(v=midYr,col='grey',lty=2)
      axis(1,at=week,labels=labs)
    } 
    for(k in 1:length(sitea)){
      for(j in 1:nyr){
        wk <- which(sitej == sitea[k] & xdata[,'year'] == years[j])
        if(length(wk) == 0)next
        lines((j-1)*365 + xdata[wk,vname],ys[wk])
      }
    }
    title(ynames[s])
  }
}

.pasteCols <- function(mm){
  tmp <- apply(mm,1,paste0,collapse='-')
  names(tmp) <- NULL
  tmp
}

.getURowCol <- function(mat){
  
  # mat is S by S
  
  rownames(mat) <- colnames(mat) <- NULL
  S <- nrow(mat)
  
  ww <- which(is.finite(mat),arr.ind=T)      #keep this
  wnames <- .pasteCols(ww)
  
  #unique rows/columns
  
  ar <- matrix(0,S,S)
  ar[ww] <- 1
  ar[ww[,c(2,1)]] <- 1
  wa <- which(ar == 1,arr.ind=T)
  wa <- wa[wa[,2] >= wa[,1],]
  wa <- wa[order(wa[,1],wa[,2]),]
  uindex <- wa
  
  un <- .pasteCols(wa)
  nu <- .pasteCols(wa[,c(2,1)])
  
  arow <- match(wnames,un)
  vrow <- match(wnames,nu)
  
  arow[is.na(arow)] <- vrow[is.na(arow)]
  rownames(ww) <- un[arow]
  wA <- cbind(arow,ww[,2])
  colnames(wA) <- c('rowA','toW')
  
  aindex <- cbind(wA,ww[,1])
  colnames(aindex)[3] <- 'fromW'
  
  Amat <- matrix(NA,nrow(uindex),S)
  Amat[wA] <- mat[ aindex[,c('fromW','toW')] ]
  rownames(Amat) <- un
  rownames(uindex) <- un
  
  list(uindex = uindex, Amat = Amat, wA = wA, aindex = aindex)
}

.invertCondKronecker <- function(sinv1, sinv2, cindex){
  
  # cindex - rows and columns to condition on in ns1*ns2
  # inverse of kronecker(s1,s2) = kronecker(sinv1,sinv2)
  
  ns1 <- nrow(sinv1)
  ns2 <- nrow(sinv2)
  n   <- ns1*ns2
  
  i1  <- c(1:n)[cindex]      #estimated
  i0  <- c(1:n)[-cindex]     #known
  
  #    sk  <- kronecker(s1,s2)
  ck  <- kronecker(sinv1,sinv2)
  sk  <- solveRcpp( .invertSigma(ck[i1,i1],REDUCT=F ) )
  
  sinverse <- ck[i0,i0] - ck[i0,i1]%*%sk%*%ck[i1,i0]
  sinverse
}

gjamFillMissingTimes <- function(xdata, ydata, edata, groups, times, 
                                 sequences=NULL, fillNA=T, fillTimes=T){
  
  # fill missing times, add initial time for prior
  # xdata, ydata, edata - x, y, effort
  # fillTimes - insert rows for missing times: integers between "times"
  # fillNA    - fill new rows in ydat with NA; otherwise fitted value
  # IMPORTANT - groups must uniquely defined 
  
  groupIndex <- xdata[,groups]
  if(is.factor(groupIndex))groupIndex <- as.character(groupIndex)
  allGroups  <- sort(unique(groupIndex))
  groupIndex <- match(groupIndex,allGroups)
  ngroups    <- length(allGroups)
  
  allTimes   <- sort(unique(xdata[,times]))
  allTimes   <- min(xdata[,times], na.rm=T):max(xdata[,times], na.rm=T)
  timeIndex  <- match(xdata[,times],allTimes)
  
  if(!fillTimes){
    timeIndex <- numeric(0)
    for(j in 1:ngroups){
      wj <- which(groupIndex == j)
      tj <- c(1:length(wj))
      timeIndex <- c(timeIndex,tj)
    }
  }
  
  xdata     <- cbind(groupIndex,timeIndex,xdata)
  timeZero <- numeric(0)
  
  if(!is.null(sequences)){
    if(!is.character(sequences) & !is.factor(sequences)){
      stop('sequences cannot be character or factor')
    }
    allSeq   <- sort(unique(xdata[,sequences]))
    seqIndex <- match(xdata[,sequences],allSeq)
    tord <- order(groupIndex, seqIndex, timeIndex, decreasing=F)
  } else{
    tord <- order(groupIndex, timeIndex, decreasing=F)
  }
  
  xtmp <- xdata[tord,]
  ytmp <- ydata[tord,]
  etmp <- edata[tord,]
  
  timeIndex  <- xtmp[,'timeIndex']
  groupIndex <- xtmp[,'groupIndex']
  
  wg <- which(colnames(xtmp) == 'groups')
  if(length(wg) == 0){
  } else{
    xtmp[,wg] <- groupIndex
  }
  wg <- which(colnames(xtmp) == 'times')
  if(length(wg) == 0){
  } else{
    xtmp[,wg] <- timeIndex
  }
  
  xtmp <- cbind(0,xtmp)
  
  notFactor <- !sapply(xtmp,is.factor)
  notChar   <- !sapply(xtmp,is.character)
  notFactor <- which(notFactor & notChar)
  
  for(j in 1:ngroups){
    
    wj <- which(xtmp$groupIndex == j)
    dj <- which( diff(xtmp[wj,times]) > 1 )
    
    timej <- xtmp[wj,'timeIndex']
    
    if(length(dj) == 0 & timej[1] == 0)next
    
    xtmp[wj,'timeIndex'] <- timej - xtmp[wj[1],'timeIndex'] + 1
    
    # initial time
    xnew <- xtmp[wj[1],]
    xnew[1,1] <- 1
    xnew[1,'timeIndex'] <- 0
    xnew[1,'timeIndex'] <- xtmp[wj[1],'timeIndex'] - 1
    xnew[1,times] <- xnew[1,times] - 1
    ynew <- ytmp[wj[1],]
    if(fillNA)ynew <- ynew + NA
    enew <- etmp[wj[1],]
    
    insert <- wj[1] - 1
    timeZero <- c(timeZero,insert+1)
    
    if(insert == 0){
      xtmp <- rbind(xnew,xtmp)
      ytmp <- rbind(ynew,ytmp)
      etmp <- rbind(enew,etmp)
    } else{
      others <- insert  + 1
      nn     <- nrow(xtmp)
      others <- others:nn
      xtmp <- rbind(xtmp[1:insert,],xnew,xtmp[others,])
      ytmp <- rbind(ytmp[1:insert,],ynew,ytmp[others,])
      etmp <- rbind(etmp[1:insert,],enew,etmp[others,])
    }
    
    wj <- which(xtmp[,groups] == allGroups[j])
    dj <- which( diff(xtmp[wj,'timeIndex']) > 1 )
    
    if(length(dj) > 0){
      
      for(k in 1:length(dj)){
        
        wj <- which(xtmp[,groups] == allGroups[j])
        dj  <- which( diff(xtmp[wj,'timeIndex']) > 1 )
        if(length(dj) == 0)break
        
        kd  <- wj[c(dj[1],dj[1]+1)]
        xdd <- xtmp[kd,]
        dt  <- diff(xtmp[kd,'timeIndex'])
        if(dt == 2)ts <- 1
        if(dt > 2) ts <- c(1:(dt-1))
        
        dm <- max( c(dt-1,1) )
        
        xnew <- xtmp[rep(kd[1],dm),]
        ynew <- ytmp[rep(kd[1],dm),]
        if(fillNA)ynew <- ynew + NA
        enew <- etmp[rep(kd[1],dm),]
        
        if(length(notFactor) > 0){
          nf    <- length(notFactor)
          nt    <- length(ts)
          xnot  <- as.matrix(xdd[,notFactor])
          slope <- matrix( apply(xnot,2,diff)/dt, nt, nf, byrow=T)
          xnew[,notFactor] <- xnew[,notFactor] + matrix(ts, nt, nf)*slope
        }
        
        xnew[,1] <- 1
        insert <- kd[1]
        others <- insert  + 1
        nn     <- nrow(xtmp)
        others <- others:nn
        xtmp <- rbind(xtmp[1:insert,],xnew,xtmp[others,])
        ytmp <- rbind(ytmp[1:insert,],ynew,ytmp[others,])
        etmp <- rbind(etmp[1:insert,],enew,etmp[others,])
      }
    }
  }
  
  timeLast <- timeZero[-1] - 1
  timeLast <- c(timeLast,nrow(xtmp))
  
  colnames(xtmp)[colnames(xtmp) == 'groupIndex'] <- 'groups'
  colnames(xtmp)[colnames(xtmp) == 'timeIndex'] <- 'times'
  
  noEffort <- which(rowSums(etmp,na.rm=T) == 0)
  noEffort <- noEffort[!noEffort %in% timeZero]
  
  rowInserts <- which(xtmp[,1] == 1)
  list(xdata = xtmp[,-1], ydata = as.matrix(ytmp), edata = etmp, 
       timeZero = timeZero, timeLast = timeLast,
       rowInserts = rowInserts, noEffort = noEffort)
}


.traitTables <- function(specs, traitTab, groups, types='CA', fill=T){
  
  # specs    - matrix of species names, at least 2 columns,
  #            for which a specByTrait table 
  #            is needed
  # traitTab - S x M spec by trait matrix of traits
  # groups   - data frame, S rows, columns are species, genus,...
  #            columns match specs matrix
  # types    - CA and CON will be means, OC will be modes
  # if(fill), then fill missing with trait means
  
  wi       <- which( !duplicated(groups[,1]) )
  groups   <- groups[wi,]  # one row per species
  traitTab <- traitTab[wi,]
  
  for(j in 2:ncol(specs)){
    specs[,j]  <- .cleanNames(as.character(specs[,j]))
    groups[,j] <- .cleanNames(as.character(groups[,j]))
  }
  
  ng        <- ncol(groups)
  
  #species traits
  ns <- nrow(specs)
  nt <- ncol(traitTab)
  
  stt <- matrix(0,ns,nt)
  
  i  <- match(specs[,1],groups[,1])
  wi <- which(is.finite(i))
  stt[wi,] <- as.matrix( traitTab[i[wi],] )
  rownames(stt) <- specs[,1]
  colnames(stt) <- colnames(traitTab)
  
  tabList <- list(stt)
  
  for(k in 2:ng){
    
    sk0    <- tabList[[k-1]]
    kn     <- rownames(sk0)
    
    allk <- unique(specs[,k])
    ii   <- match(specs[,k],allk)
    i    <- rep( ii, nt )
    j    <- rep( 1:nt, each=ns)
    
    nall <- length(allk)
    
    mm   <- matrix(0,nall,nt)
    sk0[sk0 == 0] <- NA
    
    stk <- .byGJAM(as.vector(sk0), i, j, mm,mm, fun='mean')
    colnames(stk) <- colnames(traitTab)
    rownames(stk) <- allk
    
    stk <- stk[specs[,k],]
    stk[is.finite(sk0)] <- sk0[is.finite(sk0)]
    
    tabList <- append(tabList, list(stk))
  }
  
  traitMeans <- colMeans(traitTab,na.rm=T)
  
  if(fill){
    
    ww <- which(tabList[[ng]] == 0, arr.ind=T)
    if(length(ww) > 0){
      tabList[[ng]][ww] <- traitMeans[ww[,2]]
    }
  }
  list(traitMeans = traitMeans, traitList = tabList, 
       specs = specs)
}

.combineFacLevels <- function(xfactor,fname=NULL, aname = 'reference', 
                              vminF=1){
  tmp <- as.character(xfactor)
  tmp[tmp %in% fname] <- aname
  tab <- table(tmp)
  wm  <- names(tab)[tab < vminF]
  tmp[tmp %in% wm] <- aname
  as.factor(tmp)
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}

.figure1 <- function(){
  sig <- .9
  mu  <- 3.1
  offset <- -2
  
  par(mfrow=c(1,2),bty='n',mar=c(6,5,3,.1))
  part <- c(0,2,3.3,4.9,6.6)
  w    <- seq(-1,7,length=500)
  dw   <- dnorm(w,mu,sig)
  dp   <- dw[ findInterval(part,w) ]
  pw   <- pnorm(part,mu,sig)
  pw[-1] <- diff(pw)
  
  plot(w,2*dw - .5,type='l',ylim=c(-.5,4),yaxt='n', 
       ylab= expression(paste(italic(y),'|(',italic(w),', ',bold(p),')',sep='')), 
       xlab= expression(paste(italic(w),'|(',bold(x)[i],', ',bold(beta),
                              ', ',bold(Sigma),')',sep='')), 
       xlim=c(offset,7), lwd=2)
  axis(2, at = c(0:5))
  
  db <- .15
  int <- 4
  
  polygon( c(w,rev(w)), 2*c(dw,w*0) - .5, col='grey', lwd=2)
  
  lines(c(-1,part[1]),c(0,0),lwd=2)
  
  for(j in 1:(length(part))){
    
    lines( part[j:(j+1)],c(j,j), lwd=3)
    
    ww <- which(w >= part[j] & w <= part[j+1])
    
    if(j == 3){
      w1 <- ww[1]
      w2 <- max(ww)
      arrows( mean(w[ww]), 2*max(dw[ww]) - .4, mean(w[ww]), 
              j - .4, angle=20,lwd=3, col = 'grey', length=.2)
      arrows( w[w1] - .5 , j , -.7, j , angle= 20, 
              lwd = 3, col='grey', length=.2)
      text( c(w[w1], w[w2]),c(3.3,3.3),
            expression(italic(p)[4], italic(p)[5]))
      text( w[w2] + .3,.6,expression( italic(w)[italic(is)] ))
      text( 0,3.5,expression( italic(y)[italic(is)] ))
    }
    
    coll <- 'white'
    if(j == int)coll <- 'grey'
    rect( offset, j - 1 - db, 2*pw[j] + offset, j - 1 + db, 
          col=coll, border='black', lwd=2)
  }
  
  ww <- which(w >= part[int - 1] & w <= part[int])
  abline(h = -.5, lwd = 2)
  
  title('a) Data generation',adj=0, font.main = 1, font.lab =1)
  
  plot(w,2*dw - .5,type='l',ylim=c(-.5,4), yaxt='n', 
       ylab= expression(italic(y)), 
       xlab= expression(paste(italic(w),'|(',italic(y),', ',bold(p),')',sep='')), 
       xlim=c(offset,7), lwd=2,col='grey')
  axis(2, at = c(0:5))
  
  lines(c(-1,part[1]),c(0,0),lwd=2)
  abline(h=-.5, lwd=2, col='grey')
  
  for(j in 1:(length(part))){
    
    lines( part[j:(j+1)],c(j,j), lwd=3)
    lines(part[c(j,j)],2*c(0,dp[j])-.5, col='grey')
    
    coll <- 'white'
    if(j -- int)coll <- 'grey'
    
    if(j == int){
      rect( offset, j - 1 - db, 2*pw[j] + offset, j - 1 + db,
            col='black', border='black')
    }
  }
  
  ww <- which(w >= part[int - 1] & w <= part[int])
  polygon( w[c(ww,rev(ww))], 2*c(dw[ww],ww*0) - .5, col='grey', lwd=2)
  
  arrows( mean(w[ww]),  int - 1.3,mean(w[ww]),  2*max(dw) - .5,
          angle=20,lwd=3, col = 'grey', length=.2)
  arrows( -.5,  int - 1, min(w[ww]) - .4, int - 1, angle= 20,
          lwd = 3, col='grey', length=.2)
  
  title('b) Inference',adj=0, font.main = 1, font.lab = 1)
}

.add2matrix <- function(values,xmat=NULL){
  
  #xmat   - n X ? matrix with one row, columns are integer values
  #values - length-n vector be added/slotted in to xvec
  
  if(is.null(xmat)){
    n    <- length(values)
    cc   <- sort(unique(values))
    xmat <- matrix(0,n,length(cc),dimnames = list(1:n,cc))
    xmat[ cbind( c(1:n),match(values,cc)) ] <- 1
    return(xmat)
  }
  
  n <- nrow(xmat)
  if(length(values) != n)stop('vector length must equal rows in xmat')
  
  all <- sort( unique( c(values,as.numeric(colnames(xmat))) ))
  nc       <- length(all)
  
  xnew <- matrix(0,n,nc,dimnames = list(1:n,all))
  xnew[,colnames(xmat)] <- xmat
  
  xnew[ cbind(c(1:n),match(values,all)) ] <- xnew[ cbind(c(1:n),match(values,all)) ] + 1
  xnew
}

.appendMatrix <- function(m1,m2,fill=NA,SORT=F,asNumbers=F){  
  
  # matches matrices by column names
  # asNumbers: if column heads are numbers and SORT, then sort numerically
  
  if(length(m1) == 0){
    if(is.matrix(m2)){
      m3 <- m2
    } else {
      m3 <- matrix(m2,nrow=1)
    }
    if( !is.null(names(m2)) )colnames(m3) <- names(m2)
    return(m3)
  }
  if(length(m2) == 0){
    if(!is.matrix(m1))m1 <- matrix(m1,nrow=1)
    return(m1)
  }
  if( is.vector(m1) | (length(m1) > 0 & !is.matrix(m1)) ){
    nn <- names(m1)
    if(is.null(nn))warning('cannot append matrix without names')
    m1 <- matrix(m1,1)
    colnames(m1) <- nn
  }  
  if( is.vector(m2) | (length(m2) > 0 & !is.matrix(m2)) ){
    nn <- names(m2)
    if(is.null(nn))warning('cannot append matrix without names')
    m2 <- matrix(m2,1)
    colnames(m2) <- nn
  }
  
  c1 <- colnames(m1)
  c2 <- colnames(m2)
  r1 <- rownames(m1)
  r2 <- rownames(m2)
  n1 <- nrow(m1)
  n2 <- nrow(m2)
  
  allc <-  unique( c(c1,c2) ) 
  if(SORT & !asNumbers)allc <- sort(allc)
  if(SORT & asNumbers){
    ac <- as.numeric(allc)
    allc <- as.character( sort(ac) )
  }
  
  nr <- n1 + n2
  nc <- length(allc)
  
  if(is.null(r1))r1 <- paste('r',c(1:n1),sep='-')
  if(is.null(r2))r2 <- paste('r',c((n1+1):nr),sep='-')
  new <- c(r1,r2)
  
  mat1 <- match(c1,allc)
  mat2 <- match(c2,allc)
  
  out <- matrix(fill,nr,nc)
  colnames(out) <- allc
  rownames(out) <- new
  
  out[1:n1,mat1] <- m1
  out[(n1+1):nr,mat2] <- m2
  out
}

.byIndex <- function(xx,INDICES,FUN,coerce=F,...){  
  
  # INDICES is list, each same length as  x
  
  # fun <- match.fun(FUN)
  
  nl <- length(INDICES)
  
  tmp  <-  unlist(by( as.vector(xx),INDICES,FUN,...) ) 
  nd   <- dim(tmp)
  tmp  <- array(tmp,dim=nd, dimnames=dimnames(tmp))
  
  tmp[is.na(tmp)] <- 0
  
  if(!coerce)return(tmp)
  
  dname <- dimnames(tmp)
  mk    <- rep(0,length(nd))
  
  for(k in 1:length(nd))mk[k] <- max(as.numeric(dimnames(tmp)[[k]]))
  
  wk <- which(mk > nd)
  if(length(wk) > 0){
    tnew  <- array(0,dim=mk)
    if(length(dim(tnew)) == 1)tnew <- matrix(tnew,dim(tnew),1)
    for(k in wk){
      newk <- c(1:mk[k])
      mat  <- match(dimnames(tmp)[[k]],newk)
      if(k == 1){
        tnew[mat,] <- tmp
        rownames(tnew) <- 1:nrow(tnew)
      }
      if(k == 2){
        tnew[,mat] <- tmp
        colnames(tnew) <- c(1:ncol(tnew))
      }
      tmp <- tnew
    }
  }
  tmp
}

.chains2density <- function(chainMat,labs=NULL,reverseM=F,varName=NULL,
                            cut=0){
  
  #assumes column names are varName or 'something_varname'
  #chainMat - MCMC output [samples,chains]
  
  chNames <- colnames(chainMat)
  
  if(!is.null(varName)){
    
    wc <- grep(varName,colnames(chainMat),fixed=T)
    if(length(wc) == 0)stop('varName not found in colnames(chainMat)')
    
    ww <- grep('_',colnames(chainMat),fixed=T)
    if(length(ww) > 0){
      tmp <- .splitNames(colnames(chainMat))$vnam
      wc  <- which(tmp[,2] == varName)
      if(length(wc) == 0)wc <- which(tmp[,1] == varName)
    }
    chainMat <- chainMat[,wc]
    if(!is.matrix(chainMat))chainMat <- matrix(chainMat,ncol=1)
    colnames(chainMat) <- chNames[wc]
  }
  
  nj <- ncol(chainMat)
  nd <- 512
  
  clab <- colnames(chainMat)
  if(is.null(labs) & !is.null(clab))labs <- clab
  if(is.null(labs) & is.null(clab)) labs <- paste('v',c(1:nj),sep='-')
  
  xt <- yt <- matrix(NA,nj,nd)
  rownames(xt) <- rownames(yt) <- labs
  
  xrange <- signif(range(chainMat),2)
  
  for(j in 1:nj){
    
    #   lj  <- labs[j]
    xj  <- chainMat[,j]
    tmp <- density(xj,n = nd, cut=cut, na.rm=T)
    xt[j,]  <- tmp$x
    yt[j,]  <- tmp$y
    
  }
  yymax <- max(yt,na.rm=T)
  
  if(reverseM){
    xt <- -t( apply(xt,1,rev) )
    yt <- t( apply(yt,1,rev) )
  }
  
  list(x = xt, y = yt, xrange = xrange, ymax = yymax, chainMat = chainMat)
}

.checkDesign <- function( x, intName='intercept', xflag=':', 
                          isFactor = character(0) ){  # 
  
  # xflag - indicates that variable is an interaction
  # isFactor - character vector of factor names returned if not supplied
  
  p <- ncol(x)
  
  if(ncol(x) < 3){
    return( list(VIF = 0, correlation = 1, rank = 2, p = 2, isFactor=isFactor) )
  }
  
  if(is.null(colnames(x))){
    colnames(x) <- paste('x',c(1:p),sep='_')
  }
  xrange      <- apply(x,2,range,na.rm=T)
  wi          <- which(xrange[1,] == 1 & xrange[2,] == 1)
  if(length(wi) > 0)colnames(x)[wi] <- 'intercept'
  
  wx <- grep(xflag,colnames(x))
  wi <- which(colnames(x) == 'intercept')
  wi <- unique(c(wi,wx))
  
  xname <- colnames(x)
  
  wmiss <- which(is.na(x),arr.ind=T)
  
  if(length(wmiss) > 0){
    rowTab <- table( table(wmiss[,1]) )
    colTab <- table(wmiss[,2])
  }
  
  VIF <- rep(NA,p)
  names(VIF) <- xname
  
  GETF <- F
  if(length(isFactor) > 0)GETF <- T
  
  for(k in 1:p){
    
    if(xname[k] %in% wi)next
    
    notk <- xname[xname != xname[k] & !xname %in% xname[wi]]
    ykk  <- x[,xname[k]]
    xkk  <- x[,notk,drop=F]
    
    wna <- which(is.na(ykk) | is.na(rowSums(xkk)))
    if(length(wna) > 0){
      ykk <- ykk[-wna]
      xkk <- xkk[-wna,]
    }
    
    ttt <- suppressWarnings( lm(ykk ~ xkk) )
    
    tkk <- suppressWarnings( summary(ttt)$adj.r.squared )
    VIF[k] <- 1/(1 - tkk)
    
    xu <- sort( unique(x[,k]) )
    tmp <- identical(c(0,1),xu)
    if(GETF)if(tmp)isFactor <- c(isFactor,xname[k])
  }
  
  VIF <- VIF[-wi] 
  
  corx <- cor(x[,-wi], use="complete.obs")
  if(length(wna) == 0){
    rankx <- qr(x)$rank
  } else {
    rankx <- qr(x[-wna,])$rank
  }
  corx[upper.tri(corx,diag=T)] <- NA
  
  findex <- rep(0,p)
  
  findex[xname %in% isFactor] <- 1
  
  designTable <- list('table' = rbind( round(VIF,2),findex[-wi],round(corx,2)) )
  rownames(designTable$table) <- c('VIF','factor',xname[-wi])
  
  designTable$table <- designTable$table[-3,]
  
  if(p == rankx)designTable$rank <- paste('full rank:',rankx,'= ncol(x)')
  if(p < rankx) designTable$rank <- paste('not full rank:',rankx,'< ncol(x)')
  
  list(VIF = round(VIF,2), correlation = round(corx,2), rank = rankx, p = p,
       isFactor = isFactor, designTable = designTable)
}

.fitText2Fig <- function(xx, width=T, fraction=1, cex.max=1){
  
  # returns cex to fit xx within fraction of the current plotting device
  # width - horizontal labels stacked vertically
  #!width - vertical labels plotted horizontally
  
  px <- par('pin')[1]
  py <- par('pin')[2]
  cl <- max( strwidth(xx, units='inches') )
  ch <- strheight(xx, units='inches')[1]*length(xx)  # ht of stacked vector
  
  if(width){              #horizontal labels stacked vertically
    xf <- fraction*px/cl
    yf <- fraction*py/ch
  } else {                #vertical labels plotted horizontally
    xf <- fraction*px/ch
    yf <- fraction*py/cl
  }
  
  cexx <- min(c(xf,yf))
  if(cexx > cex.max)cexx <- cex.max
  cexx
}

.cov2Dist <- function(sigma){ #distance induced by covariance
  
  n <- nrow(sigma)
  matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}

.distanceMatrix <- function(mat, DIST=F){
  
  # mat is n by m matrix
  # if DIST returns a m by m distance matrix, otherwise corr matrix
  
  if(isSymmetric(mat)){
    if( all(diag(mat) == 1) ){   #now a correlation matrix
      mmm1 <- mat
      if(DIST)mmm1 <- .cov2Dist(mat)
    } else {                      #now a covariance 
      if(DIST){
        mmm1 <- .cov2Dist( mat )
      } else {
        mmm1 <- cor(mat)
      }
    }
  } else  {     # not symmetric
    if(DIST){
      mmm1 <- .cov2Dist( cov(mat) )
    } else {
      mmm1 <- cor(mat)
    }
  }
  mmm1
}


.reorderMatrix <- function(mat, DIST, opt=NULL){
  
  # row and column order based on correlation or distance
  
  mmm <- .distanceMatrix(mat, DIST)
  imm <- .distanceMatrix(t(mat), DIST)
  
  if(is.null(opt)){
    clist <- list(PLOT=F, DIST = DIST)
  }else{
    clist <- opt
    clist$PLOT <- T
  }
  
  h1   <- .clusterPlot( imm, opt = clist)
  rord <- h1$corder
  h2   <- .clusterPlot( mmm, opt = clist )
  cord <- h2$corder
  
  list(rowOrder = rord, colOrder = cord, rowTree = h1, colTree = h2)
}

.clustMat <- function(mat, SYM){
  
  mrow <- mat
  DIST <- T
  if(SYM){
    if(all(diag(mrow) == 1)){
      DIST <- F
      mrow <- cor(mat)
      if(nrow(mrow) != nrow(mat))mrow <- cor(t(mat))
    }else{
      mrow <- .cov2Dist(cov(mat))
      if(nrow(mrow) != nrow(mat))mrow <- .cov2Dist(cov(t(mat)))
    }
  }else{
    mrow <- .cov2Dist(cov(mat))
    if(nrow(mrow) != nrow(mat))mrow <- .cov2Dist(cov(t(mat)))
  }
  list(cmat = mrow, DIST = DIST)
}

.clusterWithGrid <- function(mat1, mat2=NULL, opt, expand=1){
  
  #   layout: mat1 on left, mat2 (if given) on right
  # clusters: left & top or right & top
  #   expand: width of mat1 relative to mat2
  # if mat1/2 is symmetric can order only rows--stay symmetric
  # if DIST use distance, otherwise correlation
  
  leftClus <- rightClus <- 
    topClus1 <- topClus2 <- leftLab <- 
    rightLab <- topLab1 <- topLab2 <- lower1 <- diag1 <- lower2 <- 
    diag2 <- SYM1 <- SYM2 <- sameOrder <- FALSE
  colOrder1 <- colOrder2 <- rowOrder <- colCode1 <- colCode2 <- rowCode <- 
    slim1 <- slim2 <- horiz1 <- horiz2 <- vert1 <- vert2 <- NULL
  mainLeft <- main1 <- main2 <- ' '
  DIST1 <- DIST2 <- T
  ncluster <- 4
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  if(isSymmetric(mat1))SYM1 <- T
  
  doneLeft <- done1 <- done2 <- F
  
  nr  <- nrow(mat1)
  nc1 <- ncol(mat1)
  nc2 <- 0
  
  twoMat <- F
  if(!is.null(mat2)){
    if( min(dim(mat2)) < 2 )return()
    if(isSymmetric(mat2))SYM2 <- T
    twoMat <- T
    nc2 <- ncol(mat2)
    if(nrow(mat2) != nr)stop('matrices must have same no. rows')
  }
  cwide <- .15
  mg    <- .08
  lg    <- rg <- tg <- mg
  gg    <- .24
  
  if(leftLab) lg <- gg
  if(topLab1 | topLab2)  tg <- gg
  if(rightLab)rg <- gg
  
  xwide <- mg
  if(leftLab) xwide <- c(xwide,lg)
  if(leftClus)xwide <- c(xwide,cwide)
  
  xg <- .8
  if(twoMat){
    xg <- expand*xg*nc1/(nc1+nc2)
    xg <- c(xg,1 - xg)
  }
  xwide <- c(xwide,xg)
  
  if(rightClus)xwide <- c(xwide,cwide)
  if(rightLab) xwide <- c(xwide,rg)
  xwide <- c(xwide,mg)
  xloc <- cumsum(xwide)/sum(xwide)
  
  ywide <- c(mg,.8)
  if(topClus1 | topClus2)ywide <- c(ywide,cwide)
  if(topLab1 | topLab2) ywide <- c(ywide,tg)
  ywide <- c(ywide,mg)
  yloc  <- cumsum(ywide)/sum(ywide)
  
  if(is.null(rowCode)) rowCode  <- rep('black',nr)
  if(is.null(colCode1))colCode1 <- rep('black',nc1)
  if(is.null(colCode2))colCode2 <- rep('black',nc2)
  
  tmp   <- .clustMat(mat1, SYM1)
  m1Row <- tmp$cmat
  DIST1 <- tmp$DIST
  tmp <- .clustMat(t(mat1), SYM1)
  m1Col <- tmp$cmat
  
  
  if(is.null(rowOrder)){
    if(nrow(m1Row) != nrow(mat1))m1Row <- cor(t(mat1))
    if(ncluster > nrow(m1Row)/2)ncluster <- 2
    copt <- list( PLOT=F, DIST = DIST1, ncluster=ncluster )
    tmp  <- .clusterPlot( m1Row, copt)
    clus <- tmp$clusterIndex
    cord <- tmp$corder
    rowClust <- clus[cord]
    rowOrder <- cord
  }
  
  if(is.null(colOrder1)){
    if(SYM1){
      colOrder1 <- rowOrder
      colClust1 <- rowClust
    }else{
      copt      <- list( PLOT=F, DIST = DIST1, ncluster=ncluster )
      tmp  <- .clusterPlot( m1Col, copt)
      clus <- tmp$clusterIndex
      cord <- tmp$corder
      colClust1 <- clus[cord]
      colOrder1 <- cord
    }
  }
  
  if(twoMat){
    tmp   <- .clustMat(t(mat2), SYM2)
    m2Col <- tmp$cmat
    DIST2 <- tmp$DIST
    
    if(is.null(colOrder2)){
      if(sameOrder){
        colOrder2 <- rowOrder
        m2Col <- t(mat2)
      }else{
        copt <- list( PLOT=F, DIST = DIST2 )
        tmp  <- .clusterPlot( m2Col, copt)
        if(is.null(tmp)){
          colOrder2 <- 1:nrow(m2Col)
        }else{
          clus <- tmp$clusterIndex
          cord <- tmp$corder
          colClust2 <- clus[cord]
          colOrder2 <- cord
        }
      }
    }
  }
  
  rowLabs <- rownames(mat1)[rowOrder]
  
  
  #######################
  NEW <- add <- F
  xi <- 0:1
  yi <- 1:2
  
  ##### lab panel -- bottom to top
  if(leftLab){  
    xi <- xi + 1
    par(plt=c(xloc[xi],yloc[yi]),bty='n', new=NEW)
    
    plot(NULL,col='white',xlim=c(0,1),ylim=c(0,nr),
         xaxt='n',yaxt='n',xlab='',ylab='')
    xl  <- rep(1,nr)
    yl  <- c(1:nr)*nr/diff(par('usr')[3:4])
    cex <- .fitText2Fig(rowLabs,fraction=.96)
    text( xl,yl, rowLabs ,pos=2,cex=cex, 
          col = rowCode[rowOrder])
    NEW <- add <- T
    mtext(mainLeft,2)
    doneLeft <- T
  }
  
  #### cluster panel
  if(leftClus){
    xi <- xi + 1
    par(plt=c(xloc[xi],yloc[yi]),bty='n',  new=NEW)
    copt <- list( main=' ',cex=.2, colCode=rowCode, ncluster=ncluster,
                  LABELS = F, horiz=T, noaxis=T, DIST=DIST1 )
    tmp <- .clusterPlot( m1Row, copt)
    clus <- tmp$clusterIndex
    cord <- tmp$corder
    
    rowClust <- clus[cord]
    
    NEW <- add <- T
    if(!doneLeft)mtext(mainLeft,2)
    doneLeft <- T
  }
  
  ######## first grid plot
  
  xi <- xi + 1
  yz <- yi
  
  if(topClus1){
    yz <- yz + 1
    par(plt=c(xloc[xi],yloc[yz]),bty='n',new=NEW)
    
    copt <- list( main=' ', colCode=colCode1, DIST = DIST1,     
                  LABELS = F, horiz=F, noaxis=T, add=T )
    tmp <- .clusterPlot( m1Col ,copt)
    NEW <- add <- T
    if(!topLab1){
      mtext(main1,3)
      done1 <- T
    }
  }
  
  par(plt=c(xloc[xi],yloc[yi]), bty='n', new=NEW)
  if(is.null(slim1))slim1 = quantile(mat1,c(.01,.99)) ######
  slim1  <- signif(slim1,1)
  
  tmp    <- .colorSequence(slim1)
  scale  <- tmp$scale
  colseq <- tmp$colseq
  
  ww    <- as.matrix(expand.grid(c(1:nr),c(1:nc1)))  # reverse order
  # mt    <- t(apply(mat1[rowOrder,colOrder1],1,rev)) ###########
  mt <- mat1[rev(rowOrder),colOrder1]
  
  win <- which(ww[,1] >= ww[,2])
  
  mask <- lower.tri(mt,diag=!diag1)
  mask <- apply(mask,2,rev)
  if(lower1){
    if(min(scale) > 0 | max(scale) < 0){mt[mask] <- mean(scale)
    }else{ mt[mask] <- 0 }
  }
  
  icol <- findInterval(mt[ww],scale,all.inside=T)
  coli <- colseq[icol]
  
  xlim=c(range(ww[,2])); xlim[2] <- xlim[2] + 1
  ylim=c(range(ww[,1])); ylim[2] <- ylim[2] + 1
  
  sides <- cbind( rep(1,nrow(ww)), rep(1,nrow(ww)) )
  plot(NULL,cex=.1,xlab=' ',ylab=' ', col='white',
       xaxt='n',yaxt='n', xlim=xlim, ylim=ylim)
  
  symbols(ww[,2] + .5,nr - ww[,1] + 1 + .5,rectangles=sides,
          fg=coli,bg=coli,inches=F, xlab=' ',ylab=' ',
          xaxt='n',yaxt='n', add=T)
  
  if(!is.null(horiz1)){
    #  cut <- which(diff(horiz1[rowOrder]) != 0) + 1
    
    cut <- which(diff(rowClust) != 0) + 1
    
    ncc <- length(cut)
    for(i in 1:ncc){
      lines(c(0,cut[i]-2),cut[c(i,i)],lty=2)
    }
    text(rep(1,ncc),cut,2:(ncc+1),pos=3)
  }
  
  
  
  if(!is.null(vert1)){
    #  cut <- which(diff(vert1[colOrder1]) != 0) + .5
    
    cut <- which(diff(colClust1) != 0) + .5
    
    ncc <- length(cut)
    for(i in 1:ncc){
      lines(cut[c(i,i)],c(cut[i]+2,nc1),lty=2)
    }
    text(cut,rep(nc1,ncc),2:(ncc+1),pos=4)
  }
  
  NEW <- add <- T
  if(!doneLeft)mtext(mainLeft,2)
  doneLeft <- T
  
  if(topLab1){
    #   if(isSymmetric(mat1))colCode1 <- rowCode
    yz <- yz + 1
    par(plt=c(xloc[xi], yloc[yz]),bty='n', new=NEW)
    plot(c(0,0),c(0,0),col='white',xlim=c(1,nc1) ,ylim=c(0,1),
         xaxt='n',yaxt='n',xlab='',ylab='')
    yl <- rep(0,nc1)
    xl <- .99*c(1:nc1)*(nc1-1)/diff(par('usr')[1:2])
    #   xl <- .95*c(1:nc1)*nc1/diff(par('usr')[1:2])
    cex <- .fitText2Fig(colnames(m1Col), 
                        width=F, fraction=.95)
    text( xl - .1,yl,colnames(m1Col)[colOrder1],pos=4,cex=cex,srt=90,
          col=colCode1[colOrder1])
  }
  if(!done1)mtext(main1,3)
  
  #color scale
  par(plt=c(xloc[xi],c(.3*yloc[yi[1]],yloc[yi[1]])), bty='n', new=NEW)
  lx1 <- .3
  lx2 <- .7
  
  lx <- seq(lx1,lx2,length=length(scale))
  wx <- diff(lx)[1]
  ly <- lx*0 + .3*yloc[yi[1]]
  rx <- cbind(lx*0 + wx, ly*0 + .7*diff(yloc[yi]))
  symbols(lx,ly,rectangles=rx,fg=colseq,bg=colseq,xaxt='n',yaxt='n',
          xlab='',ylab='',xlim=c(0,1),ylim=c(0,yloc[yi[1]]))
  text(lx[1],ly[1],slim1[1],pos=2, cex=.9)    
  text(lx[length(lx)],ly[1],slim1[2],pos=4, cex=.9)  
  
  ######## 2nd grid plot
  
  if(twoMat){
    
    xi <- xi + 1
    yz <- yi
    
    if( topClus2 ){
      yz <- yz + 1
      par(plt=c(xloc[xi],yloc[yz]),bty='n',new=NEW)
      
      copt <- list( main=' ', LABELS = F,
                    colCode=colCode2, horiz=F, 
                    noaxis=T, add=T, DIST=DIST2 )
      ttt <- .clusterPlot( m2Col, copt)
      
      #    m2 <- apply(mat1[rowOrder,colOrder1],1,rev)
      
      if(!topLab2){
        mtext(main2,3)
        done2 <- T
      }
    }
    
    par(plt=c(xloc[xi],yloc[yi]), bty='n', new=T)
    if(is.null(slim2))slim2 = quantile(mat2,c(.01,.99))
    slim2  <- signif(slim2,1)
    
    tmp <- .colorSequence(slim2)
    scale  <- tmp$scale
    colseq <- tmp$colseq
    
    ww    <- as.matrix(expand.grid(c(1:nr),c(1:nc2)))  # note reverse order
    mt    <- mat2[rev(rowOrder),colOrder2]
    
    if(lower2){
      mask <- lower.tri(mt,diag=!diag1)
      mask <- apply(mask,2,rev)
      mt[mask] <- 0
      if(min(scale) > 0 | max(scale) < 0){
        mt[mask] <- mean(scale)
      }else{ mt[mask] <- 0 }
    }
    
    icol <- findInterval(mt[ww],scale,all.inside=T)
    coli <- colseq[icol]
    
    xlim=c(range(ww[,2])); xlim[2] <- xlim[2] + 1
    ylim=c(range(ww[,1])); ylim[2] <- ylim[2] + 1
    
    sides <- cbind( rep(1,nrow(ww)), rep(1,nrow(ww)) )
    plot(0,0,cex=.1,xlab=' ',ylab=' ',
         col='white',xaxt='n',yaxt='n', xlim=xlim, ylim=ylim)
    
    symbols(ww[,2] + .5,nr - ww[,1] + 1 + .5, rectangles=sides,
            fg=coli, bg=coli, inches=F, xlab=' ',ylab=' ',
            xaxt='n', yaxt='n', add=T)
    
    if(!is.null(horiz2)){
      
      #  cut <- which(diff(horiz2[rowOrder]) != 0) + 1
      
      cut <- which(diff(rowClust) != 0) + 1
      
      ncc <- length(cut)
      for(i in 1:ncc){
        xmm <- c(0,cut[i]-2)
        if(!lower2)xmm[2] <- nc2 + 1
        lines(xmm,cut[c(i,i)],lty=2)
      }
      if(lower2) text(rep(1,ncc),cut,2:(ncc+1),pos=3)
      if(!lower2)text(rep(nc2+1,ncc),cut,2:(ncc+1),pos=3)
    }
    if(!is.null(vert2)){
      cut <- which(diff(vert2[colOrder2]) != 0) + .5
      ncc <- length(cut)
      for(i in 1:ncc){
        lines(cut[c(i,i)],c(cut[i]+2,nc1),lty=2)
      }
      text(cut,rep(nc1,ncc),2:(ncc+1),pos=4)
    }
    
    if(topLab2){
      yz <- yz + 1
      par(plt=c(xloc[xi],yloc[yz]),bty='n', new=NEW)
      plot(c(0,0),c(0,0),col='white',xlim=c(1,nc2),ylim=c(0,1),
           xaxt='n',yaxt='n',xlab='',ylab='')
      yl  <- rep(0,nc2)
      xl  <- c(1:nc2)*(nc2-1)/diff(par('usr')[1:2])
      cex <- .fitText2Fig(colnames(m2Col),width=F, fraction=.95)
      text( xl - .05,yl,colnames(m2Col)[colOrder2],pos=4,cex=cex,srt=90, 
            col=colCode2[colOrder2])
    }
    if(!done2)mtext(main2,3)
  }
  
  par(plt=c(xloc[xi],c(.3*yloc[yi[1]],yloc[yi[1]])), bty='n', new=NEW)
  
  lx1 <- .3
  lx2 <- .7
  
  lx <- seq(lx1,lx2,length=length(scale))
  wx <- diff(lx)[1]
  ly <- lx*0 + .3*yloc[yi[1]]
  rx <- cbind(lx*0 + wx, ly*0 + .7*diff(yloc[yi]))
  symbols(lx,ly,rectangles=rx,fg=colseq,bg=colseq,xaxt='n',yaxt='n',
          xlab='',ylab='',xlim=c(0,1),ylim=c(0,yloc[yi[1]]))
  text(lx[1],ly[1],slim2[1],pos=2, cex=.9)    
  text(lx[length(lx)],ly[1],slim2[2],pos=4, cex=.9)  
  
  if(rightClus){
    xi <- xi + 1
    par(plt=c(xloc[xi], yloc[yi]), bty='n', mgp=c(3,1,0), new=NEW)
    mmm <- .distanceMatrix(t(mat2), DIST1)
    copt <- list( main=' ',cex=.2, REV=T,
                  LABELS = F,horiz=T, noaxis=T )
    tmp <- .clusterPlot( mmm , copt)
  }
  
  if(rightLab){
    xi <- xi + 1
    par(plt=c(xloc[xi],yloc[yi]),bty='n', new=NEW)
    plot(c(0,0),c(0,0),col='white',xlim=range(c(0,1)),ylim=c(0,nr),
         xaxt='n',yaxt='n',xlab='',ylab='')
    xl <- rep(0,nr)
    yl <- c(1:nr)*nr/diff(par('usr')[3:4])
    cex <- .fitText2Fig(rownames(m1Row),fraction=.8)
    text( xl,yl,rev( rownames(m1Row) ),pos=4,cex=cex,
          col=rev(rowCode[rowOrder]))
  }
}


.clusterPlot <- function(dmat, opt = NULL){
  
  main <- xlab <- ' '
  method <- 'complete' 
  cex <- 1; ncluster <- 2; textSize <- 1
  add <- REV <- reverse <- noaxis <- DIST <- F
  xlim <- colCode <- NULL 
  horiz <- LABELS <- PLOT <- T
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  #dmat is a correlation matrix or distance matrix
  
  getTreeLabs <- function ( tree ){ #left to right or bottom to top
    
    getL <- function(tree_node) {
      if(is.leaf(tree_node)) 
        attr(tree_node, 'label')
    }
    unlist( dendrapply(tree, getL) )
  }
  
  
  # if(!LABELS) rownames(dmat) <- colnames(dmat) <- NULL
  nr   <- nrow(dmat)
  nn   <- nrow(dmat)
  
  if(min(c(nr,nn)) < 3)return()
  
  if(DIST){
    if(!isSymmetric(dmat))dmat <- dist(dmat)
    diss <- as.dist( dmat )
  }else{
    diss <- as.dist(.cov2Dist(dmat))
  }
  
  htree  <- hclust(diss,method)
  ctmp   <- cutree(htree,k=1:ncluster)
  
  wclus <- ctmp[,ncluster]
  clusterCol <- NULL
  
  clusterIndex <- ctmp[,ncluster]
  clusterList <- character(0)
  
  notLab <- F
  if(is.null(colCode)){
    colF   <- colorRampPalette(c('black','blue','orange','brown','red'))
    mycols <- colF(ncluster)
    notLab <- T
    colCode <- mycols[ctmp[,ncluster]]
    names(colCode) <- rownames(ctmp)
  }
  col.lab <- colCode
  if(!LABELS)col.lab <- rep('white',length(colCode))
  
  colLab <- function(n) {
    
    if(is.leaf(n)) {
      a <- attributes(n)
      attr(n, "nodePar") <- c(a$nodePar, 
                              list(col = col.lab[n[1]],
                                   lab.col = col.lab[n[1]]))
    }
    n
  }
  tdendro <- as.dendrogram(htree)
  if(reverse)tdendro <- rev(tdendro)
  dL      <- dendrapply(tdendro,colLab)
  
  tlab <- getTreeLabs(tdendro)
  corder <- match(tlab,colnames(dmat))
  names(corder) <- colnames(dmat)[corder]
  
  nodePar <- list(cex = .1, lab.cex=textSize)
  leafLab         <- "textlike"
  nodePar$leaflab <-  leafLab
  
  if(!PLOT){
    return(  invisible(list( clusterList = clusterList, colCode = colCode, 
                             clusterIndex = clusterIndex,
                             corder = corder) ) )
  }
  
  if(horiz){
    if(is.null(xlim))xlim <- c(attr(dL,'height'),0)
    if(REV)xlim <- rev(xlim)
  }
  
  axes <- T
  if(noaxis)axes <- F
  new <- F
  if(add)new <- T
  
  tmp <- plot( dL,nodePar=nodePar, horiz=horiz, xlim=xlim, 
               axes = axes)
  if(!LABELS & !notLab){
    
    col <- colCode[corder]
    pvar <- par('usr')
    
    wi <- abs(diff(pvar[1:2])/10)
    hi <- abs(diff(pvar[3:4])/10)
    
    if(horiz){
      xx <- rep(pvar[2],nn)
      yy <- 1:nn
      rec <- cbind( rep(wi,nn), rep(1,nn) )
      symbols(xx,yy,rectangles=rec,fg=col, bg=col, inches=F, add=T)
    } else {
      xx <- 1:nn
      yy <- rep(pvar[3],nn)
      rec <- cbind( rep(1,nn), rep(hi,nn) )
      symbols(xx,yy,rectangles=rec,fg=col, bg=col, inches=F, add=T)
    }
  }
  
  title(main)
  
  invisible(list( clusterList = clusterList, colCode = colCode, 
                  clusterIndex = clusterIndex,
                  corder = corder) )
  
}

.colorLegend <- function(xx,yy,ytick=NULL,
                         scale=seq(yy[1],yy[2],length=length(cols)),
                         cols,labside='right', text.col=NULL,
                         bg=NULL,endLabels=NULL){  
  # xx = (x1,x2), y = (y1,y2)
  # bg is color of border
  
  nn <- length(scale)
  ys <- seq(yy[1],yy[2],length=nn)
  
  for(j in 1:(length(scale)-1)){
    rect(xx[1],ys[j],xx[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(bg))rect(xx[1],yy[1],xx[2],yy[2],border=bg,lwd=3)
  if(!is.null(ytick)){
    
    ys <- diff(yy)/diff(range(ytick))*ytick
    yt <- ys - min(ys) + yy[1]
    
    for(j in 1:length(yt)){
      lines(xx,yt[c(j,j)])
    }
  }
  if(!is.null(endLabels)){ 
    cx <- cols[c(1,nn)]
    if(!is.null(text.col))cx <- text.col
    if(!is.null(text.col))cx <- text.col
    if(labside == 'right')text(diff(xx)+c(xx[2],xx[2]),yy,endLabels,col=cx)
    if(labside == 'left')text(c(xx[1],xx[1]),yy,endLabels,pos=2,col=cx)
  }
}

.capFirstLetter <- function(xx) {     
  
  #capiltalize first letter of every word
  
  s <- unlist(strsplit(xx, " "))
  s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

.lowerFirstLetter <- function(xx){
  s <- unlist(strsplit(xx, " "))
  s <- paste(tolower(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

.colorSequence <- function(slim, colorGrad=NULL, ncol=200){  
  
  # generate color sequence with white for zero
  # slim is scale from min to max
  # used in .corPlot
  
  if(is.null(colorGrad)){
    colorSeq <- c('darkblue','darkblue','blue',
                  'green','white',
                  'yellow','red','brown','brown')
    colorGrad   <- colorRampPalette(colorSeq)
  }
  
  colseq <- colorGrad(ncol)
  
  if(slim[1] < 0 & slim[2] > 0){  #put white at zero
    dp <- slim[2] - 0
    dm <- 0 - slim[1]
    ncol <- 200
    
    colseq <- colorGrad(ncol)
    if(dp < dm)colseq <- colseq[101 + c(-100:round(dp/dm*100))]
    if(dp > dm)colseq <- colseq[ round((1 - dm/dp)*100):200 ]
    ncol  <- length(colseq)
  }
  scale <- seq(slim[1],slim[2],length.out=ncol)
  return( list(colseq = colseq, scale = scale ) )
}

.corPlot <- function(cmat,slim=NULL,PDIAG=F,plotScale=1,
                     makeColor=NULL,textSize=NULL,
                     textCol = rep('black',nrow(cmat)), 
                     CORLINES=T,tri='lower',colorGrad = NULL,
                     cex=1, SPECLABS = T, squarePlot = T,LEGEND = T,
                     widex=5.5,widey=6.5,add=F,new=T){  
  # correlation or covariance matrix
  # makeColor - list of matrices of indices for boxes
  #   names of matrices are colors
  # if(PDIAG)diag(cmat) <- 0
  # tri - 'lower','upper', or 'both'
  # colorGrad - constructed with colorRampPalette()
  # squarePlot makes symbols square
  # new means NOT NEW 
  
  if(is.null(slim))slim = quantile(cmat,c(.01,.99))
  slim  <- signif(slim,1)
  
  if(tri == 'upper')cmat[lower.tri(cmat)] <- 0
  if(tri == 'lower')cmat[upper.tri(cmat)] <- 0
  
  dy  <- nrow(cmat)
  dx  <- ncol(cmat)
  d <- dx
  xtext <- rep(c(1,100),dx/2)
  if(length(xtext) < d)xtext <- c(xtext,1)
  
  if(d < 20)xtext <- xtext*0 + 1
  
  xtext <- xtext*0 + 1
  
  if(!is.null(colorGrad)){
    ncol  <- 200
    colseq <- colorGrad(ncol)
    scale  <- seq(slim[1],slim[2],length.out=ncol)
  } else {
    tmp <- .colorSequence(slim, colorGrad)
    scale  <- tmp$scale
    colseq <- tmp$colseq
  }
  
  ww   <- as.matrix(expand.grid(c(1:dy),c(1:dx)))  # note reverse order
  
  if(tri == 'upper'){
    ww  <- ww[ww[,1] <= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }
  if(tri == 'lower'){
    ww  <- ww[ww[,1] >= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }
  
  icol <- findInterval(cmat[ww],scale,all.inside=T)
  coli <- colseq[icol]
  
  if(PDIAG)coli[ww[,1] == ww[,2]] <- 'white'
  
  ss <- max(c(dx,dy))/5/plotScale
  
  if(squarePlot).mapSetup(c(0,dx),c(0,dy),scale=ss,
                          widex=widex,widey=widey)
  
  if(squarePlot){
    symbols(ww[,2],dy - ww[,1] + 1,squares=rep(1,nrow(ww)),
            xlim=c(0,dx+4),ylim=c(0,dy+4),
            fg=coli,bg=coli,inches=F,xlab=' ',ylab=' ',xaxt='n',yaxt='n',
            add=add)
  } else {
    sides <- cbind( rep(1,nrow(ww)), rep(1,nrow(ww)) )
    symbols(ww[,2],dy - ww[,1] + 1,rectangles=sides,
            xlim=c(0,dx+4),ylim=c(0,dy+4),
            fg=coli,bg=coli,inches=F,xlab=' ',ylab=' ',xaxt='n',yaxt='n',
            add=add)
  }
  
  if(!is.null(makeColor)){
    
    for(k in 1:length(makeColor)){
      mm <- makeColor[[k]]
      if(length(mm) == 0)next
      if(tri == 'upper')mm <- mm[mm[,1] <= mm[,2],]
      if(tri == 'lower')mm <- mm[mm[,1] >= mm[,2],]
      ss <- matrix(0,dy,dx)
      ss[mm] <- 1
      wk <- which(ss[ww] == 1)
      ccc <- names(makeColor)[[k]]
      symbols(ww[wk,2],dy - ww[wk,1]+1,squares=rep(1,length(wk)),
              fg=ccc,bg=ccc,inches=F,xaxt='n',yaxt='n',add=T)
    }
  }
  
  ncolor <- length(unique(textCol))
  
  ll <- 1/d + 1
  
  if(tri == 'lower'){
    for(kk in 1:d){
      kb <- kk - .5
      ke <- d - kk + .5
      
      if(CORLINES){
        if(kk <= d)lines(c(kb,kb),c(0,ke),col='grey',lwd=1.5) #vert
        if(kk > 1){
          lines( c(.5,kb),c(ke,ke),col='grey',lwd=1.5)        #horizontal
          lines(c(kb,kb+.5),c(ke,ke+.5),col='grey',lwd=1.5)   #diagonal
        }
      }
      if(!SPECLABS & ncolor > 1){
        xp <- c(kb, kb, kb + ll + .5, kb + ll + 1.5, kb + 1)
        yp <- c(ke, ke + 1, ke + ll + 1.5, ke + ll + .5, ke)
        polygon(xp, yp, border = textCol[kk], col = textCol[kk])
      }
    }
  }
  rect(0,-1,d+1,.5,col='white',border=NA)
  
  if(is.null(textSize))textSize <- exp(-.02*ncol(cmat))
  labels   <- rev(rownames(cmat))
  if(!SPECLABS)labels <- F
  
  if(tri == 'lower' & SPECLABS)text( c(d:1)+.1*xtext, c(1:d)+.5, 
                                     rev(colnames(cmat)),pos=4,srt=45,
                                     col = rev(textCol), cex=textSize)
  
  if(tri == 'both'){
    labels   <- rev(rownames(cmat))
    par(las = 1)
    
    .yaxisHorizLabs( labels, at=c(1:length(labels)), xshift=.05,
                     col = textCol, pos=2)
    
    par(las = 0)
    
    if(SPECLABS){
      text( c(dx:1)-.1*xtext, xtext*0+dy+.8, rev(colnames(cmat)),
            pos=4, srt=55, col = rev(textCol), cex=textSize)
    } else {
      sides <- cbind( rep(1,dx),rep(1/dy,dx) )
      symbols(1:dx,rep(1+dy,dx),rectangles=sides,
              fg=textCol,bg=textCol,
              add=T)
    } 
    
  }
  
  labside <- 'left'
  
  wk <- which(scale >= slim[1] & scale <= slim[2]) 
  
  px <- par('usr')
  xs <- .01*diff(px[1:2])
  midx <- .95*mean( c(dx,px[2]) )
  
  yx <- c(.2*dy,.2*dy + .35*dy)
  
  if(LEGEND).colorLegend(c(midx-xs,midx+xs),yx,ytick=c(slim[1],0,slim[2]),
                         scale[wk],cols=colseq[wk],labside=labside,
                         endLabels=range(slim),text.col='black')
}

.cor2Cov <- function(sigvec,cormat){ 
  
  #correlation matrix and variance vector to covariance
  
  d <- length(sigvec)
  s <- matrix(sigvec,d,d)
  cormat*sqrt(s*t(s))
}

.cov2Cor <- function(covmat, covInv = NULL){  
  
  # covariance matrix to correlation matrix
  # if covInv provided, return inverse correlation matrix
  
  d    <- nrow(covmat)
  di   <- diag(covmat)
  s    <- matrix(di,d,d)
  cc   <- covmat/sqrt(s*t(s))
  
  if(!is.null(covInv)){
    dc <- diag(sqrt(di))
    ci <- dc%*%covInv%*%dc
    return(ci)
  }
  cc
}

.cov2Dist <- function(sigma){ 
  
  #distance induced by covariance
  
  n <- nrow(sigma)
  matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}

.dMVN <- function(xx,mu,smat=NULL,sinv=NULL,log=F){ 
  
  #MVN density for mean 0
  
  if(!is.matrix(xx))xx <- matrix(xx,1)
  if(!is.matrix(mu))mu <- matrix(mu,1)
  
  tmp <- try( dmvnormRcpp(xx, mu, smat, logd=log),silent=T  )     
  if( !inherits(tmp,'try-error') )return(tmp)
  
  xx <- xx - mu
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = T)$values
    logd    <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(smat),T)
    if(inherits(testv,'try-error')){
      tiny  <- min(abs(xx))/100 + 1e-5
      smat  <- smat + diag(diag(smat + tiny))
      testv <- try(chol(smat),T)
    }
    covm    <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev      <- eigen(smat, only.values = T)$values 
    logd    <- sum(log( ev ))
  }
  
  z <- -(ncol(xx) * log(2 * pi) + logd + distval)/2
  if(!log)z <- exp(z)
  z
}

.directIndirectCoeffs <- function( snames, xvector, chains, MEAN = T,
                                   factorList = NULL, keepNames, omitY,
                                   sdScaleY = F, sdScaleX, standX, 
                                   otherpar = NULL, REDUCT = F, ng, burnin,
                                   nsim = 50){
  
  # if MEAN, then use means, otherwise median
  # indirect do not change with x, can choose not to calculate
  #          - a list of vectors, one for each multilevel factor, 
  #            where hostNames appear in colnames of bchain
  #indirFrom - effect from all others
  #indirTo   - effect on all others
  
  if(is.matrix(xvector))
    stop('xvector must be a row vector with variable names')
  
  xnames <- names(xvector)
  
  N      <- otherpar$N
  r      <- otherpar$r
  bchain <- chains$bgibbs
  schain <- chains$sgibbs
  sigErrGibbs <- kchain <- NULL
  if(REDUCT){
    kchain      <- chains$kgibbs
    sigErrGibbs <- chains$sigErrGibbs
  }
  
  ns <- nsim
  simIndex <- sample(burnin:ng,ns,replace=T)
  
  if(sdScaleY){
    tmp <- .expandSigmaChains(snames, sgibbs = schain, otherpar = otherpar, 
                              simIndex = simIndex, sigErrGibbs, kchain, 
                              REDUCT)
    
    if(REDUCT)kchain <- kchain[simIndex,]
    schain <- schain[simIndex,]          # not standardized
    sigErrGibbs <- sigErrGibbs[simIndex]
  } else {
    bchain <- bchain[simIndex,]
    schain <- schain[simIndex,]
  }
  
  if(length(factorList) > 0){
    factorNames <- factorList
    for(j in 1:length(factorList)){
      tmp <- matrix( unlist(strsplit(factorList[[j]],names(factorList)[j])),
                     ncol=2,byrow=T)[,2]
      tmp[nchar(tmp) == 0] <- paste(names(factorList)[j],c(1:length(tmp)),
                                    sep='')
      factorNames[[j]] <- tmp
    }
  }
  
  S <- S1 <- length(snames)
  sindex <- c(1:S)
  knames <- snames
  
  nc <- nrow(bchain)
  
  gs <- 1:nrow(bchain)
  
  if(length(omitY) > 0){
    wob <- grep(paste(omitY,collapse="|"),colnames(bchain))
    bchain[,wob] <- 0
    sindex <- sindex[!snames %in% omitY]
    knames <- snames[sindex]
    S1     <- length(knames)
  }
  
  nspec <- length(snames)
  
  ww   <- grep(':',xnames)
  main <- xnames
  if(length(ww) > 0)main <- xnames[-ww]
  main <- main[main != 'intercept']
  int  <- unique( unlist( strsplit(xnames[ww],':') ) ) 
  
  mainEffect <- matrix(NA,nspec,length(main))
  colnames(mainEffect) <- main
  rownames(mainEffect) <- snames
  intEffect  <- dirEffect <- indEffectTo <- mainEffect
  mainSd <- dirSd <- intSd <- indSdTo <- mainEffect 
  
  maxg <- length(main)*length(sindex)*length(gs)
  pbar <- txtProgressBar(min=1,max=maxg,style=1)
  ig   <- 0
  
  for(j in 1:length(main)){
    
    ttt <- .interactionsFromGibbs(mainx=main[j], bchain=bchain,
                                  specs=snames, xmnames=names(xvector), 
                                  xx=xvector, omitY = omitY, sdScaleX=F, 
                                  standX)
    maine   <- ttt$main
    inter   <- ttt$inter  #
    indirTo <- maine*0
    direct  <- maine + inter  # already standardized for X
    
    if(MEAN){
      dmain  <- colMeans(maine)
      inte   <- colMeans(inter)
      dir    <- colMeans(direct)
    } else {
      dmain  <- apply(maine,2,median)
      inte   <- apply(inter,2,median)
      dir    <- apply(direct,2,median)
    }
    
    mainEffect[sindex,j] <- dmain
    intEffect[sindex,j]  <- inte
    dirEffect[sindex,j]  <- dir
    
    mainSd[sindex,j] <- apply(maine,2,sd)
    intSd[sindex,j]  <- apply(inter,2,sd)
    dirSd[sindex,j]  <- apply(direct,2,sd)
    
    for(g in gs){
      
      if(REDUCT){
        Z  <- matrix(schain[g,],N,r)
        ss <- .expandSigma(sigErrGibbs[g], S, Z = Z, kchain[g,], 
                           REDUCT = T)[sindex,sindex]
        if(sdScaleY)cc <- .cov2Cor(ss)
      } else {
        ss <- .expandSigma(schain[g,], S = S, REDUCT = F)[sindex,sindex]
        if(sdScaleY)cc <- .cov2Cor(ss)
      }
      
      for(s in 1:length(sindex)){
        
        if(REDUCT){
          si <- invWbyRcpp(sigErrGibbs[g], Z[kchain[g,sindex[-s]],])
          if(sdScaleY){
            dc <- diag(sqrt(diag(ss)))[-s,-s]
            ci <- dc%*%si%*%dc
          }
        } else {
          si <- solveRcpp(ss[-s,-s])
          if(sdScaleY)ci <- solveRcpp(cc[-s,-s])
        }
        
        if(!sdScaleY){
          sonk <- ss[drop=F,s,-s]
          e2   <- sonk%*%si%*%direct[g,-s]
        } else {
          sonk <- cc[drop=F,s,-s]
          e2   <- sonk%*%ci%*%direct[g,-s]      # correlation scale
        }
        indirTo[g,s] <- e2
        
        ig <- ig + 1
        setTxtProgressBar(pbar,ig)
        
      } ##############
    }
    
    if(MEAN){
      indirectTo   <- colMeans(indirTo[gs,])
    } else {
      indirectTo   <- apply(indirTo[gs,],2,median)
    }
    indEffectTo[sindex,j]   <- indirectTo
    indSdTo[sindex,j]       <- apply(indirTo[gs,],2,sd)
  } ######################################
  
  if(!is.null(keepNames)){
    wk <- which(rownames(mainEffect) %in% keepNames)
    mainEffect <- mainEffect[wk,]
    intEffect <- intEffect[wk,]
    dirEffect <- dirEffect[wk,]
    indEffectTo   <- indEffectTo[wk,]
    mainSd    <- mainSd[wk,]
    dirSd     <- dirSd[wk,]
    indSdTo   <- indSdTo[wk,]
  }
  
  list(mainEffect = mainEffect, intEffect = intEffect, dirEffect = dirEffect,
       indEffectTo = indEffectTo, mainSd = mainSd, dirSd = dirSd,
       intSd = intSd, indSdTo = indSdTo)
}

.interactionsFromGibbs <- function(mainx,bchain,specs,xmnames=names(xx),
                                   xx=colMeans(xx), omitY=NULL, sdScaleX, 
                                   standX){
  
  # returns main effects and interactions for variable named main
  # xx are values of covariates to condition on
  # mainx is the name of a main effect
  
  if(length(omitY) > 0){
    wob <- numeric(0)
    for(k in 1:length(omitY)){
      wob <- c(wob, grep(omitY[k],colnames(bchain)))
    }
    bchain[,wob] <- 0
    specs <- specs[!specs %in% omitY]
  }
  
  ww   <- grep(':',xmnames)
  int  <- unique( unlist( strsplit(xmnames[ww],':') ) ) 
  int  <- int[int != mainx]
  
  xj <- paste(mainx,specs,sep='_')
  wj <- which(colnames(bchain) %in%  xj)
  if(length(wj) == 0){
    xj <- paste(specs,mainx,sep='_')
    wj <- which(colnames(bchain) %in%  xj)
  }
  
  maine <- bchain[,xj]
  inter <- maine*0
  
  m1 <- paste(mainx,':',sep='')
  m2 <- paste(':',mainx,sep='')
  i1 <- grep( m1,xmnames )
  i2 <- grep( m2,xmnames )
  
  if(sdScaleX)maine <- maine*standX[mainx,'xsd']  #standardize main effect
  
  if( length(i1) > 0 ){
    
    ww <- match(unlist( strsplit(xmnames[i1],m1) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i1)){
      xi <- paste(xmnames[i1[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi <- paste(specs,xmnames[i1[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      xik   <- xx[ox[kk]]
      bik   <- bchain[,xi]
      if(sdScaleX){
        xik <- (xik - standX[ox[kk],'xmean'])/standX[ox[kk],'xsd']
        bik <- bik*standX[mainx,'xsd']*standX[ox[kk],'xsd']
      }
      inter <- inter + bik*xik
    }
  }
  
  if( length(i2) > 0 ){
    
    ww <- match(unlist( strsplit(xmnames[i2],m2) ),xmnames)
    ox <- xmnames[ww[is.finite(ww)]]
    for(kk in 1:length(i2)){
      xi <- paste(xmnames[i2[kk]],specs,sep='_')
      wi <- which(colnames(bchain) %in%  xi)
      if(length(wi) == 0){
        xi    <- paste(specs,xmnames[i2[kk]],sep='_')
        wi <- which(colnames(bchain) %in%  xi)
      }
      xik   <- xx[ox[kk]]
      bik   <- bchain[,xi]
      if(sdScaleX){
        xik <- (xik - standX[ox[kk],'xmean'])/standX[ox[kk],'xsd']
        bik <- bik*standX[mainx,'xsd']*standX[ox[kk],'xsd']
      }
      inter <- inter + bik*xik
    }
  }
  list(main = maine, inter = inter)
}

.stackedBoxPlot <- function( stackList, stackSd=character(0),
                             ylim=NULL,sortBy = NULL, barnames=NULL,
                             col=rep(NULL,length(stackList)),
                             border=rep(NA,length(stackList)),
                             decreasing=T, nsd=1.96, cex=1,
                             legend=NULL, scaleLegend=.1){
  
  # sortBy - if length 1 indicates which variable in stackList to sort by
  #        - if a vector it is the order to plot
  # nds    - no. standard deviations for whiskers
  
  nn  <- length(stackList)
  ord <- c(1:length(stackList[[1]]))
  nx  <- length(ord)
  
  xx <- 0:(nx-1)
  
  if(is.null(ylim)){
    
    ymax <- ymin <- 0
    
    for(j in 1:nn){
      ymax <- ymax + max( c(0,stackList[[j]]),na.rm=T )
      ymin <- ymin + min( c(0,stackList[[j]]),na.rm=T )
    }
    
    ylim <- c(ymin,ymax)
    
    yscale <- diff(ylim,na.rm=T)*.4
    ylim[1] <- ylim[1] - yscale
    ylim[2] <- ylim[2] + yscale
  }
  
  if(!is.null(sortBy)){
    
    if(length(sortBy) > 1){
      ord <- sortBy
    } else {
      ord <- order( stackList[[sortBy]], decreasing = decreasing)
    }
    if(!is.null(barnames))barnames <- barnames[ord]
  }
  
  dy   <- diff(ylim)
  xlim <- c(0,1.2*length(ord))
  
  add <- F
  
  offset <- offsetPos <- offsetNeg <- rep(0,length(stackList[[1]]))
  
  if(is.null(col))col <- c(1:nn)
  
  for(j in 1:nn){
    
    xj <- stackList[[j]][ord]
    names(xj) <- NULL
    
    wp <- which(xj > 0)     # increase
    wn <- which(xj < 0)     # decrease
    
    offset[wp] <- offsetPos[wp]
    offset[wn] <- offsetNeg[wn]
    
    hj <- xj 
    
    barplot(height= hj,offset=offset,xlim=xlim,ylim=ylim,
            col=col[j],border=border[j],add=add)
    
    ww <- grep(names(stackList)[j],names(stackSd))
    if(length(ww) > 0){
      xz <- xx + .5
      xz <- xz*1.2
      
      tall <-  nsd*stackSd[[ww]]
      y1   <-  hj + offset + tall
      y2   <-  hj + offset - tall
      
      for(i in 1:length(ord)){
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=6,col='white')
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=6,col='white')
        
        lines(xz[c(i,i)],c(y1[i],y2[i]),lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y1[c(i,i)],lwd=2,col=col[j])
        lines(c(xz[i]-.1,xz[i]+.1),y2[c(i,i)],lwd=2,col=col[j])
      }
    }
    
    if(j == 1)add <- T
    
    offsetPos[wp] <- offsetPos[wp] + hj[wp]
    offsetNeg[wn] <- offsetNeg[wn] + hj[wn]
    
    if(j == nn & !is.null(barnames)){
      
      xall <- par('usr')[1:2]
      xtic <- c(1:nx)*(diff(xall) - 1)/nx - .8
      
      yy <- offsetPos + .2*dy
      pos <- yy*0 + 1
      wl <- which(abs(offsetNeg) < abs(offsetPos))
      yy[wl] <- offsetNeg[wl] - .2*dy
      pos[wl] <- 4
      text(xtic,yy,barnames,srt=90.,pos=pos,cex=cex)
    }
  } 
  
  if(!is.null(legend)){
    
    dy <- diff(ylim)*scaleLegend
    dx <- 1.2
    x1 <- length(ord)*.02 + 1
    y1 <- ylim[1]
    pos <- 4
    if(legend == 'topright'){
      x1  <- length(ord)
      y1  <- ylim[2]
      dy  <- -dy
      dx <- -dx
      pos <- 2
    }
    if(legend == 'topleft'){
      y1  <- ylim[2]
      dy  <- -dy
    }
    if(legend == 'bottomright'){
      x1  <- length(ord)
      dx <- -dx
      pos <- 2
    }
    for(j in 1:length(stackList)){
      y2 <- y1 + dy
      rect(x1,y1,x1 + 1,y2,col=col[j],border=border[j])
      y1 <- y2
      colj <- col[j]
      if(colj == 'white')colj <- border[j]
      text(x1 + dx,y1 - dy/2,names(stackList)[[j]],col=colj,pos=pos,cex=cex)
    }
  }
  
  invisible( ord )
}  

.getScoreNorm <- function(x,mu,xvar){  #Gneiting/ Raftery proper scoring rule
  
  #outcome x, prediction mean variance (mu, xvar)
  
  - ( (x - mu)^2)/xvar - log(xvar)
  
}

.gjamBaselineHist <- function(y1, bins=NULL, nclass=20){
  
  # add histogram to base of current plot
  
  if(!is.null(bins)){
    hh <- hist(y1,breaks=bins,plot=F)
  } else {
    hh <- hist(y1,nclass=nclass,plot=F)
  }
  
  xvals <- rep(hh$breaks,each=2)
  yvals <- rep(hh$density,each=2)
  
  nb    <- length(hh$breaks)
  yvals <- c( 0, yvals, 0)
  
  rbind(xvals,yvals)
}

.gjamCensorSetup <- function(y,w,z,plo,phi,wm,censorMat){
  
  nc <- ncol(censorMat)
  br <- numeric(0)
  nk <- length(wm)
  n  <- nrow(y)
  
  zk <- y[,wm]*0
  blast <- -Inf
  
  for(j in 1:nc){
    
    valuej <- censorMat[1,j]
    bj     <- censorMat[2:3,j]
    names(bj) <- paste('c-',names(bj),j,sep='')
    
    if(j > 1){
      if(censorMat[2,j] < censorMat[3,j-1] )
        stop('censor intervals must be unique')
      if(bj[1] == br[length(br)])bj <- bj[2]
    }
    br <- c(br,bj)
    nb <- length(br)
    
    zk[ y[,wm] > blast & y[,wm] < bj[1] ] <- nb - 2
    zk[ y[,wm] == valuej ] <- nb - 1
    blast <- br[length(br)]
  }
  
  if(nc == 1){
    zk[zk == 0] <- 2
    br <- c(br,Inf)
  }
  
  zk[zk == 0] <- 1
  br <- matrix(br,nk,length(br),byrow=T)
  
  censk    <- which(y[,wm] %in% censorMat[1,])
  z[,wm]   <- zk
  
  tmp   <- .gjamGetCuts(z,wm)
  cutLo <- tmp$cutLo
  cutHi <- tmp$cutHi
  
  plo[,wm] <- br[cutLo]
  phi[,wm] <- br[cutHi]
  
  ww <- which(plo[,wm,drop=F] == -Inf,arr.ind=T)
  if(length(ww) > 0){
    if(length(wm) == 1){
      mm <- w[,wm]
    }else{
      mm <- apply(w[,wm],2,max)
    }
    plo[,wm][ww] <- -10*mm[ww[,2]]
  }
  
  tmp <-  .tnorm(nk*n,plo[,wm],phi[,wm],w[,wm],1)
  
  w[,wm][censk] <- tmp[censk]
  
  imat <- w*0                    #location in full matrix
  imat[,wm][censk] <- 1
  censValue <- which(imat == 1)
  
  list(w = w, z = z, cutLo = cutLo, cutHi = cutHi, plo = plo, phi = phi,
       censValue = censValue, breakMat = br)
}

.gjamCuts2theta <- function(tg,ss){   # variance to correlation scale
  
  if(length(ss) == 1)return(tg/sqrt(ss))
  nc   <- ncol(tg)
  sr   <- nrow(ss)
  tg/matrix( sqrt(diag(ss)),sr,nc)
}

.gjamGetCuts <- function(zz,wk){
  
  nk <- length(wk)
  n  <- nrow(zz)
  
  cutLo <- cbind( rep(1:nk,each=n), as.vector(zz[,wk]) )
  cutHi <- cbind( rep(1:nk,each=n), as.vector(zz[,wk]) + 1 )
  
  list(cutLo = cutLo, cutHi = cutHi)
}

.gjamGetTypes <- function(typeNames=NULL){
  
  TYPES <- c('CON','PA','CA','DA','CAT','FC','CC','OC')
  FULL  <- c('continuous','presenceAbsence','contAbun','discAbun',
             'categorical','fracComp','countComp','ordinal')
  LABS  <- c('Continuous','Presence-absence','Continuous abundance',
             'Discrete abundance', 'Categorical','Fractional composition',
             'Count composition','Ordinal')
  
  if(is.null(typeNames)){
    names(FULL) <- TYPES
    return( list(typeCols = NULL, TYPES = TYPES, typeFull = FULL, 
                 labels = LABS ) )
  }
  
  S        <- length(typeNames)
  typeCols <- match(typeNames,TYPES)
  
  ww <- which(is.na(typeCols))
  if(length(ww) > 0)stop( paste('type code error',typeNames[ww],sep=' ') )
  
  list(typeCols = typeCols, TYPES = TYPES, typeFull = FULL[typeCols],
       typeNames = typeNames, labels = LABS[typeCols])
}

.gjamHoldoutSetup <- function(holdoutIndex,holdoutN,n){
  
  #holdout samples
  if(length(holdoutIndex) > 0)holdoutN <- length(holdoutIndex)
  if(holdoutN > (n/5))stop('too many holdouts')
  
  inSamples <- c(1:n)
  if(holdoutN > 0){
    if(length(holdoutIndex) == 0)holdoutIndex <- sort( sample(n,holdoutN) )
    inSamples <- inSamples[-holdoutIndex]
  }
  nIn <- length(inSamples)
  
  list(holdoutIndex = holdoutIndex, holdoutN = holdoutN, 
       inSamples = inSamples, nIn = nIn)
}

.gjamMissingValues <- function(x, y, factorList, typeNames){
  
  n <- nrow(x)
  xnames <- colnames(x)
  
  # missing values in x
  xmiss  <- which(!is.finite(x),arr.ind=T)
  nmiss  <- nrow(xmiss)
  missX  <- missX2 <- xprior <- yprior <- numeric(0)
  
  xbound <- apply(x,2,range,na.rm=T)
  
  if(nmiss > 0){         #initialize missing values with means
    xmean    <- colMeans(x,na.rm=T)
    x[xmiss] <- xmean[xmiss[,2]]
    xprior   <- x[xmiss]
    nmiss    <- nrow(xmiss)
    fmiss    <- signif(100*nmiss/length(x[,-1]),2)
    print( paste(nmiss,' values (',fmiss,'%) missing in x imputed'), sep='' )
    missX <- missX2 <- rep(0,nmiss)
  }
  
  # rare y
  tmp  <- gjamTrimY(y,minObs=0,OTHER=F)
  wzo  <- which(tmp$nobs == 0)
  if(length(wzo) > 0){
    stop( ' remove from ydata types never present:',
          paste0(names(wzo),collapse=', '))
  }
  fobs <- tmp$nobs/n
  wlo  <- which(fobs < .01)
  if(length(wlo) > 0){
    flo <- paste0(names(fobs)[wlo],collapse=', ')
    cat(paste('\nPresent in < 1% of obs:',flo,'\n') )
  }
  
  # missing values in y
  ymiss <- which(!is.finite(y),arr.ind=T)
  mmiss <- nrow(ymiss)
  missY <- missY2 <- numeric(0)
  
  if(mmiss > 0){         #initialize missing values with means by TYPEs
    ymean    <- colMeans(y,na.rm=T)
    y[ymiss] <- ymean[ymiss[,2]]
    yprior   <- jitter(y[ymiss])
    fmiss    <- round(100*mmiss/length(y),1)
    mmiss <- nrow(ymiss)
    missY <- missY2 <- rep(0,mmiss)
    print( paste(mmiss,' values (',fmiss,'%) missing in y imputed'), sep='' )
  }
  
  disTypes <- c('DA','CC','OC')
  wdd    <- which(disTypes %in% typeNames)
  if(length(wdd) > 0){
    www <- which( typeNames[ymiss[,2]] %in% disTypes )
    yprior[www] <- floor(yprior[www])
  }
  
  if(nmiss > 0){
    x[xmiss] <- xprior
    
    print(factorList)
    
    if(length(factorList) > 0){
      for(k in 1:length(factorList)){
        wm <- which(xnames[ xmiss[,2] ] == factorList[[k]][1])
        if(length(wm) == 0)next
        wk <- sample(length(factorList[[k]]),length(wm),replace=T)
        xtmp <- x[xmiss[wm,1],factorList[[k]],drop=F]*0
        
        xtmp[ cbind(1:nrow(xtmp),wk) ] <- 1
        x[xmiss[wm,1],factorList[[k]]] <- xtmp
      }
    }
  }
  
  if(mmiss > 0)y[ymiss] <- yprior
  
  list(xmiss = xmiss, xbound = xbound, missX = missX, missX2 = missX2,
       ymiss = ymiss, missY = missY, xprior = xprior, yprior = yprior,
       x = x, y = y)
}

.gjamPlotPars <- function(type='CA',y1,yp,censm=NULL){
  
  if(!is.matrix(y1))y1 <- matrix(y1)
  if(!is.matrix(yp))yp <- matrix(yp)
  
  n       <- nrow(y1)
  nk      <- ncol(y1)
  nbin    <- NULL
  nPerBin <- max( c(10,n*nk/15) )
  breaks  <- NULL
  xlimit  <- range(y1,na.rm=T)
  ylimit  <- range(yp,na.rm=T)
  vlines  <- NULL
  wide    <- NULL
  MEDIAN  <- T
  LOG     <- F
  yss     <- quantile(as.vector(y1),.5, na.rm=T)/mean(y1,na.rm=T)
  
  if(type == 'CA'){
    wpos <- length( which(y1 > 0) )
    nPerBin <-  max( c(10,wpos/15) )
  }
  if(type %in% c('PA', 'CAT')){
    breaks  <- c(-.05,.05,.95,1.05)
    wide    <- rep(.08,4)
    nPerBin <- NULL
    ylimit  <- c(0,1)
    xlimit <- c(-.1,1.1)
  } 
  if(type == 'OC'){
    
    breaks  <- seq(min(y1,na.rm=T)-.5,max(y1,na.rm=T) + .5,by=1)
    wide    <- 1/max(y1)
    nPerBin <- NULL
    ylimit  <- range(yp,na.rm=T)
    xlimit  <- c( min(floor(y1),na.rm=T), max(ceiling(y1),na.rm=T) )
  } 
  if(type == 'DA')MEDIAN <- F
  if(type %in% c('DA','CA')){
    if(yss < .8){
      xlimit <- range(y1,na.rm=T)
      xlimit[2] <- xlimit[2] + 1
      LOG <- T
    }
  }
  if(type %in% c('FC','CC')){
    MEDIAN <- F
    nPerBin <- round( n*nk/50,0 )
  } 
  if(type == 'CC'){
    xlimit[2] <- xlimit[2] + 1
    if(yss <  1){
      LOG <- T
      xlimit[1] <- ylimit[1] <- 1
    }
  }
  
  if( !is.null(censm) ){
    
    cc  <- censm$partition
    vlines  <- numeric(0)
    breaks  <- NULL
    nPerBin <- n*nk/15
    xlimit  <- range(y1,na.rm=T)
    ylimit  <- quantile(yp,c(.01,.99),na.rm=T)
    
    if(ncol(cc) > 1){
      cm     <- unique( as.vector(cc[-1,]) )
      vlines <- cm[is.finite(cm)]
      breaks <- vlines
      nbin   <- nPerBin <- NULL
      uncens <- cbind(cc[3,-ncol(cc)],cc[2,-1])
      wu     <- which( uncens[,1] != uncens[,2] )
      for(m in wu){
        sm <- seq(uncens[m,1],uncens[m,2],length=round(10/length(wu),0))
        if(type == 'DA') sm <- c(uncens[m,1]:uncens[m,2])
        breaks <- c(breaks,sm)
      }
      if(max(cc[1,]) < Inf){
        breaks <- c(breaks, seq(max(breaks),(max(y1,na.rm=T)+1),length=12) )
      } else {
        breaks <- c(breaks,max(y1,na.rm=T) + 1)
      }
      breaks <- sort( unique(breaks) )
    }
  }
  
  if(LOG){
    xlimit[1] <- ylimit[1] <- quantile(y1[y1 > 0],.001, na.rm=T)
    w0     <- which(y1 == 0)
    y1[w0] <- ylimit[1]
    w0     <- which(yp == 0)
    yp[w0] <- ylimit[1]
    #  nPerBin <- nPerBin/2
    ylimit[2] <- max(yp,na.rm=T)
  }
  
  list( y1 = y1, yp = yp, nbin=nbin, nPerBin=nPerBin, vlines=vlines,
        xlimit=xlimit,ylimit=ylimit,breaks=breaks,wide=wide,LOG=LOG,
        POINTS=F,MEDIAN=MEDIAN )
}

.gjamPredictTraits <- function(w,specByTrait,traitTypes){
  
  M  <- nrow(specByTrait)
  tn <- rownames(specByTrait)
  
  ww <- w
  ww[ww < 0] <- 0
  
  tt <- ww%*%t(specByTrait)
  # wf <- grep('FC',traitTypes)
  # if(length(wf) > 0){
  #   w0 <- which(tt[,wf] < 0)
  #   tt[tt[,wf] < 0,wf] <- 0
  #   tsum <- colSums(tt)
  #   tt   <- sweep(tt,1,tsum,'/')
  # }
  tt
}

.initW <- function(tw, x, yy, minw = -ncol(yy), cat=F){
  
  # initialize w for y = 0
  
  X <- x
  X[,-1] <- jitter(X[,-1],factor=1)
  
  XX  <- crossprod(X)
  IXX <- solveRcpp(XX)
  
  for(j in 1:50){
    
    bb <- IXX%*%crossprod(X,tw)
    muw <- X%*%bb
    
    tw[yy == 0] <- muw[yy == 0]    #neg values 
    tw[yy == 0 & tw > 0] <- 0      #no bigger than zero
  }
  tw[tw < minw] <- minw
  # }
  tw
}


.gjamSetup <- function(typeNames, x, y, breakList=NULL, holdoutN, holdoutIndex,
                       censor=NULL, effort=NULL, maxBreaks=100){
  
  Q <- ncol(x)
  n <- nrow(y)
  S <- ncol(y)
  
  effMat <- effort$values
  
  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  typeCols <- tmp$typeCols
  allTypes <- unique(typeCols)
  
  cuts <- cutLo <- cutHi <- numeric(0)
  minOrd <- maxOrd <- breakMat <- numeric(0)
  
  ordShift <- NULL
  
  ordCols  <- which(typeNames == 'OC')
  disCols  <- which(typeNames == 'DA')
  compCols <- which(typeNames == 'CC')
  corCols  <- which(typeNames %in% c('PA','OC','CAT'))
  catCols  <- which(typeNames == c('CAT'))
  
  CCgroups <- attr(typeNames,'CCgroups')
  if(length(CCgroups) == 0)CCgroups <- rep(0,S)
  ngroup <- max(CCgroups)
  
  FCgroups <- attr(typeNames,'FCgroups')
  if(length(FCgroups) == 0)FCgroups <- rep(0,S)
  fgroup <- max(FCgroups)
  
  CATgroups <- attr(typeNames,'CATgroups')
  if(length(CATgroups) == 0)CATgroups <- rep(0,S)
  cgroup <- max(CATgroups)
  
  wo <- grep('others',colnames(y))
  if(length(wo) > 0){
    colnames(y)[wo] <- .replaceString(colnames(y)[wo],'others','other')
  }
  
  other <- grep('other',colnames(y))
  
  colnames(y) <- .cleanNames(colnames(y))
  
  w  <- y 
  if(!is.null(effort))w <- w/effort$values
  
  maxy <- apply(w,2,max,na.rm=T)
  miny <- apply(w,2,min,na.rm=T)
  miny[miny > -maxy] <- -maxy[miny > -maxy]
  maxy[maxy < 0] <- -maxy[maxy < 0]
  maxy <- matrix(maxy, n, S, byrow=T)
  
  z  <- w*0
  z[y == 0] <- 1
  z[y > 0]  <- 2
  plo <- phi <- y*0
  plo[z == 1] <- -2*maxy[z == 1]
  phi[z == 2] <- 2*maxy[z == 2]
  
  censorCON <- censorCA <- censorDA <- numeric(0) # values to be sampled
  sampleW  <- y*0
  sampleW[is.na(sampleW)] <- 1
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    
    if( typeFull[wk[1]] == 'presenceAbsence' ){       
      
      sampleW[,wk] <- 1
      plo[,wk][z[,wk] == 1] <- -10
      phi[,wk][z[,wk] == 2] <- 10
      
      w[,wk] <- .tnorm(nk*n,plo[,wk],phi[,wk],0,1)
      br <- c(-Inf,0,Inf)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('PA',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat     <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'continuous' ){      
      
      sampleW[,wk] <- 0
      z[,wk]   <- 1
      
      if( !is.null(censor) & 'CON' %in% names(censor) ){
        
        wc     <- which(names(censor) == 'CON')
        bc     <- censorCON <- numeric(0)
        
        for(m in wc){
          
          wm     <- censor[[m]]$columns
          cp     <- censor[[m]]$partition
          for(ii in 1:ncol(cp)){
            wmm <- which( y[,wm] == cp[1,ii] | (y[,wm] > cp[2,ii] & y[,wm] < cp[3,ii]) )
            mmm <- cp[2,ii]
            if(mmm == -Inf)mmm <- cp[3,ii] - 10
            censor[[m]]$partition[2,ii] <- mmm
            plo[,wm][wmm] <- mmm
            phi[,wm][wmm] <- cp[3,ii]
          }
          
          tmp    <- .gjamCensorSetup(y,w,z,plo,phi,wm,censorMat=
                                       censor[[m]]$partition)
          z[,wm] <- tmp$z[,wm]
          censorCON <- c(censorCON,tmp$censValue)
          bt       <- tmp$breakMat
          colnames(bt) <- as.character(c(1:ncol(bt)))
          rownames(bt) <- paste('CA',wm,sep='-')
          
          bc <- .appendMatrix(bc,bt,SORT=T,asNumbers=T)
        }
      }
      
      br <- c(-Inf,-Inf,Inf)
      br  <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CON',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'contAbun' ){       
      
      phi[,wk] <- 5*maxy[,wk]
      plo[,wk] <- -phi[,wk]
      
      wy1 <- which(y[,wk] > 0)
      
      w[,wk][wy1] <- y[,wk][wy1]
      
      wy0 <- which(y[,wk] == 0)
      phi[,wk][wy0] <- 0
      
      w[,wk] <- .initW(w[,wk], x, y[,wk], minw = -max(y[,wk],na.rm=T)*5)
      w[,wk][wy1] <- y[,wk][wy1]
      
      br <- c(-Inf,0,Inf)
      br  <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CA',wk,sep='-')
      
      sampleW[,wk][y[,wk] == 0] <- 1
      
      if( !is.null(censor) & 'CA' %in% names(censor) ){
        
        wc     <- which(names(censor) == 'CA')
        bc     <- censorCA <- numeric(0)
        
        for(m in wc){
          
          wm     <- censor[[m]]$columns
          cp     <- censor[[m]]$partition
          for(ii in 1:ncol(cp)){
            wmm <- which(y[,wm] == cp[1,ii] | (y[,wm] > cp[2,ii] & y[,wm] < cp[3,ii]) )
            plo[,wm][wmm] <- cp[2,ii]
            phi[,wm][wmm] <- cp[3,ii]
          }
          
          tmp    <- .gjamCensorSetup(y,w,z,plo,phi,wm,censorMat=
                                       censor[[m]]$partition)
          z[,wm] <- tmp$z[,wm]
          censorCA <- c(censorCA,tmp$censValue)
          bt       <- tmp$breakMat
          colnames(bt) <- as.character(c(1:ncol(bt)))
          rownames(bt) <- paste('CA',wm,sep='-')
          
          bc <- .appendMatrix(bc,bt,SORT=T,asNumbers=T)
        }
        
        mm <- match(rownames(bc),rownames(br))
        
        if(is.na(min(mm)))stop('error in censor list, check for conflicts')
        
        bb <- br[-mm,]
        tmp <- .appendMatrix(bc,bb,SORT=T,asNumbers=T)
        o   <- as.numeric( matrix( unlist(strsplit(rownames(tmp),'-')),
                                   ncol=2,byrow=T)[,2] ) 
        br <- tmp[drop=F,order(o),]
      }
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'discAbun' ){
      
      plo[,wk] <- (y[,wk] - .5)/effMat[,wk]
      phi[,wk] <- (y[,wk] + .5)/effMat[,wk]
      plo[,wk][y[,wk] == 0] <- -5*maxy[,wk][y[,wk] == 0] 
      phi[,wk][y[,wk] == maxy[,wk]] <- 5*maxy[,wk][y[,wk] == maxy[,wk]]
      
      sampleW[,wk] <- 1
      
      disCols <- wk
      
      z[,wk] <- y[,wk] + 1
      w[,wk] <- .tnorm(nk*n,plo[,wk],phi[,wk],w[,wk],1)
      
      n <- nrow(y)
      S <- ncol(y)
      
      br <- c(-Inf,seq(0,(max(y[,wk])-1)),Inf)
      if(length(br) > maxBreaks){
        #     warning('breaks created')
        br <- c(br[1:maxBreaks],Inf)
      }
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('DA',wk,sep='-')
      
      if( !is.null(censor) & 'DA' %in% names(censor) ){
        
        wc     <- which(names(censor) == 'DA')
        bc     <- censorDA <- numeric(0)
        
        for(m in wc){
          wm     <- censor[[m]]$columns
          tmp    <- .gjamCensorSetup(y,w,z,plo,phi,wm,
                                     censorMat=censor[[m]]$partition)
          w[,wm] <- tmp$w[,wm]
          z[,wm] <- tmp$z[,wm]
          plo[,wm] <- tmp$plo[,wm]
          phi[,wm] <- tmp$phi[,wm]
          censorDA <- c(censorDA,tmp$censValue)
          bt       <- tmp$breakMat
          colnames(bt) <- as.character(c(1:ncol(bt)))
          rownames(bt) <- paste('DA',wm,sep='-')
          
          bc <- .appendMatrix(bc,bt,SORT=T,asNumbers=T)
        }
        mm <- match(rownames(bc),rownames(br))
        
        bb <- br[-mm,]
        tmp <- .appendMatrix(bc,bb,SORT=T,asNumbers=T)
        o   <- as.numeric( matrix( unlist(strsplit(rownames(tmp),'-')),
                                   ncol=2,byrow=T)[,2] ) 
        br <- tmp[order(o),]
      }
      
      
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'fracComp' ){
      
      wss <- which(y[,wk] == 0 | y[,wk] == 1)
      sampleW[,wk][wss] <- 1
      
      for(i in 1:fgroup){
        
        if(fgroup == 1){
          wki <- wk
        } else {
          wki <- which(typeCols == k & FCgroups == i)
        }
        nki <- length(wki)
        yki  <- y[,wki]
        
        lo <- plo[,wki]
        hi <- phi[,wki]
        lo[lo < -200/S] <- -200/S
        hi[hi > 3]  <- 3
        plo[,wki] <- lo
        phi[,wki] <- hi
        
        w[,wki] <- .initW(w[,wki], x, yki, minw = -200/S)
      }
      
      br <- c(-1,0,1)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('FC',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] %in% c('countComp','categorical')){
      
      sampleW[,wk] <- 1
      
      ntt <- ngroup
      if(typeFull[wk[1]] == 'categorical')ntt <- cgroup
      
      for(i in 1:ntt){
        
        if(ntt == 1){
          wki <- wk
        } else {
          wki <- which( typeCols == k )
          wki <- wki[ CCgroups[wki] == i | CATgroups[wki] == i ]
        }
        nki <- length(wki)
        yki  <- y[,wki]
        
        if( wki[1] %in% catCols ){
          lo  <- hi <- yki*0
          lo[yki == 0] <- -100
          hi[yki == 0] <- 0
          hi[yki == 1] <- 100
          mu <- yki*0
          mu[lo == 0] <- 20
          mu[hi == 0] <- -20
        } else {
          ee <- rowSums(yki)  + 1
          lo <- (yki - .5)/ee          
          hi <- (yki + .5)/ee
          lo[lo < 0]  <- -20/S
          mu <- yki/ee
        }
        
        z[,wki] <- yki + 1
        
        plo[,wki] <- lo
        phi[,wki] <- hi
        
        tmp <- matrix( .tnorm(nki*n,as.vector(lo),as.vector(hi), as.vector(mu), sig=5),n,nki )
        
        tt <- tmp
        if( !wki[1] %in% catCols ){
          tt[tt < 0] <- 0
          tsum <- rowSums(tt)
          tt   <- sweep(tt,1,tsum,'/')
          tt[tmp < 0] <- tmp[tmp < 0]
        }
        
        #     w[,wki] <- .initW(tt,x,y[,wki], minw = -100, cat=T)
        w[,wki] <- tt
      }
      
      br <- c(-1,0,1)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CC',wk,sep='-')
      rownames(br) <- paste(colnames(y)[wk],rownames(br),sep='_')
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'ordinal' ){
      
      miny <- ordShift <- apply(y[,wk,drop=F],2,min)
      
      y[,wk] <- y[,wk] - matrix(miny,n,nk,byrow=T)   #min value is zero 
      
      nc <- apply(y[,wk,drop=F],2,max)
      
      sampleW[,wk] <- 1
      
      # more than one obs needed in last cell to estimate partition
      ii  <- list(spec = as.vector(matrix(c(1:nk),n,nk,byrow=T)), 
                  ss = as.vector(y[,wk,drop=F]))
      ctmp <- .byIndex(as.vector(y[,wk,drop=F])*0+1,ii,sum)
      
      ncc <- nc + 1
      if(max(ncc) > ncol(ctmp))ncc <- nc
      
      maxOne <- which(ctmp[ cbind(1:nk,ncc) ] == 1)
      
      if(length(maxOne) > 0){
        
        for(m in 1:length(maxOne)){
          mc <- wk[maxOne[m]]
          y[y[,mc] == nc[maxOne[m]],mc] <- nc[maxOne[m]] - 1
        }
        nc <- apply(y[,wk,drop=F],2,max)
      }
      
      ncut <- max(y[,wk,drop=F])
      crow <- c(0:ncut)
      cuts <- t( matrix(crow,(ncut+1),nk) )
      cuts[ cbind((1:nk),nc+1) ] <- Inf
      
      call <- t( apply(cuts,1,cumsum) )
      cuts[call == Inf] <- Inf
      cuts <- cbind(-Inf,cuts)
      if(!is.matrix(cuts))cuts <- matrix(cuts,1)
      
      tmp   <- .gjamGetCuts(y + 1,wk)
      cutLo <- tmp$cutLo
      cutHi <- tmp$cutHi
      
      ss   <- seq(0,(nk-1)*n,by=n)
      wh <- as.vector( outer(holdoutIndex,ss,'+') )
      c1 <- cutLo
      if(length(wh) > 0)c1 <- cutLo[-wh,]
      
      otab <- .byIndex(c1[,1]*0 + 1,INDICES=list('i'=c1[,1],
                                                 'j'=c1[,2]),sum,coerce=T)
      oo <- cbind(0,t( apply(otab,1,cumsum) ))
      wo <- which(oo == 0,arr.ind=T)
      wo[,2] <- as.numeric(colnames(otab))[wo[,2]]
      minOrd <- .byIndex(wo[,2],wo[,1],max)
      
      oo <- cbind(0,t( apply( t(apply(otab,1,rev)),1,cumsum) ))
      wo <- which(oo == 0,arr.ind=T)
      maxOrd <- ncut - .byIndex(wo[,2],wo[,1],max) + 2
      
      plo[,wk] <- cuts[cutLo]
      phi[,wk] <- cuts[cutHi]
      
      z[,wk] <- y[,wk] + 1
      w[,wk] <- matrix( .tnorm(nk*n,plo[,wk],phi[,wk],y[,wk],1),n,nk )
      
      colnames(cuts) <- c(1:ncol(cuts))
      rownames(cuts) <- paste('OC',wk,sep='-')
      rownames(cuts) <- paste(colnames(y)[wk],rownames(cuts),sep='_')
      breakMat <- .appendMatrix(breakMat,cuts,SORT=T,asNumbers=T)
    }
  }
  sord <- .splitNames(rownames(breakMat))$vnam[,1]
  yord <- match(colnames(y),sord)
  breakMat <- breakMat[yord,]
  
  sampleW[censorCON] <- 1
  sampleW[censorCA] <- 1
  sampleW[censorDA] <- 1
  
  wCols <- which(colSums(sampleW) > 0)
  wRows <- which(rowSums(sampleW) > 0)
  
  attr(sampleW,'type') <- 'cols'
  attr(sampleW,'index') <- wCols
  if( sum(sampleW) == 0)attr(sampleW,'type') <- 'none'
  if( sum(sampleW) > 0 & (length(wRows) < length(wCols)) ){
    attr(sampleW,'type') <- 'rows'
    attr(sampleW,'index') <- wRows
  }
  
  ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), 
             discrete_class = as.vector(z))
  classBySpec <- .byIndex(as.vector(z)*0+1,ii,sum)
  rownames(classBySpec) <- colnames(y)
  
  ncc <- min(20,ncol(classBySpec))
  nrr <- min(20,nrow(classBySpec))
  
  list(w = w, z = z, y = y, other = other, cuts = cuts, 
       cutLo = cutLo, cutHi = cutHi, ordShift = ordShift,
       plo = plo, phi = phi, ordCols=ordCols, disCols = disCols, 
       compCols = compCols, corCols = corCols,
       classBySpec = classBySpec, breakMat = breakMat, 
       minOrd = minOrd, maxOrd = maxOrd, sampleW = sampleW,
       censorCA = censorCA, censorDA = censorDA, censorCON = censorCON )
}

.gjamTrueVest <- function(chains,true,typeCode,allTypes,xlim=NULL,ylim=NULL,
                          label=NULL,colors=NULL,add=F,legend=T){
  
  true   <- as.vector(true)
  ntypes <- length(allTypes)
  
  if(is.null(ylim))ylim <- range(chains,na.rm=T)
  if(is.null(xlim))xlim <- range(true,na.rm=T)
  
  if(!is.matrix(chains)){
    chains <- matrix(chains,ncol=1)
    bCoeffTable <- c(mean(chains),sd(chains),quantile(chains,c(.025,.975)),true)
    bCoeffTable <- matrix(bCoeffTable,1)
  } else {
    bCoeffTable <- .processPars(chains,xtrue=true )
  }
  
  if(is.null(colors)){
    colors <- 1
    if(ntypes > 1)colors <- typeCode
  }
  if(length(colors) == 1) colors <- rep(colors,ntypes)
  
  .predVsObs(true,p=chains,xlab='true',xlim=xlim,ylim=ylim,ylab='estimated',
             colors=colors,add=add)
  
  if(ntypes > 1 & legend)legend('topleft',allTypes,text.col=colors,bty='n')
  if(!is.null(label)).plotLabel(label,above=T)
  
  invisible( bCoeffTable )
}

#.gjamUpdateBetaNoPrior <- function(WIX,IXX,sg,...){
#  matrix( .rMVN(1,as.vector(WIX),kronecker(sg,IXX)),nrow(IXX),ncol(WIX) )
#}

.conditionalMVN <- function(xx, mu, sigma, cdex, p=ncol(mu)){  
  # xx, mu are matrices
  # cdex conditional for these variables
  # gdex condition on these variables
  
  if(ncol(xx) != ncol(sigma))stop('ncol(xx) != ncol(sigma)')
  if(ncol(mu) != ncol(sigma))stop('ncol(mu) != ncol(sigma)')
  if(max(cdex) > ncol(mu))stop('max(cdex) > ncol(mu)')
  
  gdex <- (1:p)[-cdex] - 1
  cdex <- cdex - 1
  condMVNRcpp(cdex, gdex, xx, mu, sigma) 
}

.byGJAM <- function(x, i, j, summat=matrix(0,max(i),max(j)), 
                    totmat=summat, fun='mean'){  #
  
  nn <- length(x)
  if( nn != length(i) | nn != length(j) )
    stop('vectors unequal in byFunctionRcpp')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )
    stop('matrix too small')
  
  ww <- which(is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    i <- i[-ww]
    j <- j[-ww]
  }
  
  frommat <- cbind(i,j,x)
  
  nr  <- nrow(frommat)
  
  maxmat <- summat*0 - Inf
  minmat <- summat*0 + Inf
  
  tmp <- byRcpp(nr, frommat, totmat, summat, minmat, maxmat)
  
  if(fun == 'sum')return(tmp$sum)
  if(fun == 'mean'){
    mu <- tmp$sum/tmp$total
    mu[is.na(mu)] <- 0
    return(mu)
  }
  if(fun == 'min'){
    return( tmp$min )
  }
  tmp$max
}

.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(1000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat))){
  
  #lo, hi must be same dimensions as muvec,avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) )  
  r[,whichSample] <- a[,whichSample]
  r
}

.whichFactor <- function(dframe){
  
  if(!is.data.frame(dframe))return(character(0))
  
  tmp <- model.frame(data = dframe)
  ym <- attr( attributes(tmp)$terms, 'dataClasses' )
  
  which(ym == 'factor')
}

.xpredSetup <- function(Y, x, bg, isNonLinX, factorObject, intMat, standMat, standMu,
                        notOther, notStandard){
  
  isFactor   <- factorObject$isFactor
  factorList <- factorObject$factorList
  linFactor  <- numeric(0)
  
  Q      <- ncol(x)
  if(Q == 1){
    return( list(linFactor = linFactor, xpred = x, px = 1, 
                 lox = 1, hix = 1) )
  }
  
  # initialize predicted X
  
  xpred  <- x
  n      <- nrow(x)
  
  xnames <- colnames(x)
  SO     <- length(notOther)
  
  px <- 1:Q
  if(length(isNonLinX) > 0)px <- px[-isNonLinX]
  px <- px[!xnames[px] %in% isFactor]
  px <- px[px != 1]
  
  ii <- grep(':',xnames,fixed=T)
  i2 <- grep('^2',xnames,fixed=T)
  
  qx <- c( 1, ii, i2)
  qx <- c(1:Q)[-qx]
  bx <- bg[drop=F,qx,notOther]
  cx <- crossprod(t(bx))
  if(length(cx) == 1){
    cx <- 1/(cx*1.01)
  } else {
    diag(cx) <- .0000001 + diag(cx) 
    cx <- solve(cx)
  }
  
  xx <- (Y[,notOther] - matrix(bg[1,notOther],n,SO,byrow=T))%*%t(bx)%*%cx
  colnames(xx) <- xnames[qx]
  scol      <- colnames(xx)[!colnames(xx) %in% notStandard]
  xx[,scol] <- sweep(xx[,scol,drop=F],2,colMeans(xx[,scol,drop=F]),'-')
  xx[,scol] <- sweep(xx[,scol,drop=F],2,apply(xx[,scol,drop=F],2,sd),'/')
  xpred[,qx] <- xx
  
  xpred[xpred < -3] <- -3
  xpred[xpred > 3] <- 3
  xpred[!is.finite(xpred)] <- 0
  
  if(length(intMat) > 0){
    for(k in 1:nrow(intMat)){
      xpred[,intMat[k,1]] <- xpred[,intMat[k,2]]*xpred[,intMat[k,3]]
    }
  }
  
  
  if(length(isFactor) > 0){
    xpred[,isFactor] <- 0
    
    for(k in 1:length(factorList)){
      kf  <- lf <- factorList[[k]]
      
      if( !is.null(isNonLinX) ){
        xin <- xnames[isNonLinX]
        lf  <- kf[!kf %in% xin]
      }
      if(length(lf) == 0)next
      lf  <- match(lf,xnames)
      ww  <- which(is.finite(lf))
      
      wt  <- colSums(x[,c(1,lf)])   #random, but weighted by prevalence
      wt  <- wt/sum(wt)
      sk  <- sample(c(1,lf),n, replace=T, prob=wt)
      xpred[ cbind(c(1:n),sk) ] <- 1
      
      if(length(ww) == 0)next
      lf <- c(1,lf)   # intercept is reference
      linFactor <- append(linFactor, list(lf))
    }
  }
  
  lox       <- apply(x,2 ,min)
  hix       <- apply(x,2,max)
  
  lox[isFactor] <- -3
  hix[isFactor] <- 3
  if(length(intMat) > 0){
    lox[intMat[,1]] <- -3
    hix[intMat[,1]] <- 3
  }
  
  ws        <- which(notStandard %in% xnames)
  if(length(ws) == 0){
    notStandard <- NULL
  } else {
    notStandard <- notStandard[ws]
    lox[notStandard] <- standMu[notStandard,1] - 3*standMat[notStandard,1]
    hix[notStandard] <- standMu[notStandard,1] + 3*standMat[notStandard,1]
  }
  
  list(linFactor = linFactor, xpred = xpred, px = px, lox = lox, hix = hix)
}

.blockDiag <- function(mat1,mat2){
  
  #creates block diagional
  
  if(length(mat1) == 0)return(mat2)
  
  namesc <- c(colnames(mat1),colnames(mat2))
  namesr <- c(rownames(mat1),rownames(mat2))
  
  nr1 <- nrow(mat1)
  nr2 <- nrow(mat2)
  nc1 <- ncol(mat1)
  nc2 <- ncol(mat2)
  nr  <- nr1 + nr2
  nc  <- nc1 + nc2
  
  new <- matrix(0,nr,nc)
  new[ 1:nr1, 1:nc1 ] <- mat1
  new[ (nr1+1):nr, (nc1+1):nc ] <- mat2
  colnames(new) <- namesc
  rownames(new) <- namesr
  new
}

.getContrasts <- function(facK, fnames){
  
  # D - x to z
  # L - beta to alpha
  # facK - name of factor
  # fnames - character of factor levels
  
  ff <- paste(facK,fnames,sep='')
  
  Q  <- length(fnames)
  cc <- diag(Q)
  cc[1,] <- -1
  dd <- cc
  dd[1] <- 1
  cc[,1] <- 1
  colnames(cc) <- colnames(dd) <- c('intercept',ff[-1])
  rownames(cc) <- rownames(dd) <- fnames
  L <- t( solve(cc) )
  rownames(cc) <- rownames(L) <- rownames(dd) <- ff
  list(C = cc, D = dd, L = L)
}

.getUnstandX <- function(xx, standRows, xmu, xsd, intMat){
  # design to unstandard scale
  
  xUnstand <- xx
  xUnstand[,standRows] <- t( xmu[standRows] + 
                               t(xx[,standRows,drop=F])*xsd[standRows] )
  
  if(length(intMat) > 0){
    for(j in 1:nrow(intMat)){
      xUnstand[,intMat[j,1]] <- xUnstand[,intMat[j,2]] * xUnstand[,intMat[j,3]] 
    }
  }
  S2U <- ginv(crossprod(xUnstand))%*%t(xUnstand)  # (X'X){-1}X'
  rownames(S2U) <- colnames(xx)
  list(xu = xUnstand, S2U = S2U)
}

.getStandX <- function(xx, standRows, xmu=NULL, xsd=NULL, intMat=NULL){
  
  xstand <- xx
  if(is.null(xmu))xmu <- colMeans(xx[,standRows,drop=F],na.rm=T)
  if(is.null(xsd))xsd <- apply(xx[,standRows,drop=F],2,sd,na.rm=T)
  
  xstand[,standRows] <- t( (t(xx[,standRows]) - xmu)/xsd )
  
  if(length(intMat) > 0){
    for(j in 1:nrow(intMat)){
      xstand[,intMat[j,1]] <- xstand[,intMat[j,2]] * xstand[,intMat[j,3]] 
    }
  }
  list(xstand = xstand, xmu = xmu, xsd = xsd)
}

.getHoldLoHi <- function(yh, wh, pl, ph, eff, ymax, typeNames, cutg, ordCols){
  
  # update plo, phi for holdouts, yh is prediction
  
  allTypes <- unique(typeNames)
  
  for(k in 1:length(allTypes)){
    
    tk <- allTypes[k]
    wk <- which(typeNames == tk)
    
    if(tk == 'CON')next
    
    if(tk == 'PA'){
      pl[,wk][yh[,wk] == 0] <- -10
      pl[,wk][yh[,wk] == 1] <- 0
      ph[,wk][yh[,wk] == 0] <- 0
      ph[,wk][yh[,wk] == 1] <- 10
    }
    if(tk == 'CA'){
      ym <- max(ymax[wk])
      pl[,wk][yh[,wk] == 0] <- -5*ym
      pl[,wk][yh[,wk] > 0]  <- 0
      ph[,wk][yh[,wk] == 0] <- 0
      ph[,wk][yh[,wk] > 0] <- 5*ym
    }
    if(tk == 'DA'){
      ym <- max(ymax[wk])
      ee <- 1
      if(!is.null(eff))ee <- eff[,wk]
      pl[,wk] <- (yh[,wk] - .5)/ee
      ph[,wk] <- (yh[,wk] + .5)/ee
      pl[,wk][yh[,wk] == 0] <- -5*ym
      pl[,wk][yh[,wk] == ym] <- 5*ym
    }
    if(tk == 'FC'){
      pl[,wk][yh[,wk] == 0] <- -5
      pl[,wk][yh[,wk] > 0]  <- 0
      pl[,wk][yh[,wk] > 1]  <- 1
      ph[,wk][yh[,wk] == 0] <- 0
      ph[,wk][yh[,wk] > 0]  <- 1
      ph[,wk][yh[,wk] == 1] <- 5
    }
    if(tk == 'CC'){
      ym <- rowSums(yh[,wk,drop=F])
      ee <- matrix(ym,nrow(yh),length(wk))
      pl[,wk] <- (yh[,wk] - .5)/ee
      ph[,wk] <- (yh[,wk] + .5)/ee
      pl[,wk][yh[,wk] == 0] <- -5
      pl[,wk][yh[,wk] == ee] <- 5
    }
  }
  
  list(pl = pl, ph = ph)
}

#.setupReduct <- function(modelList, S, Q, n){
#  REDUCT <- F
#  N <- r <- NULL
#  rl <- NULL
#  if( 'REDUCT' %in% names(modelList) ){
#    rl <- list(N = NULL, r = NULL )
#    if(!modelList$REDUCT)return( rl )
#  }
# npar <- (S+1)/2 + Q
#  
##  ratio <- 1/5
#  N <- min(c(5, S))
#  r <- N - 1
#  if(npar/n > ratio){
#    N <- ceiling( ( ratio*n - Q )/5 )
#    if(N > 25)N <- 25
#    if(N < 4)N  <- 4
#    r <- ceiling( N/2 )
#  }
#  if( 'reductList' %in% names(modelList) ){
#   REDUCT <- T
#    rl <- modelList$reductList
#    N <- rl$N
#    r <- rl$r
#    if(N >= S){
#      N <- S - 1
#      warning(' dimension reduction requires reductList$N < no. responses ')
#    }
#  }
#  if( !'reductList' %in% names(modelList) ){
#    rl <- list(r = r, N = N, alpha.DP = S)
#  }
#  rl
#}

##change J.Clark 

.setupReduct <- function(modelList, S, Q, n){
 if((is.null(modelList$reductList$DRtype)) || ((modelList$reductList$DRtype=="basic"))){
    N <- r <- rl <- NULL
    if( 'REDUCT' %in% names(modelList) | 'reductList' %in% names(modelList) ){
      rl <- list(N = NULL, r = NULL )
      if(('REDUCT' %in% names(modelList))&&!modelList$REDUCT)return( rl ) ##REDUCT = F in modelList overrides automatic dimension reduction.
    #}
    #automatic mode for dimension reduction
    if(n < 2*S | S > 200){
      N  <- round( S/3 )
      if(N > 25)N <- 25
      if(N <= 4)N <- 4
      r  <- ceiling( N/2 )
    } 
     else{
       N <- modelList$reductList$N
       r <- modelList$reductList$r
     }
      rl <- list(r = r, N = N, alpha.DP = S)
      warning( 'dimension reduction' )
    }

 }else{
   if(modelList$reductList$DRtype %in% c("1","2","3")){
     rl<- modelList$reductList
   } else stop("Incorrectly specified DRtype")
 }
  rl
}

##Change for PY

.getTimeIndex <- function(timeList, other, notOther, xdata, x, xl, y, w ){
  
  Q <- ncol(x)
  n <- nrow(x)
  xnames <- colnames(x)
  snames <- colnames(y)
  
  times    <- timeList$times
  if(is.null(times))
    stop(' column name "times" needed for time-series model' )
  timeZero <- which(xdata[,times] == 0)
  if(length(timeZero) == 0)stop(' must have time zero in xdata[,time] ')
  timeLast <- timeZero - 1
  timeLast <- timeLast[-1]
  timeLast <- c(timeLast,nrow(xdata))
  
  ix <- 1:n
  t1 <- ix[-timeZero]
  t0 <- t1 - 1
  t2 <- t1 + 1
  tindex <- cbind(t0,t1,t2)
  S  <- ncol(y)
  
  tindex <- tindex[!tindex[,'t1'] %in% timeLast,]
  
  i1 <- seq(1,nrow(tindex),by=2)
  i2 <- seq(2,nrow(tindex),by=2)
  i1 <- i1[i1 < max(i2)]
  
  maxTime <- max(xdata$times)
  inSamples <- tindex[,2]
  
  # beta
  loBeta <- hiBeta <- NULL
  if('betaPrior' %in% names(timeList)){
    loBeta <- timeList$betaPrior$lo
    hiBeta <- timeList$betaPrior$hi
    beta   <- (loBeta + hiBeta)/2
    beta[is.na(beta)] <- 0
  } else{
    beta <- matrix(0,Q,S)
    rownames(beta) <- colnames(x)
    BPRIOR <- F
  }
  
  tmp <- .betaPrior(beta, notOther, loBeta, hiBeta)
  bg <- tmp$beta; loB <- tmp$loB; hiB <- tmp$hiB
  wB <- tmp$wB; BPRIOR <- tmp$BPRIOR
  bg[is.nan(bg)] <- 0
  
  tmp <- .getPattern(bg[,notOther], wB)
  Brows <- tmp$rows
  Bpattern <- tmp$pattern
  bg[!is.finite(bg)] <- 0
  
  # alpha
  
  alphaPrior <- NULL
  
  if( 'alphaPrior' %in% names(timeList) ){
    loAlpha <- timeList$alphaPrior$lo
    hiAlpha <- timeList$alphaPrior$hi
  } else{
    alpha <- diag(NA,S)
    diag(alpha) <- -1
  }
  
  tmp <- .alphaPrior(w, tindex, timeList$alphaPrior)
  Amat <- tmp$Amat; loAmat <- tmp$loAmat; hiAmat <- tmp$hiAmat
  wA <- tmp$wA; Umat <- tmp$Umat; umat <- tmp$umat
  uindex <- tmp$uindex; aindex <- tmp$aindex
  
  U <- nrow(Amat)
  Umat <- matrix(0,n,U)
  wz   <- w
  wz[wz < 0] <- 0
  
  Umat <- wz[,uindex[,1]]*wz[,uindex[,2]] 
  
  tmp   <- .getPattern(loAmat, wA)
  Arows <- tmp$rows
  Apattern <- tmp$pattern
  Amat[!is.finite(Amat)] <- 0
  
  # lambda
  if('lambdaPrior' %in% names(timeList)){
    lprior <- timeList$lambdaPrior
  } else{
    lprior <- timeList$betaPrior
  }
  tmp <- .lambdaPrior(lprior, w, xl, tindex, xnames, 
                      snames, other, notOther)
  Lmat <- tmp$Lmat; loLmat <- tmp$loLmat; hiLmat <- tmp$hiLmat
  wL <- tmp$wL; gindex <- tmp$gindex; Vmat <- tmp$Vmat
  
  ltmp <- matrix(NA,nrow(Lmat),length(notOther))
  ltmp[wL] <- 1
  
  tmp <- .getPattern(ltmp, wL)
  Lrows <- tmp$rows
  Lpattern <- tmp$pattern
  Lmat[!is.finite(Lmat)] <- 0
  
  list(Lmat = Lmat, Lpattern = Lpattern, wL = wL, gindex = gindex,
       Vmat = Vmat, Lrows = Lrows, loLmat = loLmat, hiLmat = hiLmat,
       Arows = Arows, Amat = Amat, Apattern = Apattern, wA = wA, 
       Umat = Umat, uindex = uindex,loAmat = loAmat, hiAmat = hiAmat,
       aindex = aindex, Brows = Brows, bg = bg, Bpattern = Bpattern, wB = wB, 
       loB = loB, hiB = hiB, timeZero = timeZero, 
       timeLast = timeLast, maxTime = maxTime, inSamples = inSamples, 
       tindex = tindex[,1:2], i1 = i1, i2 = i2)
}

.checkYfactor <- function(ydata, typeNames){
  
  yordNames <- NULL
  
  wf <- which( sapply(ydata,is.factor) )
  
  if(length(wf) > 0){
    
    if(!all(typeNames[wf] == 'CAT'))
      stop('factors in ydata must be CAT data')
    
  }
  return( list(ydata = ydata, yordNames = yordNames) )
  
  # disabled:
  
  yordNames <- vector('list', length=length(wf))
  names(yordNames) <- names(ydata)[wf]
  if(all(!is.ordered(ydata[[j]])))
    warning('OC responses as factors must be ordered')
  
  for(j in wf){
    
    jlev <- attr(ydata[[j]],'levels')
    if('NA' %in% jlev)jlev <- jlev[jlev != 'NA']
    yordNames[[j]] <- jlev
    yj <- as.numeric(ydata[[j]])
    yj[ydata[[j]] == 'NA'] <- NA
    ydata[,j] <- yj
  }
  
  list(ydata = ydata, yordNames = yordNames)
}

.buildEffort <- function(y, effort, typeNames){
  
  S <- length(typeNames)
  
  effMat <- y*0 + 1
  effMat[is.na(effMat)] <- 1
  
  if( is.null(effort) ){
    effort <- list(columns = 1:S, values = effMat)
  } else {
    effMat[,effort$columns] <- effort$values
    effort$values <- effMat
    if(!is.null(colnames(effort$values)))colnames(effMat) <- .cleanNames(colnames(effMat))
  }
  effort$columns <- 1:S
  
  we <- which(effort$values == 0 | is.na(effort$values))
  if(length(we) > 0){
    effort$values[we] <- 1
    if( any(c('DA','CC') %in% typeNames) )
      warning('missing or zero values in effort')
  }
  effMat[,!typeNames == 'DA'] <- 1
  effMat[effMat == 0] <- 1
  effort <- list(columns = effort$columns, values = effMat)
  
  effort
}

.setupFactors <- function(xdata, xnames, factorObject){
  
  factorList <- factorObject$factorList
  contrast   <- factorObject$contrast
  
  Q  <- length(xnames)
  if(Q == 1){
    return( list(dCont = matrix(1)) )
  }
  q1 <- Q - 1
  fnames <- xnames
  findex <- character(0)
  nfact  <- length(factorList)
  
  if(nfact > 0){              # exclude main effects of factors
    findex <- sort( unique( unlist(factorList) ) )
    fnames <- fnames[!fnames %in% findex]
  }
  
  tmp <- diag(length(fnames))
  rownames(tmp) <- colnames(tmp) <- fnames
  if(length(tmp) < 2){
    eCont <- frow <- intercept <- numeric(0)
  } else {
    eCont <- tmp[drop=F,-1,]
    frow  <- rep(0,nrow(eCont))
    intercept <- rep(0,nrow(eCont))
  }
  dCont <- lCont <- eCont
  
  if(nfact > 0){
    
    for(k in 1:nfact){
      
      cm <- contrast[[k]]
      colnames(cm) <- factorList[[k]]
      rownames(cm) <- paste(names(factorList)[[k]],rownames(cm),sep='')
      
      facK <- names(factorList)[[k]]
      
      wx <- match(facK,colnames(xdata))
      
      fnames <- as.character( levels(xdata[[wx]]) ) 
      mm     <- .getContrasts(facK, fnames)
      D  <- mm$D                      # for Z <- x%*%D; 
      L  <- mm$L                      # for A <- L%*%bg; 
      C  <- mm$C                      # L <- solve(t(C)); C = solve(t(L))
      
      if(length(eCont) > 1){
        eCont <- .blockDiag(eCont,cm)
        dCont <- .blockDiag(dCont,D[,-1,drop=F])
        lCont <- .blockDiag(lCont,L[,-1,drop=F])
        ec    <- nrow(lCont)
        bc    <- ec - nrow(L) + 1
        lCont[bc:ec,1] <- L[,1]
        dCont[bc,1] <- -1
      } else {
        eCont <- cbind(0,cm)
        colnames(eCont)[1] <- 'intercept'
        dCont <- D
        lCont <- L
      }
      nr2   <- nrow(cm)
      nc2   <- ncol(cm)
      intercept <- c(intercept,rep(1,nr2))
      
      frow <- c(frow,rep(k,nr2))
    }
    
    eCont[,1] <- intercept
  }
  
  eCont <- eCont[drop=F,,xnames]
  
  dCont <- t(dCont[drop=F,,xnames])
  dCont[1,] <- abs(dCont[1,])
  lCont <- lCont[drop=F,,xnames]
  
  q1 <- nrow(eCont)                # level names only
  fnames   <- rownames(eCont)
  facList2 <- factorList
  if(nfact > 0){
    for(j in 1:nfact){
      wj <- which(names(xdata) == names(factorList)[j])
      facList2[[j]] <- levels(xdata[[wj]])
    }
  }
  
  fmat <- matrix(0,q1,q1)
  colnames(fmat) <- rownames(fmat) <- fnames
  
  findex <- match(findex,xnames)
  
  list(factorList = factorList, facList2 = facList2, fmat = fmat, fnames = fnames,
       q1 = q1, lCont = lCont, dCont = dCont, eCont = eCont, findex = findex)
}



gjamSensitivity <- function(output, group=NULL, nsim=100){
  
  REDUCT <- F
  
  standRows   <- output$inputs$standRows
  factorBeta  <- output$inputs$factorBeta
  notOther    <- output$inputs$notOther
  standMat    <- output$inputs$standMat
  notStandard <- output$modelList$notStandard
  ng          <- output$modelList$ng
  burnin      <- output$modelList$burnin
  x           <- output$inputs$x
  y           <- output$inputs$y
  beta        <- output$parameters$betaMu
  snames      <- colnames(y)
  xnames      <- colnames(x)
  Q <- length(xnames)
  S <-  length(snames)
  S1 <- length(notOther)
  
  bgibbs    <- output$chains$bgibbs
  sgibbs    <- output$chains$sgibbs
  if('kgibbs' %in% names(output$chains)){
    REDUCT <- T
    kgibbs      <- output$chains$kgibbs   
    sigErrGibbs <- output$chains$sigErrGibbs
    N <- output$modelList$reductList$N
    r <- output$modelList$reductList$r
  }
  
  jj <- sample(burnin:ng,nsim,replace=T)
  i  <- 1
  
  for(j in jj){
    
    bg <-  matrix(bgibbs[j,],Q,S)
    rownames(bg) <- xnames
    colnames(bg) <- snames
    
    if(!REDUCT){
      sg <- .expandSigma(sgibbs[j,], S = S, REDUCT = F) 
      si <- solveRcpp( sg ) 
      
    } else {
      Z  <- matrix(sgibbs[j,],N,r)
      sg <- .expandSigma(sigErrGibbs[j], S, Z = Z, kgibbs[j,], REDUCT = T)
      si <- invWbyRcpp(sigErrGibbs[j], Z[kgibbs[j,],])
    }
    
    tmp <- .contrastCoeff(beta=bg[,notOther], 
                          notStand = notStandard[notStandard %in% xnames], 
                          sigma = sg[notOther,notOther],
                          sinv = si[notOther,notOther],
                          stand = standMat, factorObject=factorBeta,
                          conditional = group)
    if(i == 1){
      fmat <- matrix(0,nsim,ncol(tmp$sens))
    }
    
    fmat[i,] <- diag(tmp$sens)
    i <- i + 1
  }
  colnames(fmat) <- colnames(tmp$sens)
  fmat
}

.factorCoeffs2Zero <- function(factorObject, snames, priorObject){
  
  zero  <- numeric(0)
  
  for(k in 1:factorObject$nfact){
    
    wk <- grep('_',factorObject$missFacSpec[[k]])
    
    if(length(wk) > 0){
      sx <- .splitNames(factorObject$missFacSpec[[k]])$vnam
      ij <- cbind(match(sx[,2],rownames(priorObject$lo)),match(sx[,1],snames))
      zero <- rbind(zero,ij)
    }
  }
  zero
}

.gjam <- function(formula, xdata, ydata, modelList){
  
  holdoutN      <-  0
  holdoutIndex  <- numeric(0)
  modelSummary  <- betaPrior  <- traitList <- effort <- NULL
  specByTrait   <- traitTypes <- breakList <- notStandard <- NULL
  censor <- censorCA <- censorDA <- CCgroups <- FCgroups <- intMat <- NULL
  reductList <- y0 <- N  <- r <- otherpar <- pg <- NULL
  ng     <- 2000
  burnin <- 500
  REDUCT <- TRAITS <- FULL <- F
  PREDICTX <- T
  lambdaPrior <- betaPrior <- NULL
  
  RANDOM <- F              # random group intercepts
  
  TIME <- F
  timeList <- timeZero <- timeLast <- timeIndex <- groupIndex <- 
    rowInserts <- Lmat <- Amat <- beta <- NULL
  
  ematAlpha <- .5
  
  
  #alpha.DP <- ncol(ydata)          # not needed
  
  
  #if(alpha.DP == 1) #no more correct now
  if(ncol(ydata) == 1)
    stop('multivariate model: at least 2 columns needed in ydata')
  
  for(k in 1:length(modelList))assign( names(modelList)[k], modelList[[k]] )
  
  if('CCgroups' %in% names(modelList))attr(typeNames,'CCgroups')  <- CCgroups
  if('FCgroups' %in% names(modelList))attr(typeNames,'FCgroups')  <- FCgroups
  if('CATgroups' %in% names(modelList))attr(typeNames,'CATgroups') <- CATgroups
  
  if(!is.null(timeList)){
    if("betaPrior" %in% names(timeList)){
      colnames(timeList$betaPrior$lo) <- 
        colnames(timeList$betaPrior$hi) <- 
        .cleanNames(colnames(timeList$betaPrior$lo))
    }
    if("lambdaPrior" %in% names(timeList)){
      colnames(timeList$lambdaPrior$lo) <- colnames(timeList$lambdaPrior$hi) <- 
        .cleanNames(colnames(timeList$lambdaPrior$lo))
    }
    
    for(k in 1:length(timeList))assign( names(timeList)[k], timeList[[k]] )
    TIME <- T
    REDUCT <- T
    BPRIOR <- T
    holdoutN      <-  0
    holdoutIndex  <- numeric(0)
  }
  
  if(!is.null(traitList)){
    TRAITS <- T
    for(k in 1:length(traitList))assign( names(traitList)[k], traitList[[k]] )
    
    stt <- .replaceString(colnames(specByTrait),'_','')
    colnames(specByTrait) <- stt
    colnames(plotByTrait) <- stt
    colnames(traitList$specByTrait) <- stt
    colnames(traitList$plotByTrait) <- stt
    modelList$traitList <- traitList
  }
  
  if(burnin >= ng) stop( 'burnin must be < no. MCMC steps, ng' )
  if('censor' %in% names(modelList)){
    for(k in 1:length(censor)){
      if( nrow(censor[[k]]$partition) != 3 )
        stop('censor matrix: 3 rows for value, lo, hi')
      rownames(censor[[k]]$partition) <- c('value','lo','hi')
    }
  }
  
  if(missing(xdata)) xdata <- environment(formula)
  
  S <- ncol(ydata)
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  if(length(typeNames) != S) 
    stop('typeNames must be one value or no. columns in y')
  
  ############### factors in y
  
  tmp <- .checkYfactor(ydata, typeNames)
  ydata <- tmp$ydata; yordNames <- tmp$yordNames
  
  if(TRAITS){
    if(!all( typeNames %in% c('CC','FC') ) )
      stop('trait prediction requires composition data (CC or FC)')
    if(nrow(plotByTrait) != nrow(ydata))
      stop('nrow(plotByTrait) must equal nrow(ydata)')
    if(ncol(plotByTrait) != length(traitTypes))
      stop('ncol(plotByTrait) must equal length(traitTypes)')
    if(ncol(plotByTrait) != length(traitTypes))
      stop('ncol(plotByTrait) must equal length(traitTypes)')
    ii <- identical(rownames(specByTrait),colnames(ydata))
    if(!ii){
      ww <- match(colnames(ydata),rownames(specByTrait) )
      if( is.finite(min(ww)) ){
        specByTrait <- specByTrait[ww,]
      } else {
        stop( 'rownames(specByTrait) must match colnames(ydata)' )
      }
    }
    if(typeNames[1] == 'CC'){
      ytmp <- round(ydata,0)
      ytmp[ytmp == 0 & ydata > 0] <- 1
      ydata <- ytmp
      rm(ytmp)
    }
  }
  
  tmp <- .buildYdata(ydata, typeNames)
  y   <- tmp$y
  ydataNames <- tmp$ydataNames
  typeNames  <- tmp$typeNames
  CCgroups   <- tmp$CCgroups
  FCgroups   <- tmp$FCgroups
  CATgroups  <- tmp$CATgroups
  if(TRAITS) rownames(specByTrait) <- colnames(y)
  
  S <- ncol(y)
  n <- nrow(y)
  
  cat("\nObservations and responses:\n")
  print(c(n, S))
  
  tmp    <- .buildEffort(y, effort, typeNames)
  effort <- tmp
  effMat <- effort$values
  modelList$effort <- effort
  re <- floor( diff( range(log10(effMat),na.rm=T) ) )
  if(re > 2)
    message(paste('sample effort > ', re, ' orders of magnitude--consider units near 1',sep='') )
  
  
  tmp      <- .gjamGetTypes(typeNames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  
  tmp <- .gjamXY(formula, xdata, y, typeNames, notStandard)
  x      <- tmp$x; y <- tmp$y; snames <- tmp$snames
  xdata  <- tmp$xdata; xnames <- tmp$xnames
  interBeta   <- tmp$interaction 
  factorBeta  <- tmp$factorAll
  designTable <- tmp$designTable;    xscale <- tmp$xscale
  predXcols   <- tmp$predXcols
  standMat    <- tmp$standMat;      standMu <- tmp$standMu  
  standRows   <- tmp$standRows;    
  xdataNames  <- tmp$xdataNames
  notStandard <- tmp$notStandard[tmp$notStandard %in% xnames]
  
  factorLambda <- interLambda <- NULL
  
  if(!is.null(lambdaPrior)){
    
    lformula <- attr(lambdaPrior$lo,'formula')
    
    tmp <- .gjamXY(lformula, xdata, y, typeNames, notStandard)
    xl   <- tmp$x
    mm <- match(colnames(xl),colnames(xdata))
    wm <- which(is.finite(mm))
    if(length(wm) > 0){
      xdata[,mm[wm]] <- xl[,wm]
    }
    
    xlnames <- tmp$xnames
    interLambda   <- tmp$interaction
    factorLambda <- tmp$factorAll
    
    designTable <- list(beta = designTable, lambda = tmp$designTable)
    
    standMatL    <- tmp$standMat;      standMuL <- tmp$standMu
    standRowsL    <- tmp$standRows;    
    notStandardL <- tmp$notStandard[tmp$notStandard %in% xlnames]
  }
  
  modelList     <- append(modelList, list('formula' = formula,
                                          'notStandard' = notStandard))
  
  Q <- ncol(x)
  
  tmp <- .gjamMissingValues(x, y, factorBeta$factorList, typeNames)
  xmiss  <- tmp$xmiss;   xbound <- tmp$xbound; 
  ymiss  <- tmp$ymiss;   missY <- tmp$missY
  xprior <- tmp$xprior;  yprior <- tmp$yprior
  nmiss  <- nrow(xmiss); mmiss  <- nrow(ymiss)
  x  <- tmp$x; y <- tmp$y
  
  if(TIME){
    tmp <- .gjamMissingValues(xl, y, factorLambda$factorList, typeNames)
    xlmiss  <- tmp$xmiss;   xlbound <- tmp$xbound; 
    xlprior <- tmp$xprior
    nlmiss  <- nrow(xmiss)
    xl <- tmp$x
  }
  reductList <- .setupReduct(modelList, S, Q, n) ##########
  N <- reductList$N; r <- reductList$r ; K_pr <-  reductList$K; 
  PY_var <- reductList$V
  DRtype <- reductList$DRtype
  if(is.null(DRtype)){DRtype<-"basic"
  alpha.DP <- S
  }else{
    if(!(DRtype %in% c("1","2","3"))){stop("The type of dimension reduction is not valid")} 
  }
  if((!is.null(PY_var)&&!DRtype %in% c("3"))){
      stop("Variance specified only for fixed PY")
  }
  ## change for prior K(for "basic" no K assumed)
  if(is.null(K_pr)&(DRtype %in% c("1","2","3"))){stop("Prior number of groups not specified, if you don't need the prior information choose basic version")} 
  ## Change parameters for Ga(shape, rate)
  if(DRtype=="1") { 
    gamma_pars<- compute_gamma_parameters(fun=function(x) simulatuion_function_DPM(x,funct=functionDPM,ns=30000,Sn=S,N_tr = N), K=K_pr)
    rate <- gamma_pars$nu2
    shape <-gamma_pars$nu1
    alpha.DP <- 1
    cat(c(rate,shape),"\n rate and shape \n") 
  }
  if(DRtype=="2") { 
    discount.PY<-reductList$sigma_py; 
    alpha.PY<-reductList$alpha_py; 
    N<- reductList$N
    Precomp_matrix <- reductList$Precomp_mat
    cat(c(alpha.PY,discount.PY),"\n alpha and sigma \n")
    ptr_logv_comp_mat <- create_xptr("log_v_pdf_comp_mat")
  }
  if(DRtype=="3") { 
    if(!(is.null(PY_var))) {
      py_params <-  compute_fixed_parameters_PY_2d(K_pr,PY_var,S)
      sigma_py<-py_params$sigma
      alpha.DP<-py_params$alpha
    } else{
      discount.PY<-reductList$sigma_py 
   # alpha.DP<-compute_fixed_parameters_1d(fun= function(x) functionPY(x, S,sigma_py=sigma_py),K=K_pr)
    alpha.PY<- compute_parameters_SB_1d(K_pr,S,S,10^4)
    }
    N_eps<-floor(.compute_tau_mean(discount.PY,alpha.PY,0.1) + 2*.compute_tau_var(discount.PY,alpha.PY,0.1))
    N<- max(N_eps,30)
    if (N <= S){
      N=S}
    reductList$N<- N
    cat(c(alpha.PY,discount.PY),"\n alpha and sigma \n")
  }
  #the last values of the parameters are the starting points of the chains (for DRtype 1)
  
  if(!is.null(reductList$N))REDUCT <- T
  
  
  tmp <- .gjamHoldoutSetup(holdoutIndex, holdoutN, n)
  holdoutIndex <- tmp$holdoutIndex; holdoutN <- tmp$holdoutN
  inSamples    <- tmp$inSamples;         nIn <- tmp$nIn
  
  tmp <- .gjamSetup(typeNames, x, y, breakList, holdoutN, holdoutIndex,
                    censor=censor, effort=effort) 
  w <- tmp$w; z <- tmp$z; y <- tmp$y; other <- tmp$other; cuts <- tmp$cuts
  cutLo       <- tmp$cutLo; cutHi <- tmp$cutHi; plo <- tmp$plo; phi <- tmp$phi
  ordCols     <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
  conCols     <- which(typeNames == 'CON')
  classBySpec <- tmp$classBySpec; breakMat <- tmp$breakMat
  minOrd      <- tmp$minOrd; maxOrd <- tmp$maxOrd; censorCA <- tmp$censorCA
  censorDA    <- tmp$censorDA; censorCON <- tmp$censorCON; 
  ncut <- ncol(cuts);  corCols <- tmp$corCols
  catCols     <- which(attr(typeNames,'CATgroups') > 0)
  sampleW     <- tmp$sampleW
  ordShift    <- tmp$ordShift
  
  sampleW[censorCA] <- 1
  sampleW[censorDA] <- 1
  sampleW[censorCON] <- 1
  sampleWhold <- tgHold <- NULL
  wHold <- NULL
  wmax  <- apply(y/effMat,2,max,na.rm=T)
  pmin  <- -2*abs(wmax)
  
  if(mmiss > 0){
    phi[ ymiss ] <- wmax[ ymiss[,2] ]
    plo[ ymiss ] <- pmin[ ymiss[,2] ]
    sampleW[ ymiss ] <- 1
  }
  
  ploHold <- phiHold <- NULL
  
  if(holdoutN > 0){
    sampleWhold <- sampleW[holdoutIndex,]  #to predict X
    sampleW[holdoutIndex,] <- 1
    tgHold  <- cuts
    wHold   <- w[drop=F,holdoutIndex,]
    
    ploHold <- plo[drop=F,holdoutIndex,]   # if LOHI: updated to current yp
    phiHold <- phi[drop=F,holdoutIndex,]
  }
  
  byCol <- byRow <- F
  if(attr(sampleW,'type') == 'cols')byCol <- T
  if(attr(sampleW,'type') == 'rows')byRow <- T
  indexW <- attr(sampleW,'index')
  
  notCorCols <- c(1:S)
  if(length(corCols) > 0)notCorCols <- notCorCols[-corCols]
  
  ############ 'other' columns
  sigmaDf  <- nIn - Q + S - 1
  sg <- diag(.1,S)
  SO <- S
  
  notOther <- c(1:S)
  sgOther  <- NULL
  if(length(other) > 0){                     
    notOther   <- notOther[!notOther %in% other]
    SO         <- length(notOther)
    sg[other,] <- sg[,other] <- 0
    sgOther    <- matrix( cbind(other,other),ncol=2 )
    sg[sgOther] <- .1
  }
  
  ############## prior on beta
  loB <- hiB <- NULL
  beta <- bg <- matrix(0,Q,S)
  rownames(beta) <- colnames(x)
  BPRIOR <- F
  
  if( !is.null(betaPrior) ){
    colnames(betaPrior$lo) <- .cleanNames(colnames(betaPrior$lo))
    colnames(betaPrior$hi) <- .cleanNames(colnames(betaPrior$hi))
    loB <- betaPrior$lo
    hiB <- betaPrior$hi
    
    bg <- (loB + hiB)/2
    bg[is.nan(bg)] <- 0
    
    wB <- which(!is.na(t(loB[,notOther])), arr.ind=T)[,c(2,1)]
    wB <- rbind(wB, which(!is.na(t(hiB[,notOther])), arr.ind=T)[,c(2,1)])
    colnames(wB) <- c('row','col')
    
    tmp <- .betaPrior(bg, notOther, loB, hiB)
    bg <- tmp$beta; loB <- tmp$loB; hiB <- tmp$hiB
    wB <- tmp$wB; BPRIOR <- tmp$BPRIOR
    bg[is.nan(bg)] <- 0
    
    tmp <- .getPattern(bg[,notOther], wB)
    Brows <- tmp$rows
    Bpattern <- tmp$pattern
    BPRIOR <- T
    bg[!is.finite(bg)] <- 0
  }
  
  zeroBeta <- .factorCoeffs2Zero(factorBeta, snames, betaPrior)  # max zero is missing factor level
  zeroLambda <- NULL
  
  ############### time 
  if( TIME ){
    
    BPRIOR <- T
    
    tmp <- .getTimeIndex(timeList, other, notOther, xdata, x, xl, y, w)
    Lmat   <- tmp$Lmat; Lpattern <- tmp$Lpattern;  wL <- tmp$wL
    Vmat   <- tmp$Vmat;    Lrows <- tmp$Lrows; gindex <- tmp$gindex
    loLmat <- tmp$loLmat; hiLmat <- tmp$hiLmat; Arows <- tmp$Arows
    Amat   <- tmp$Amat; Apattern <- tmp$Apattern; wA <- tmp$wA
    Umat   <- tmp$Umat;   uindex <- tmp$uindex
    loAmat <- tmp$loAmat; hiAmat <- tmp$hiAmat; aindex <- tmp$aindex
    Brows  <- tmp$Brows;      bg <- tmp$bg; Bpattern <- tmp$Bpattern
    wB     <- tmp$wB;        loB <- tmp$loB; hiB <- tmp$hiB
    timeZero <- tmp$timeZero; timeLast <- tmp$timeLast
    maxTime  <- tmp$maxTime; inSamples <- tmp$inSamples 
    tindex   <- tmp$tindex; sindex <- tmp$sindex; i1 <- tmp$i1; i2 <- tmp$i2
    
    if(is.null(loB))BPRIOR <- F
    
    Unew <- Umat
    Vnew <- Vmat
    mua  <- mub <- mug <- muw <- w*0
    
    zeroLambda <- .factorCoeffs2Zero(factorLambda, snames, lambdaPrior)
    timeList$lambdaPrior$hi[zeroLambda] <- lambdaPrior$hi[zeroLambda] <- 0
    timeList$betaPrior$hi[zeroBeta]     <- betaPrior$hi[zeroBeta] <- 0
    
    standMatLmat <- Lmat*0
    notStandardLmat <- numeric(0)
    
    if(length(standRowsL) > 0){
      csl <- paste('_',names(standRowsL),sep='')
      for(j in 1:length(csl)){
        wj <- grep(csl[j],rownames(Lmat))
        standMatLmat[wj,] <- standMatL[standRowsL[j],]
        notStandardLmat <- c(notStandardLmat,wj)
      }
    }
  } 
  
  if(byCol){
    inw <- intersect( colnames(y)[indexW], colnames(y)[notOther] )
    indexW <- match(inw,colnames(y)[notOther])
  }
  
  IXX <- NULL
  if(nmiss == 0){
    XX    <- crossprod(x)
    IXX <- chol2inv(chol( XX ) )
  }
  
  
  updateBeta <- .betaWrapper(REDUCT, TIME, BPRIOR, notOther, IXX, 
                             betaLim=max(wmax)/2)
  
  ############ dimension reduction
  
  inSamp <- inSamples
  if(TIME)inSamp <- tindex[,1]     # index for x
  
  CLUST <- T   # dirichlet 
  if(DRtype=="basic") .param.fn <- .paramWrapper(REDUCT, inSamp, SS=length(notOther))
  if(DRtype=="1") .param.fn <- .paramWrapper_1(REDUCT, inSamp, SS=length(notOther))
  if(DRtype=="2") .param.fn <- .paramWrapper_2(REDUCT, inSamp, SS=length(notOther))
  if(DRtype=="3") .param.fn <- .paramWrapper_3(REDUCT, inSamp, SS=length(notOther))
  
  
  sigmaerror <- .1
  
  if(DRtype=="basic")  otherpar   <- list(S = S, Q = Q, sigmaerror = sigmaerror, 
                                          Z = NA, K =rep(1,S), sigmaDf = sigmaDf)
   if(DRtype=="1")  otherpar   <- list(S = S, Q = Q, sigmaerror = sigmaerror, 
                                      Z = NA, K =rep(1,S), sigmaDf = sigmaDf,alpha.DP=alpha.DP,rate=rate,shape=shape, alpha.DP_vec=alpha.DP)
  if(DRtype=="2") otherpar   <- list(S = S, Q = Q, sigmaerror = sigmaerror, 
                                     Z = NA, K =rep(1,S), sigmaDf = sigmaDf,alpha.PY=alpha.PY,discount.PY=discount.PY, matrixCnk = Precomp_matrix, fun_pointer = ptr_logv_comp_mat)
  if(DRtype=="3")  otherpar   <- list(S = S, Q = Q, sigmaerror = sigmaerror, 
                                      Z = NA, K =rep(1,S), sigmaDf = sigmaDf,alpha.PY=alpha.PY,discount.PY=discount.PY)
  
  
  sigErrGibbs <- rndEff <- NULL
  
  yp <- y
  wmax <- ymax <- apply(y,2,max)
  wmax <- wmax/effMat
  
  if(REDUCT){
    cat( paste('\nDimension reduced from',S,'X',S,'->',N,'X',r,'responses\n') )
    otherpar$N <- N; otherpar$r <- r; otherpar$sigmaerror <- 0.1
    otherpar$Z <- rmvnormRcpp(N,rep(0,r),1/S*diag(r))
    otherpar$D <- .riwish(df = (2 + r + N), 
                          S = (crossprod(otherpar$Z) +
                                 2*2*diag(rgamma(r,shape=1,rate=0.001))))
    otherpar$K <- sample(1:N,length(notOther),replace=T)
    if(DRtype=="basic"){ otherpar$alpha.DP <- alpha.DP
    otherpar$pvec     <- .sampleP(N=N, avec=rep(alpha.DP/N,(N-1)),
                                  bvec=((N-1):1)*alpha.DP/N, K=otherpar$K)
    }

    if(DRtype=="1")  {otherpar$alpha.DP <- alpha.DP #initial point for alpha
    otherpar$alpha.DP_vec=alpha.DP
    otherpar$alpha.DP <- alpha.DP
    otherpar$pvec<- .sampleP(N=N, avec=rep(alpha.DP/N,(N-1)),
                             bvec=((N-1):1)*alpha.DP/N, K=otherpar$K)
    otherpar$rate<-rate
    otherpar$shape<-shape
    alpha.DP_g<-rep(0,ng)
    pk_g<-matrix(1,ng,N)
    }
    if(DRtype=="2")  {
      otherpar$discount.PY <-discount.PY
      otherpar$alpha.PY <- alpha.PY
      otherpar$pvec      <-  .sampleP_PYM(N = N, alpha_val = alpha.PY, sigma_val = discount.PY, K = otherpar$K, Mat =Precomp_matrix,  func = ptr_logv_comp_mat)  
      
      otherpar$matrixCnk <- Precomp_matrix
      otherpar$fun_pointer <- ptr_logv_comp_mat
      pk_g<-matrix(1,ng,N)
    }
    
    if(DRtype=="3")  {
    otherpar$discount.PY <-discount.PY
    otherpar$alpha.PY <- alpha.PY
    otherpar$pvec     <- .sampleP(N=N, avec=rep(1-discount.PY,(N-1)),
                                  bvec=(1:(N-1))*discount.PY + alpha.PY, K=otherpar$K)
    pk_g<-matrix(1,ng,N)
    }
  
    kgibbs <- matrix(1,ng,S)
    sgibbs <- matrix(0,ng, N*r)
    nnames <- paste('N',1:N,sep='-')
    rnames <- paste('r',1:r,sep='-')
    colnames(sgibbs) <- .multivarChainNames(nnames,rnames)
    sigErrGibbs <- rep(0,ng)   
    
    rndEff <- w*0
    
  } else {
    Kindex <- which(as.vector(lower.tri(diag(S),diag=T)))
    nK     <- length(Kindex)
    sgibbs <- matrix(0,ng,nK)
    colnames(sgibbs) <- .multivarChainNames(snames,snames)[Kindex] # half matrix
  }
  
  out <- .param.fn(CLUST=T, x, beta = bg[,notOther], Y = w[,notOther], otherpar)  
  sg[notOther,notOther]    <- out$sg
  otherpar      <- out$otherpar
  
  muw <- w
  
  if(!TIME){
    
    Y <- w[inSamp,notOther]
    sig <- sg[notOther,notOther]
    
    if(REDUCT){
      Y <- Y - rndEff[inSamp,notOther]
      sig <- sigmaerror
    }
    
    bg[,notOther] <- updateBeta(X = x[inSamp,], Y, sig, beta = bg[,notOther],
                                loB, hiB)
    muw <- x%*%bg
    
  }else{
    
    mua <- Umat%*%Amat
    mug <- Vmat%*%Lmat
    Y <- w - mua - mug - rndEff
    
    if(REDUCT){
      sig <- sigmaerror
    }else{ sig <- sg[notOther,notOther] }
    
    bg[,notOther] <- updateBeta(X = x[tindex[,2],], Y = Y[tindex[,2],notOther], 
                                sig = sig, beta = bg[,notOther], 
                                lo = loB[,notOther], hi = hiB[,notOther], 
                                rows=Brows, pattern=Bpattern)
    mub <- x%*%bg
    muw <- mub + mug + mua 
    wpropTime <- .001 + .1*abs(w)
  }
  
  sg[other,] <- sg[,other] <- 0
  diag(sg)[other]          <- 1
  rownames(bg)  <- xnames
  rownames(sg)  <- colnames(sg) <- colnames(bg) <- snames
  colnames(x)   <- xnames
  
  
  
  ############ ordinal data
  
  cutg <- tg <- numeric(0)
  
  if('OC' %in% typeCode){
    tg       <- cutg <- cuts
    cnames   <- paste('C',1:ncut,sep='-')
    nor      <- length(ordCols)
    cgibbs   <- matrix(0,ng,(ncut-3)*nor)
    colnames(cgibbs) <- as.vector( outer(snames[ordCols],
                                         cnames[-c(1,2,ncut)],paste,sep='_') )
    tmp   <- .gjamGetCuts(y+1,ordCols)
    cutLo <- tmp$cutLo
    cutHi <- tmp$cutHi
    plo[,ordCols] <- tg[cutLo]                                        
    phi[,ordCols] <- tg[cutHi]
    lastOrd <- ncol(tg)
  }
  
  ############ setup w
  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  typeCols <- tmp$typeCols
  allTypes <- unique(typeCols)
  Y <- w
  
  LOHI <- F
  if(!LOHI & holdoutN > 0){
    minlo <- apply(plo,2,min)
    minlo[minlo > 0] <- 0
    maxhi <- apply(phi,2,max)
  }
  
  if(!TIME){
    .updateW <- .wWrapper(REDUCT, RANDOM, S, effMat, corCols, notCorCols, typeNames, 
                          typeFull, typeCols, 
                          allTypes, holdoutN, holdoutIndex, censor, 
                          censorCA, censorDA, censorCON, notOther, sampleW, 
                          byRow, byCol,
                          indexW, ploHold, phiHold, sampleWhold, inSamp)
  }else{
    
    .updateW <- .wWrapperTime(sampleW, y, timeZero, i1, i2, tindex, gindex,
                              uindex, notOther, n, S, REDUCT, RANDOM)
    Y <- w - mua -  mug - rndEff
  }
  
  ycount <- rowSums(y)
  if('CC' %in% typeCode)ycount <- rowSums(y[,compCols])
  
  ############ X prediction
  
  tmp <- .xpredSetup(Y, x, bg, interBeta$isNonLinX, factorBeta, 
                     factorBeta$intMat, 
                     standMat, standMu, notOther, notStandard) 
  factorBeta$linFactor <- tmp$linFactor; xpred <- tmp$xpred; px <- tmp$px
  lox <- tmp$lox; hix <- tmp$hix
  
  priorXIV  <- diag(1e-5,ncol(x))
  priorX    <- colMeans(x)
  priorX[abs(priorX) < 1e-10] <- 0
  
  linFactor <- NULL
  
  ################## random groups
  
  if('random' %in% names(modelList)){
    
    RANDOM <- T
    rname  <- modelList$random
    randGroupTab <- table( as.character(xdata[,rname]) )
    
    wss <- names(randGroupTab[randGroupTab <= 2])
    if(length(wss) > 0){
      xdata[,rname] <- .combineFacLevels(xdata[,rname], fname=wss, 
                                         aname = 'rareGroups', vminF=1)
      randGroupTab <- table( as.character(xdata[,rname]) )
    }
    
    randGroups <- names( randGroupTab )
    G <- length(randGroups)
    
    groupIndex  <- match(as.character(xdata[,rname]),randGroups)
    rmm <- matrix(groupIndex,length(groupIndex), S)
    smm <- matrix(1:S, length(groupIndex), S, byrow=T)
    
    randGroupIndex <- cbind( as.vector(smm), as.vector(rmm) )
    colnames(randGroupIndex) <- c('species','group')
    xdata[,rname] <- as.factor(xdata[,rname])
    alphaRandGroup <- matrix(0, S, G)
    rownames(alphaRandGroup) <- snames
    colnames(alphaRandGroup) <- randGroups
    Cmat <- var(w[,notOther]/2)
    Cmat <- Cmat + diag(.1*diag(Cmat))
    Cprior <- Cmat
    CImat <- solve(Cprior)
    Ckeep <- diag(S)
    
    alphaRanSums <- alphaRandGroup*0
    groupRandEff <- w*0
    
    Kindex <- which(as.vector(lower.tri(diag(S),diag=T)))
    nK     <- length(Kindex)
    alphaVarGibbs <- matrix(0,ng,nK)
    colnames(alphaVarGibbs) <- .multivarChainNames(snames,snames)[Kindex] # half matrix
  }
  
  
  ################################## XL prediction: variables in both
  
  if(TIME){
    
    tmp <- .xpredSetup(Y, xl, lambdaPrior$lo, 
                       interLambda$isNonLinX, factorLambda, interLambda$intMat, standMatL, 
                       standMuL, notOther, notStandardL) 
    factorLambda$linFactor <- tmp$linFactor
    lox <- c(lox,tmp$lox[!names(tmp$lox) %in% names(lox)])
    hix <- c(hix,tmp$lox[!names(tmp$hix) %in% names(hix)])
    
    ################ or
    xpred <- cbind(xpred,xl[,!colnames(xl) %in% colnames(x)])
    Qall <- ncol(xpred) - 1
    
    intMat <- numeric(0)
    if( length(interBeta$intMat) > 0 ){
      intMat <- match(xnames[interBeta$intMat],colnames(xpred))
      intMat <- matrix(intMat,nrow(interBeta$intMat),3)
    }
    if( length(interLambda$intMat) > 0){
      ib <- match(xlnames[interLambda$intMat],colnames(xpred))
      ib <- matrix(ib,nrow(interLambda$intMat),3)
      intMat <- rbind(intMat,ib)
    }
    
    linFactor <- numeric(0)
    lf <- factorBeta$linFactor
    if( length(lf) > 0 ){
      for(k in 1:length(lf)){
        kf <- match(xnames[lf[[k]]],colnames(xpred))
        linFactor <- append(linFactor,list(kf))
      }
    }
    lf <- factorLambda$linFactor
    if( length(lf) > 0 ){
      for(k in 1:length(lf)){
        kf <- match(xlnames[lf[[k]]],colnames(xpred))
        linFactor <- append(linFactor,list(kf))
      }
    }
  }
  
  ############  contrasts, predict F matrix
  
  tmp <- .setupFactors(xdata, xnames, factorBeta)
  ff  <- factorBeta[names(factorBeta) != 'factorList']
  factorBeta <- append(ff,tmp)
  
  ############ E matrix
  emat <- matrix(0,S,S)
  colnames(emat) <- rownames(emat) <- snames
  lo <- hi <- lm <- hm <- ess <- emat
  
  fmat <- factorBeta$fmat
  fnames <- rownames( factorBeta$lCont )
  q2 <- nrow(fmat)
  
  if(TIME){
    tmp <- .setupFactors(xdata, xlnames, factorLambda)
    ff <- factorLambda[names(factorLambda) != 'factorList']
    if(length(tmp) > 0)factorLambda <- append(ff,tmp)
    factorLambda$LCONT <- rep(TRUE, factorLambda$nfact)
    flnames <- rownames( factorLambda$lCont )
    
    ############ E matrix TIME
    ematL <- matrix(0,S,S)
    colnames(ematL) <- rownames(ematL) <- snames
    essL <- ematL
  }
  
  ############ sp richness
  richness <- richFull <- NULL
  RICHNESS <- F
  
  inRichness <- which(!typeNames %in% c('CON','CAT','OC'))
  inRichness <- inRichness[!inRichness %in% other]
  if(length(inRichness) > 2)RICHNESS  <- T
  
  wrich <- y*0 
  wrich[,inRichness] <- 1
  wrich[ymiss] <- 0
  
  presence <- w*0
  
  covx <- cov(x)
  
  ############ sums
  predx  <- predx2 <- xpred*0
  yerror <- ypred  <- ypred2 <- wpred  <- wpred2 <- ymissPred <- ymissPred2 <- y*0
  sumDev <- 0   #for DIC
  sMean  <- sg*0
  ntot   <- 0
  
  if(nmiss > 0){
    xmissSum <- xmissSum2 <- rep(0,nmiss)
  }
  
  if(TIME)predxl <- predxl2 <- xl*0
  
  ############ gibbs chains
  
  q2 <- length(fnames)
  fSensGibbs <- matrix(0,ng,q2)
  colnames(fSensGibbs) <- fnames
  
  bFacGibbs <- matrix(0,ng,q2*SO)
  colnames(bFacGibbs) <- .multivarChainNames(fnames,snames[notOther])
  
  bgibbs <- matrix(0,ng,S*Q)
  colnames(bgibbs) <- .multivarChainNames(xnames,snames)
  bgibbsUn <- bgibbs                   # unstandardized
  
  covE <- cov( x%*%factorBeta$dCont )  # note that x is standardized
  
  if(TRAITS){
    
    specTrait <- specByTrait[colnames(y),]
    tnames    <- colnames(specTrait)
    M         <- ncol(specTrait)
    specTrait <- t(specTrait)
    
    tpred  <- tpred2 <- matrix(0,n,M)
    
    missTrait <- which(is.na(specTrait),arr.ind=T)
    if(length(missTrait) > 0){
      traitMeans <- rowMeans(specTrait,na.rm=T)
      specTrait[missTrait] <- traitMeans[missTrait[,2]]
      warning( paste('no. missing trait values:',nrow(missTrait)) )
    }
    
    bTraitGibbs <- matrix(0,ng,M*Q)
    colnames(bTraitGibbs) <- .multivarChainNames(xnames,tnames)
    
    bTraitFacGibbs <- matrix(0,ng,q2*M)
    colnames(bTraitFacGibbs) <- .multivarChainNames(fnames,tnames)
    
    mgibbs <- matrix(0,ng,M*M)
    colnames(mgibbs) <- .multivarChainNames(tnames,tnames)
  }
  
  if(TIME){
    
    yy <- y*0
    yy[rowInserts,] <- 1
    ymiss <- which(yy == 1, arr.ind=T)
    rm(yy)
    mmiss <- length(ymiss)
    
    covL <- cov( xl%*%factorLambda$dCont )  # note x is standardized
    
    ggibbs <- matrix(0,ng,nrow(wL))
    colnames(ggibbs) <- rownames(wL)
    
    wnames <- apply(wA,1,paste0,collapse='-')  #locations in Amat, not alpha
    alphaGibbs <- matrix(0,ng,nrow(wA))
    colnames(alphaGibbs) <- wnames
    
    nl <- nrow(lambda)
    lgibbs <- matrix(0,ng,length(lambda[,notOther]))
    colnames(lgibbs) <- .multivarChainNames(xlnames,snames[notOther])
    
    gsensGibbs <- matrix(0,ng,nl)
    colnames(gsensGibbs) <- rownames(lambda)
    
    asensGibbs <- matrix(0,ng,nrow(Amat))
    colnames(asensGibbs) <- rownames(Amat)
    
    ni  <- length(i1)
    
    nA <- nrow(wA)
    nL <- nrow(wL)
    
    spA <- rep(.001, nA)
    spL <- rep(.01, nL)
    g1 <- 1
    gcheck <- c(50, 100, 200, 400, 800)
    tinyg <- 1e-6
  }
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  # unstandardize
  
  tmp <- .getUnstandX(x, standRows, standMu[,1],standMat[,1],
                      interBeta$intMat)
  S2U      <- tmp$S2U
  xUnstand <- tmp$xu
  
  if(TIME){
    tmp <- .getUnstandX(xl, standRowsL, standMuL[,1],standMatL[,1],
                        interLambda$intMat)
    S2UL      <- tmp$S2U
    xlUnstand <- tmp$xu
  }
  
  if(REDUCT){
    rndTot <- w*0 
  }
  notPA <- which(!typeNames == 'PA' & !typeNames == 'CON')
  
  
  if(length(y) < 10000 | FULL) FULL <- T
  
  if(FULL){
    ygibbs <- matrix(0,ng,length(y))
  }
  if(RICHNESS){
    ypredPres <- ypredPres2 <- ypredPresN <- y*0
    shannon   <- rep(0,n)
  }
  
  for(g in 1:ng){ ########################################################
    
    if(REDUCT){
      
      #   if(g > burnin)CLUST <- F
      
      Y <- w[,notOther]
      if(RANDOM)Y <- Y - groupRandEff[,notOther] 
      if(TIME)  Y <- Y - mua[,notOther] - mug[,notOther] 
      
      tmp <- .param.fn(CLUST=T, x, beta = bg[,notOther], Y = Y, otherpar)
      sg[notOther,notOther] <- tmp$sg
      otherpar            <- tmp$otherpar
      rndEff[,notOther]   <- tmp$rndEff
      sigmaerror          <- otherpar$sigmaerror
      kgibbs[g,notOther]  <- otherpar$K
      sgibbs[g,]          <- as.vector(otherpar$Z)
      sigErrGibbs[g]      <- sigmaerror
      if(DRtype=="1") {alpha.DP_g[g]<- otherpar$alpha.DP
      pk_g[g,]<-otherpar$pvec}
      if(DRtype=="3")  {pk_g[g,]<-otherpar$pvec}
      if(DRtype=="2")  {pk_g[g,]<-otherpar$pvec}
      if(length(corCols) > 0){
        if(max(diag(sg)[corCols]) > 5){  #overfitting covariance
          stop(
            paste('\noverfitted covariance, reductList$N = ',N, 
                  'reductList$r = ',r, '\nreduce N, r\n')
          )
        }
      }
      
      sg[sgOther]         <- .1*sigmaerror
      
      sinv <- .invertSigma(sg[notOther,notOther],sigmaerror,otherpar,REDUCT)
      sdg  <- sqrt(sigmaerror)
      
      if(!TIME){
        Y <- w[inSamp,notOther] - rndEff[inSamp,notOther]
        if(RANDOM)Y <- Y - groupRandEff[inSamp,notOther]
        bg[,notOther] <- updateBeta(X = x[inSamp,], Y, 
                                    sig = sigmaerror, beta = bg[,notOther],
                                    lo=loB[,notOther], hi=hiB[,notOther])
        muw[inSamp,] <- x[inSamp,]%*%bg
        
      } else {
        
        mua  <- Umat%*%Amat 
        mug  <- Vmat%*%Lmat
        
        Y <- w[,notOther] - mua[,notOther] - mug[,notOther] - rndEff[,notOther]
        if(RANDOM)Y <- Y - groupRandEff[,notOther]
        bg[,notOther] <- updateBeta(X = x[tindex[,2],], Y = Y[tindex[,2],], 
                                    sig = sigmaerror, beta = bg[,notOther],
                                    rows = Brows, pattern = Bpattern,
                                    lo=loB[,notOther], hi=hiB[,notOther])
        mub <- x%*%bg
        Y   <- w - mub - mua - rndEff
        if(RANDOM)Y <- Y - groupRandEff
        
        Lmat[,notOther] <- updateBeta(X = Vmat[tindex[,2],], 
                                      Y = Y[tindex[,2],notOther], sig=sigmaerror, 
                                      beta = Lmat[,notOther],
                                      rows = Lrows, pattern = Lpattern, 
                                      lo=loLmat, hi=hiLmat, ixx=F)
        
        #     Lmat[,notOther] <- .updateBetaMet(X = Vmat[tindex[,2],], 
        #                                       Y[tindex[,2],notOther], 
        #                                       B = Lmat[,notOther],
        #                           lo=loLmat, hi=hiLmat, loc = wL, REDUCT, 
        #                           sig=sigmaerror,sp=spL)
        mug  <- Vmat%*%Lmat
        Y    <- w - mub - mug - rndEff
        if(RANDOM)Y <- Y - groupRandEff
        Amat <- updateBeta(X = Umat[tindex[,2],], Y[tindex[,2],], sig=sigmaerror, 
                           rows = Arows, pattern = Apattern, 
                           beta = Amat,
                           lo=loAmat, hi=hiAmat, ixx=F)
        #     Amat <- .updateBetaMet(X = Umat[tindex[,2],], Y[tindex[,2],notOther], 
        #                                       B = Amat,
        #                                       lo=loAmat, hi=hiAmat, loc = wA, REDUCT, 
        #                                      sig=sigmaerror,sp=rexp(nA,1/spA))
        mua <- Umat%*%Amat
        
        #     if(g %in% gcheck){
        #       g2   <- g - 1
        #       spA <- apply(alphaGibbs[g1:g2,],2,sd)/2 + tinyg
        #       spL <- apply(ggibbs[g1:g2,],2,sd)/2 + tinyg
        #       if(g < 200)g1 <- g
        #     }
        
        muw <- mub + mug + mua + rndEff
      }
      
    } else {
      Y <- w[inSamp,notOther]
      if(RANDOM)Y <- Y - groupRandEff[inSamp,notOther]
      bg[,notOther] <- updateBeta(X = x[inSamp,], Y, 
                                  sig = sg[notOther,notOther], 
                                  beta = bg[,notOther], 
                                  lo=loB, hi=hiB)
      
      muw[inSamp,] <- x[inSamp,]%*%bg
      
      SS   <- crossprod(w[inSamp,] - muw[inSamp,])
      SI   <- solveRcpp(SS[notOther,notOther])
      sinv <- .rwish(sigmaDf,SI)
      
      sg[notOther,notOther] <- solveRcpp(sinv)
      sgibbs[g,] <- sg[Kindex]
    }
    
    # muw does not include rndEff or groupRandEff
    
    alphaB <- .sqrtRootMatrix(bg,sg,DIVIDE=T)
    
    if( 'OC' %in% typeCode ){
      tg <- .updateTheta(w,tg,cutLo,cutHi,ordCols,
                         holdoutN,holdoutIndex,minOrd,maxOrd) # var scale
      cutg <- .gjamCuts2theta(tg,ss = sg[ordCols,ordCols]) # corr scale
      breakMat[ordCols,1:lastOrd] <- cutg
      cgibbs[g,] <- as.vector( cutg[,-c(1,2,ncut)] )
      
      plo[,ordCols] <- cutg[cutLo]
      phi[,ordCols] <- cutg[cutHi]
    }
    
    if(RANDOM){
      
      cw <- w - muw
      if(REDUCT){
        cw <- cw - rndEff
        v  <- 1/sigmaerror*.byGJAM(as.vector(cw), randGroupIndex[,1], 
                                   randGroupIndex[,2], alphaRandGroup*0, 
                                   fun='sum')[notOther,]
        sinv <- diag(1/sigmaerror, SO)
      }else{
        v <- .byGJAM(as.vector(cw), randGroupIndex[,1], 
                     randGroupIndex[,2], alphaRandGroup*0, fun='sum')[notOther,]
        v <- sinv%*%v
      }
      
      alphaRandGroup[notOther,] <- randEffRcpp(v, randGroupTab, 
                                               sinv, CImat)
      if(length(other) > 0)alphaRandGroup[other,] <- 0
      if(g < 100){
        alphaRandGroup[notOther,] <- 
          sweep( alphaRandGroup[notOther,], 2, 
                 colMeans(alphaRandGroup[notOther,]), '-')
      }
      SS  <- crossprod(t(alphaRandGroup[notOther,]))
      SS  <- S*SS + Cmat
      
      testv <- try( chol(SS) ,T)
      if( inherits(testv,'try-error') ){
        tiny  <- .1*diag(SS)
        SS  <- SS + diag(diag(SS + tiny))
      }
      
      Ckeep[notOther,notOther] <- .riwish( df = S*G + 1, SS )
      CImat <- solveRcpp(Ckeep[notOther,notOther])
      
      alphaVarGibbs[g,] <- Ckeep[Kindex]
      groupRandEff <- t(alphaRandGroup)[groupIndex,]
    }
    
    if(TIME){
      
      #    muw does not include groupRandEff
      tmp <- .updateW(w,plo,phi,wpropTime,xl,yp,Lmat,Amat,mub,rndEff, groupRandEff,
                      sdg,muw,Umat,Vmat,sinv)
      w <- tmp$w; muw <- tmp$muw; yp <- tmp$yp; Umat <- tmp$Umat; Vmat <- tmp$Vmat
      
      groups <- NULL
      
      for(k in allTypes){
        
        wk <- which(typeCols == k)
        nk <- length(wk)
        wo <- which(wk %in% notOther)
        wu <- which(typeCols[notOther] == k)
        wp <- w[, wk, drop=F]
        yp <- yp[, wk, drop=F]
        
        if(typeFull[wk[1]] == 'countComp')groups <- CCgroups
        if(typeFull[wk[1]] == 'fracComp')groups  <- FCgroups
        if(typeFull[wk[1]] == 'categorical')groups <- CATgroups
        
        glist <- list(wo = wo, type = typeFull[wk[1]], yy = y[,wk,drop=F], 
                      wq = wp, yq = yp, cutg = cutg, 
                      censor = censor, censorCA = censorCA, 
                      censorDA = censorDA, censorCON  = censorCON, 
                      eff = effMat[,wk,drop=F], groups = groups, 
                      k = k, typeCols = typeCols, notOther = notOther, 
                      wk = wk, sampW = sampleW[,wk])
        tmp <- .gjamWLoopTypes( glist )
        w[,wk]  <- tmp[[1]]
        yp[,wk] <- tmp[[2]]
      }
      
      #predict X
      
      ww <- w
      ww[ww < 0] <- 0
      mua  <- Umat%*%Amat
      mug  <- Vmat%*%Lmat
      muw <- mua + mub + mug + rndEff
      
      xtmp <- xpred
      xtmp[,-1] <- .tnorm(n*Qall,-3,3,xpred[,-1],.1)
      
      # factors
      if( length(linFactor) > 0 ){
        
        for(k in 1:length(linFactor)){
          
          mm  <- linFactor[[k]]
          wcol <- sample(mm,n,replace=T)
          xtmp[,mm[-1]] <- 0
          xtmp[ cbind(1:n, wcol) ] <- 1
          
        }
      }
      
      if(length(intMat) > 0){     #  interactions
        xtmp[,intMat[,1]] <- xtmp[,intMat[,2]]*xtmp[,intMat[,3]]
      }
      
      ae     <- mua + rndEff
      Vnow   <- Vmat
      mubNow <- xpred[,xnames]%*%bg
      mubNew <- xtmp[,xnames]%*%bg
      
      Vnow[tindex[,2],] <- ww[tindex[,1],gindex[,'colW']]*
        xpred[tindex[,2],xlnames][,gindex[,'rowG']]
      Vnow[timeZero+1,] <- ww[timeZero,gindex[,'colW']]*
        xpred[timeZero+1,xlnames][,gindex[,'rowG']]
      mugNow <- Vnow%*%Lmat
      muNow  <- mubNow + mugNow + ae
      
      Vnew[tindex[,2],] <- ww[tindex[,1],gindex[,'colW']]*
        xtmp[tindex[,2],xlnames][,gindex[,'rowG']]
      Vnew[timeZero+1,] <- ww[timeZero,gindex[,'colW']]*
        xtmp[timeZero+1,xlnames][,gindex[,'rowG']]
      mugNew <- Vnew%*%Lmat
      muNew  <- mubNew + mugNew + ae
      
      if(REDUCT){
        pnow <- dnorm(w[,notOther],muNow[,notOther],sdg,log=T)
        pnew <- dnorm(w[,notOther],muNew[,notOther],sdg,log=T)
        a1   <- exp( rowSums(pnew - pnow) )
      }else{
        pnow <- .dMVN(w[tindex[,2],notOther],muNow,sinv=sinv,log=T) 
        pnew <- .dMVN(w[tindex[,2],notOther],muNew,sinv=sinv,log=T) 
        a1   <- exp(pnew - pnow)
      }
      z    <- runif(length(a1),0,1)
      za   <- which(z < a1)
      if(length(za) > 0){
        xpred[za,] <- xtmp[za,]
        Vmat[za,] <- Vnew[za,]
        muw[za,]  <- muNew[za,]
        mub[za,]  <- mubNew[za,]
        mug[za,]  <- mugNew[za,]
      }
      
      if(nlmiss > 0)xl[xlmiss] <- xpred[xmiss]
      
      if(nmiss > 0){
        
        x[xmiss] <- xpred[xmiss]
        
        tmp    <- .getUnstandX(x, standRows, standMu[,1],
                               standMat[,1], intMat)            
        S2U    <- tmp$S2U
        XX     <- crossprod(x)
        IXX    <- solveRcpp(XX)
      }
      
      ggibbs[g,]     <- Lmat[wL]
      alphaGibbs[g,] <- Amat[wA]
      
    } else{ #############not TIME
      
      tmp   <- .updateW( rows=1:n, x, w, y, bg, sg, alpha=alphaB, 
                         cutg, plo, phi, rndEff, groupRandEff, 
                         sigmaerror, wHold )
      w     <- tmp$w
      yp    <- tmp$yp
      plo   <- tmp$plo
      phi   <- tmp$phi
      wHold <- tmp$wHold    #values for w if not held out
      
      
      Y <- w[,notOther]
      if(holdoutN > 0) Y[holdoutIndex,] <- wHold[,notOther]  # if w not held out
      if(RANDOM)Y <- Y - groupRandEff[,notOther]
      
      if(nmiss > 0){
        
        x[xmiss] <- .imputX_MVN(x,Y,bg[,notOther],xmiss,sinv,xprior=xprior,
                                xbound=xbound)[xmiss]
        tmp      <- .getUnstandX(x, standRows, standMu[,1],
                                 standMat[,1], intMat)            
        S2U    <- tmp$S2U
        XX     <- crossprod(x)
        IXX    <- solveRcpp(XX)
      }
      
      if( PREDICTX & length(predXcols) > 0){
        
        if( length(interBeta$isNonLinX) > 0 ){
          
          xpred <- .predictY2X_nonLinear(xpred, yy=Y,bb=bg[,notOther],
                                         ss=sg[notOther,notOther],
                                         priorIV = priorXIV,priorX=priorX,
                                         factorObject = factorBeta, interObject = interBeta,
                                         lox, hix)$x
        }
        
        if( length(px) > 0 ){
          wn <- which(!is.finite(xpred),arr.ind=T)
          if(length(wn) > 0){
            tmp <- matrix(priorX,Q,nrow(wn))
            xpred[wn[,1],] <- t(tmp)
          }
          xpred[,px] <- .predictY2X_linear(xpred, yy=Y, bb=bg[,notOther],
                                           ss=sg[notOther,notOther], sinv = sinv,
                                           priorIV = priorXIV, 
                                           priorX=priorX,predCols=px, 
                                           REDUCT=REDUCT, lox, hix)[,px]
          wn <- which(!is.finite(xpred),arr.ind=T)
          if(length(wn) > 0){
            tmp <- matrix(priorX,Q,nrow(wn))
            xpred[wn[,1],] <- t(tmp)
          }
        }
        
        if( length(factorBeta$linFactor) > 0 ){
          
          # predict all factors
          xtmp <- xpred
          xtmp[,factorBeta$findex] <- 
            .predictY2X_linear(xpred, yy=Y, 
                               bb=bg[,notOther],
                               ss=sg[notOther,notOther], sinv = sinv,
                               priorIV = priorXIV, 
                               priorX=priorX,predCols=factorBeta$findex, 
                               REDUCT=REDUCT, lox, hix)[,factorBeta$findex]
          for(k in 1:length(factorBeta$linFactor)){
            
            mm  <- factorBeta$linFactor[[k]]
            tmp <- xtmp[,mm]
            
            tmp[,1] <- 0
            ix  <- apply(tmp,1,which.max)   
            
            tmp <- tmp*0
            tmp[ cbind(1:n,ix) ] <- 1
            tmp <- tmp[,-1,drop=F]
            xpred[,mm[-1]] <- tmp
          }
        }
        xpred[,1] <- 1
      }
    }
    
    setTxtProgressBar(pbar,g)
    
    bgu <- bg                    # unstandardize beta
    if(length(standRows) > 0){
      if(TIME){
        bgu <- S2U%*%mub
        lambda[ gindex[,c('rowG','colW')]] <- Lmat[wL]
        lambdas <- S2UL%*%mug      # unstandardized lambda
        lgibbs[g,] <- lambdas[,notOther]
      }else{
        bgu <- S2U%*%x%*%bg
      }
    }
    
    bgibbsUn[g,] <- bgu          # unstandardized
    bgibbs[g,]   <- bg           # standardized
    
    # Fmatrix centered for factors, 
    # bg is standardized by x, bgu is unstandardized
    
    tmp <- .contrastCoeff(beta=bg[,notOther], 
                          notStand = notStandard[notStandard %in% xnames], 
                          sigma = sg[notOther,notOther], sinv = sinv,
                          stand = standMat, factorObject=factorBeta )
    agg   <- tmp$ag
    beg   <- tmp$eg
    fsens <- tmp$sens
    
    fSensGibbs[g,]  <- sqrt(diag(fsens))
    bFacGibbs[g,] <- agg       # stand for X and W, centered for factors
    
    if(TRAITS){
      Atrait <- bg%*%t(specTrait[,colnames(yp)])  # standardized
      Strait <- specTrait[,colnames(yp)]%*%sg%*%t(specTrait[,colnames(yp)])
      bTraitGibbs[g,] <- Atrait
      mgibbs[g,] <- Strait
      
      minv <- ginv(Strait)
      
      tmp <- .contrastCoeff(beta=Atrait, 
                            notStand = notStandard[notStandard %in% xnames], 
                            sigma = Strait, sinv = minv,
                            stand = standMat, factorObject=factorBeta )
      tagg   <- tmp$ag
      bTraitFacGibbs[g,] <- tagg # stand for X and W, centered for factors
    }
    
    if(TIME){
      
      tmp <- .contrastCoeff(beta=lambda[,notOther], 
                            notStand = notStandardL[notStandardL %in% xlnames], 
                            sigma = sg[notOther,notOther],sinv = sinv,
                            stand=standMatL, factorObject=factorLambda)
      lgg   <- tmp$ag
      leg   <- tmp$eg
      lsens <- tmp$sens
      
      lss <- sqrt(diag(lsens))
      
      if(g == 1){
        if( !all(names(lss) %in% colnames(gsensGibbs)) )
          colnames(gsensGibbs) <- names(lss)
      }
      
      gsensGibbs[g,names(lss)] <- lss
      
      alpha[ aindex[,c('toW','fromW')] ] <- Amat[wA]
      asens <- Amat[,notOther]%*%sinv%*%t(Amat[,notOther])
      asens <- sqrt(diag(asens))
      asensGibbs[g,] <- asens
    }
    
    if(FULL)ygibbs[g,] <- as.vector(yp)
    
    if(g > burnin){
      
      ntot   <- ntot + 1
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      
      tmp <- .dMVN(w[,notOther], muw[,notOther], sg[notOther,notOther], log=T)
      
      sumDev <- sumDev - 2*sum(tmp) 
      yerror <- yerror + (yp - y)^2
      
      fmat <- fmat + fsens
      
      sMean  <- sMean + sg
      
      wpred  <- wpred + w
      wpred2 <- wpred2 + w^2
      
      if(RICHNESS){
        
        yy <- yp
        
        if('PA' %in% typeNames){
          wpa <- which(typeNames[inRichness] == 'PA')
          yy[,inRichness[wpa]] <- round(yp[,inRichness[wpa]]) #######
        }
        
        if(length(notPA) > 0){
          w0 <- which(yy[,notPA] <= 0)
          w1 <- which(yy[,notPA] > 0)
          yy[,notPA][w0] <- 0
          yy[,notPA][w1] <- 1
        }
        
        shan <- sweep(yy[,inRichness], 1, rowSums(yy[,inRichness]), '/')
        shan[shan == 0] <- NA
        shan <- -rowSums(shan*log(shan),na.rm=T)
        shannon <- shannon + shan
        
        wpp <- which(yy > 0)
        ypredPres[wpp]  <- ypredPres[wpp] + yp[wpp]
        ypredPres2[wpp] <- ypredPres2[wpp] + yp[wpp]^2
        ypredPresN[wpp] <- ypredPresN[wpp] + 1
        
        presence[,inRichness] <- presence[,inRichness] + yy[,inRichness]
        ones <- round(rowSums(yy[,inRichness]))
        more <- round(rowSums(yy[,inRichness]*wrich[,inRichness,drop=F]))
        richFull <- .add2matrix(ones,richFull)
        richness <- .add2matrix(more,richness)  # only for non-missing
      }
      
      if(RANDOM){
        alphaRanSums <- alphaRanSums + alphaRandGroup
      }
      
      if(mmiss > 0){
        ymissPred[ymiss]  <- ymissPred[ymiss] + y[ymiss]
        ymissPred2[ymiss] <- ymissPred2[ymiss] + y[ymiss]^2
      }
      if(nmiss > 0){
        xmissSum  <- xmissSum + x[xmiss]
        xmissSum2 <- xmissSum2 + x[xmiss]^2
      }
      
      if(PREDICTX & length(predXcols) > 0){
        predx  <- predx + xpred
        predx2 <- predx2 + xpred^2
      }
      
      wa0 <- which(colSums(agg) != 0)
      ess[notOther[wa0],notOther[wa0]]  <- 
        t(agg[,wa0,drop=F])%*%covE%*%agg[,wa0,drop=F] 
      if(TIME){
        wa0 <- which(colSums(lgg) != 0)
        ess[notOther[wa0],notOther[wa0]]  <- 
          ess[notOther[wa0],notOther[wa0]] +
          t(lgg[,wa0,drop=F])%*%covL%*%lgg[,wa0,drop=F] 
      }
      
      emat[notOther[wa0],notOther[wa0]] <- 
        emat[notOther[wa0],notOther[wa0]] + 
        .cov2Cor( ess[notOther[wa0],notOther[wa0]] )
      
      lo[ ess < 0 ] <- lo[ ess < 0 ] + 1
      hi[ ess > 0 ] <- hi[ ess > 0 ] + 1
      
      ess[notOther,notOther] <- ginv(ess[notOther,notOther])
      
      lm[ ess < 0 ] <- lm[ ess < 0 ] + 1  # neg values
      hm[ ess > 0 ] <- hm[ ess > 0 ] + 1  # pos values
      
      if(REDUCT){
        rndTot <- rndTot + rndEff
      }
      
      if(TRAITS){
        yw     <- sweep(yp,1,rowSums(yp),'/')
        yw[yw <= 0]   <- 0
        yw[is.na(yw)] <- 0
        Ttrait <- .gjamPredictTraits(yw,specTrait[,colnames(yp)], traitTypes)
        tpred  <- tpred + Ttrait
        tpred2 <- tpred2 + Ttrait^2
      }
    }
  }     
  
  ################# end gibbs loop ####################
  
  
  otherpar$S <- S 
  otherpar$Q <- Q
  otherpar$snames <- snames
  otherpar$xnames <- xnames
  
  presence <- presence/ntot
  
  if(RICHNESS){
    missRows <- sort(unique(ymiss[,1]))
    richNonMiss <- richness/ntot            #only non-missing plots
    yr  <- as.matrix(ydata[,inRichness]) 
    yr[yr > 0] <- 1
    yr <- rowSums(yr,na.rm=T)
    vv  <- matrix(as.numeric(colnames(richNonMiss)),n,
                  ncol(richNonMiss),byrow=T)
    rmu <- rowSums( vv * richNonMiss )/rowSums(richNonMiss)
    
    rsd <- sqrt( rowSums( vv^2 * richNonMiss )/rowSums(richNonMiss) - rmu^2)
    
    vv  <- matrix(as.numeric(colnames(richFull)),n,ncol(richFull),byrow=T)
    rfull <- rowSums( vv * richFull )/rowSums(richFull)
    rfull[missRows] <- NA
    rmu <- rowSums(presence)
    
    shan <- sweep(y[,inRichness], 1, rowSums(y[,inRichness]), '/')
    shan[shan == 0] <- NA
    shanObs <- -rowSums(shan*log(shan),na.rm=T)
    
    richness <- cbind(yr, rmu, rsd, rfull, shanObs, shannon/ntot )
    colnames(richness) <- c('obs','predMu','predSd','predNotMissing',
                            'H_obs', 'H_pred')
    if(TIME)richness[timeZero,] <- NA
    
    ypredPresMu  <- ypredPres/ypredPresN   #predictive mean and se given presence
    ypredPresMu[ypredPresN == 0] <- 0
    yvv <- ypredPres2/ypredPresN - ypredPresMu^2
    yvv[!is.finite(yvv)] <- 0
    ypredPresSe <- sqrt(yvv)
  }
  
  if('OC' %in% typeNames){
    ordMatShift <- matrix(ordShift,n,length(ordCols),byrow=T)
    onames <- snames[ordCols]
    wb <- match(paste(onames,'intercept',sep='_'), colnames(bgibbs))
    bgibbs[,wb] <- bgibbs[,wb] + matrix(ordShift,ng,length(ordCols),byrow=T)
    bgibbsUn[,wb] <- bgibbsUn[,wb] + matrix(ordShift,ng,length(ordCols),byrow=T)
    y[,ordCols] <- y[,ordCols] + ordMatShift
  }
  
  if(mmiss > 0){
    ymissPred[ymiss]  <- ymissPred[ymiss]/ntot
    yd <- ymissPred2[ymiss]/ntot - ymissPred[ymiss]^2
    yd[!is.finite(yd)| yd < 0] <- 0
    ymissPred2[ymiss] <- sqrt(yd)
    
    if('OC' %in% typeNames){
      ymissPred[,ordCols] <- ymissPred[,ordCols] + ordMatShift
    }
  }
  
  xunstand    <- .getUnstandX(x, standRows, standMu[,1],
                              standMat[,1], interBeta$intMat)$xu
  
  rmspeBySpec <- sqrt( colSums(yerror)/ntot/n )
  rmspeAll    <- sqrt( sum(yerror)/ntot/n/S )
  
  sMean <- sMean/ntot
  
  if(TIME){
    
    xtime <- xpred*0
    xtime[,xnames] <- x
    xtime[,xlnames] <- xl
    
    xlunstand    <- .getUnstandX(xl, standRowsL, standMuL[,1],
                                 standMatL[,1], interLambda$intMat)$xu
    xtimeUn <- xtime*0
    xtimeUn[,xnames]  <- xunstand
    xtimeUn[,xlnames] <- xlunstand
    
    loL <- hiL <- lambdaMuUn <- lambdaSeUn <- lambda*0
    tmp1 <- colMeans(ggibbs[burnin:ng,])    #unstandardized
    tmp2 <- apply(ggibbs[burnin:ng,],2,sd)
    lambdaMuUn[ gindex[,c('rowG','colW')] ] <- tmp1
    lambdaSeUn[ gindex[,c('rowG','colW')] ] <- tmp2
    loL[gindex[,c('rowG','colW')] ] <- loLmat[wL]
    hiL[gindex[,c('rowG','colW')] ] <- hiLmat[wL]
    
    loA <- hiA <- alphaMu <- alphaSe <- matrix(0,S,S)
    tmp1 <- colMeans(alphaGibbs[burnin:ng,])    #unstandardized
    tmp2 <- apply(alphaGibbs[burnin:ng,],2,sd)
    alphaMu[ aindex[,c('toW','fromW')] ] <- tmp1
    alphaSe[ aindex[,c('toW','fromW')] ] <- tmp2
    loA[ aindex[,c('toW','fromW')] ] <- loAmat[wA]
    hiA[ aindex[,c('toW','fromW')] ] <- hiAmat[wA]
    
    gsensMu <- colMeans(gsensGibbs[burnin:ng,]) 
    gsensSd <- apply(gsensGibbs[burnin:ng,],2,sd)
    asensMu <- colMeans(asensGibbs[burnin:ng,])
    asensSd <- apply(asensGibbs[burnin:ng,],2,sd)
  }
  
  tmp <- .chain2tab(bgibbs[burnin:ng,], snames, xnames)
  betaStandXmu <- tmp$mu
  betaStandXTable <- tmp$tab
  
  tmp <- .chain2tab(bgibbsUn[burnin:ng,], snames, xnames)
  betaMu <- tmp$mu
  betaTable <- tmp$tab
  
  tmp <- .chain2tab(bFacGibbs[burnin:ng,], snames[notOther], rownames(agg))
  betaStandXWmu <- tmp$mu
  betaStandXWTable <- tmp$tab
  
  tmp <- .chain2tab(fSensGibbs[burnin:ng,,drop=F])
  sensTable <- tmp$tab[,1:4]
  
  yMu <- ypred/ntot
  y22 <- ypred2/ntot - yMu^2
  y22[y22 < 0] <- 0
  ySd <- sqrt(y22)
  
  cMu <- cuts
  cSe <- numeric(0)
  
  wMu <- wpred/ntot
  wpp <- pmax(0,wpred2/ntot - wMu^2)
  wSd <- sqrt(wpp)
  
  if('OC' %in% typeNames){
    yMu[,ordCols] <- yMu[,ordCols] + ordMatShift
    wMu[,ordCols] <- wMu[,ordCols] + ordMatShift
  }
  
  meanDev <- sumDev/ntot
  
  tmp <- .dMVN(wMu[,notOther],x%*%betaMu[,notOther],
               sMean[notOther,notOther], log=T)
  pd  <- meanDev - 2*sum(tmp )
  DIC <- pd + meanDev
  
  yscore <- colSums( .getScoreNorm(y[,notOther],yMu[,notOther],
                                   ySd[,notOther]^2),na.rm=T )  # gaussian w
  xscore <- xpredMu <- xpredSd <- NULL
  standX <- xmissMu <- xmissSe <- NULL
  
  if(RANDOM){
    ns <- 500
    simIndex <- sample(burnin:ng,ns,replace=T)
    tmp <- .expandSigmaChains(snames, alphaVarGibbs, otherpar, simIndex=simIndex,
                              sigErrGibbs, kgibbs, REDUCT=F)
    alphaRandGroupVarMu <- tmp$sMu
    alphaRandGroupVarSe <- tmp$sSe
    alphaRandByGroup <- alphaRanSums/ntot
    
  }
  
  if(PREDICTX){
    xpredMu <- predx/ntot
    xpredSd <- predx2/ntot - xpredMu^2
    xpredSd[xpredSd < 0] <- 0
    xpredSd <- sqrt(xpredSd)
    
    xrow <- standRows
    xmu  <- standMu[,1]
    xsd  <- standMat[,1]
    
    if(TIME){
      xrow <- c(standRows, standRowsL) 
      ww   <- !duplicated(names(xrow))
      xrow <- names(xrow)[ww]
      xmu  <- c(standMu[xrow,1], standMuL[xrow,1])
      xsd  <- c(standMat[xrow,1],standMatL[xrow,1])
      #  xrow <- names(xrow)[ww]
      #  xrow <- match(xrow,colnames(xpredMu))
      #  names(xrow) <- colnames(xpredMu)[xrow]
    }
    
    xpredMu <- .getUnstandX(xpredMu, xrow, xmu, xsd, intMat)$xu
    xpredSd[,xrow] <- xpredSd[,xrow]*matrix( xsd[xrow], n, length(xrow),
                                             byrow=T ) 
    
    if(TIME){
      if(Q == 2)xscore <- mean( .getScoreNorm(xtime[,2],
                                              xpredMu[,2],xpredSd[,2]^2) )
      if(Q > 2)xscore <- colMeans(.getScoreNorm(xtime[,-1],
                                                xpredMu[,-1],xpredSd[,-1]^2) )
    }else{
      if(Q == 2)xscore <- mean( .getScoreNorm(x[,2],
                                              xpredMu[,2],xpredSd[,2]^2) )
      if(Q > 2)xscore <- colMeans(.getScoreNorm(x[,-1],
                                                xpredMu[,-1],xpredSd[,-1]^2) )
    }
    
    if(TIME){
      wz <- wMu
      wz[wz < 0] <- 0
      Vmat[tindex[,2],] <- wz[tindex[,2], 
                              gindex[,'colW']]*xl[tindex[,2], gindex[,'colX']]
      Vmat[timeZero,]   <- wz[timeZero, 
                              gindex[,'colW']]*xl[timeZero, gindex[,'colX']]
      Umat <- wz[,uindex[,1]]*wz[,uindex[,2]] 
      
      Amat[ aindex[,c('rowA','fromW')] ] <- alphaMu[ aindex[,c('toW','fromW')] ]
      Lmat[ gindex[,c('rowL','colW')] ] <- lambdaMuUn[ gindex[,c('rowG','colW')] ]
      
      muw <- x%*%betaMu[,notOther] + Vmat%*%Lmat[,notOther] + Umat%*%Amat[,notOther]
      
      tmp <- .dMVN(wMu[,notOther],muw[,notOther],
                   sMean[notOther,notOther], log=T )
      pd  <- meanDev - 2*sum(tmp )
      DIC <- pd + meanDev
    }
  }
  
  if(nmiss > 0){
    xmissMu <- xmissSum/ntot
    xmissSe <- sqrt( xmissSum2/ntot - xmissMu^2 )
  }
  
  if(length(standRows) > 0){                #unstandardize
    standX <- cbind(standMu[,1],standMat[,1])
    colnames(standX) <- c('xmean','xsd')
    rownames(standX) <- rownames(standMat)
  }
  
  # betaSens, sigma and R
  
  ns <- 500
  simIndex <- sample(burnin:ng,ns,replace=T)
  
  tmp <- .expandSigmaChains(snames, sgibbs, otherpar, simIndex=simIndex,
                            sigErrGibbs, kgibbs, REDUCT)
  corMu <- tmp$rMu; corSe <- tmp$rSe
  sigMu  <- tmp$sMu; sigSe  <- tmp$sSe
  
  whichZero <- which(lo/ntot < ematAlpha & 
                       hi/ntot < ematAlpha,arr.ind=T) #not different from zero
  whConZero <- which(lm/ntot < ematAlpha & 
                       hm/ntot < ematAlpha,arr.ind=T)
  
  ematrix  <- emat/ntot
  fmatrix  <- fmat/ntot
  
  tMu <- tSd <- tMuOrd <- btMu <- btSe <- stMu <- stSe <- numeric(0)
  
  if(TRAITS){
    
    tMu <- tpred/ntot
    tSd <- sqrt(tpred2/ntot - tMu^2)
    wo  <- which(traitTypes == 'OC')    #predict ordinal scores
    M   <- ncol(tMu)
    
    if(length(wo) > 0){
      tMuOrd <- tMu*0
      for(j in wo)tMuOrd[,j] <- round(tMu[,j],0) - 1
      tMuOrd <- tMuOrd[,wo]
    }
    
    tmp <- .chain2tab(bTraitGibbs[burnin:ng,], tnames, xnames) #standardized
    betaTraitXMu <- tmp$mu
    betaTraitXTable <- tmp$tab
    
    tmp <- .chain2tab(mgibbs[burnin:ng,], tnames, tnames) 
    varTraitMu <- tmp$mu
    varTraitTable <- tmp$tab
    
    tmp <- .chain2tab(bTraitFacGibbs[burnin:ng,], tnames, rownames(tagg) )
    betaTraitXWmu <- tmp$mu
    betaTraitXWTable <- tmp$tab
  }
  
  if('OC' %in% typeNames){
    nk  <- length(ordCols)
    nc  <- ncut - 3
    
    os <- rep(ordShift,nc)
    
    cgibbs <- cgibbs + matrix(os,ng,length(os),byrow=T)
    
    tmp <- .processPars(cgibbs)$summary
    cMu <- matrix(tmp[,'estimate'],nk,nc)
    cSe <- matrix(tmp[,'se'],nk,ncut-3)
    cMu <- cbind(ordShift,cMu)
    cSe <- cbind(0,cSe)
    colnames(cMu) <- colnames(cSe) <- cnames[-c(1,ncut)]
    rownames(cMu) <- rownames(cSe) <- snames[ordCols]
    breakMat[ordCols,c(2:(2+(ncol(cMu))-1))] <- cMu
  }
  
  if('PA' %in% typeNames){
    zMu <- yMu
    zSd <- ySd
  }
  ##To change line after progress bar
  cat('\n')
  # outputs
  if(length(reductList) == 0)reductList <- list(N = 0, r = 0)
  reductList$otherpar <- otherpar
  
  modelList$effort    <- effort;      modelList$formula <- formula
  modelList$typeNames <- typeNames;    modelList$censor <- censor
  modelList$effort    <- effort; modelList$holdoutIndex <- holdoutIndex
  modelList$REDUCT    <- REDUCT;       modelList$TRAITS <- TRAITS
  modelList$ematAlpha <- ematAlpha; modelList$traitList <- traitList
  modelList$reductList <- reductList; modelList$ng <- ng
  modelList$burnin <- burnin
  
  inputs <- list(xdata = xdata, x = xunstand, standX = standX,
                 standMat = standMat, standRows = standRows, y = y, 
                 notOther = notOther, other = other, breakMat = breakMat, 
                 designTable = designTable, classBySpec = classBySpec, 
                 factorBeta = factorBeta, interBeta = interBeta,
                 linFactor = linFactor, intMat = intMat, RANDOM = RANDOM)
  missing <- list(xmiss = xmiss, xmissMu = xmissMu, xmissSe = xmissSe, 
                  ymiss = ymiss, ymissMu = ymissPred, ymissSe = ymissPred2)
  parameters <- list(betaMu = betaMu, betaTable = betaTable, 
                     betaStandXmu = betaStandXmu, 
                     betaStandXTable = betaStandXTable,
                     betaStandXWmu =  betaStandXWmu,
                     betaStandXWTable = betaStandXWTable,
                     corMu = corMu, corSe = corSe, 
                     sigMu = sigMu, sigSe = sigSe, 
                     ematrix = ematrix, fmatrix = fmatrix,
                     whichZero = whichZero, whConZero = whConZero,
                     wMu = wMu, wSd = wSd, sensTable = sensTable)
  prediction <- list(presence = presence, xpredMu = xpredMu, xpredSd = xpredSd,
                     ypredMu = yMu, ypredSd = ySd, richness = richness)
  chains <- list(sgibbs = sgibbs, bgibbs = bgibbs, bgibbsUn = bgibbsUn,
                 fSensGibbs = fSensGibbs, bFacGibbs = bFacGibbs) 
  fit <- list(DIC = DIC, yscore = yscore, 
              xscore = xscore, rmspeAll = rmspeAll,
              rmspeBySpec = rmspeBySpec)
  if(FULL)chains <- append(chains, list(ygibbs = ygibbs))
  if(RANDOM){
    parameters <- append(parameters,
                         list( randGroupVarMu = alphaRandGroupVarMu,
                               randGroupVarSe = alphaRandGroupVarSe,
                               randByGroup = alphaRandByGroup) )
  }
  if(RICHNESS){
    prediction <- append(prediction, 
                         list(yPresentMu = ypredPresMu, yPresentSe = ypredPresSe))
  }
  if(REDUCT) {
    parameters <- append(parameters, list(rndEff = rndTot/ntot))#, specRand = specRand))
    
    if(DRtype=="basic") chains <- append(chains,list(kgibbs = kgibbs, sigErrGibbs = sigErrGibbs))
    if(DRtype=="1") chains <- append(chains,list(kgibbs = kgibbs, sigErrGibbs = sigErrGibbs,alpha.DP_g=alpha.DP_g, pk_g=pk_g))
    if(DRtype=="2")  chains <- append(chains,list(kgibbs = kgibbs, sigErrGibbs = sigErrGibbs, pk_g=pk_g))
    if(DRtype=="3") chains <- append(chains,list(kgibbs = kgibbs, sigErrGibbs = sigErrGibbs, pk_g=pk_g))
    
  }
  
  if('OC' %in% typeNames){
    parameters <- c(parameters,list(cutMu = cMu, cutSe = cSe))
    chains <- c(chains,list(cgibbs = cgibbs))
    modelList <- c(modelList,list(yordNames = yordNames))
  }
  
  if(TRAITS){
    parameters <- c(parameters,
                    list(betaTraitXMu = betaTraitXMu, 
                         betaTraitXTable = betaTraitXTable,
                         varTraitMu = varTraitMu, 
                         varTraitTable = varTraitTable,
                         betaTraitXWmu = betaTraitXWmu,
                         betaTraitXWTable = betaTraitXWTable))
    prediction <- c(prediction, list(tMuOrd = tMuOrd, tMu = tMu, tSe = tSd))
    chains <- append( chains,list(bTraitGibbs = bTraitGibbs,
                                  bTraitFacGibbs = bTraitFacGibbs,
                                  mgibbs = mgibbs) ) 
  }
  if(TIME){
    inputs <- c(inputs, list(xtime = xtime, timeZero = timeZero,
                             interLambda = interLambda, 
                             factorLambda = factorLambda))
    chains <- c(chains, list(ggibbs = ggibbs, alphaGibbs = alphaGibbs,
                             gsens = gsensGibbs, asens = asensGibbs))
    parameters <- c(parameters, 
                    list(lambdaMuUn = lambdaMuUn, lambdaSeUn = lambdaSeUn, 
                         lambdaLo = loL, lambdaHi = hiL,
                         alphaMu = alphaMu, alphaSe = alphaSe,
                         alphaLo = loA, alphaHi = hiA,
                         gsensMu = gsensMu, gsensSe = gsensSd, 
                         asensMu = asensMu, asensSe = asensSd,
                         aindex = aindex, wA = wA, unidex = uindex))
  }
  
  chains     <- chains[ sort( names(chains) )]
  fit        <- fit[ sort( names(fit) )]
  inputs     <- inputs[ sort( names(inputs) )]
  missing    <- missing[ sort( names(missing) )]
  modelList  <- modelList[ sort( names(modelList) )]
  parameters <- parameters[ sort( names(parameters) )]
  prediction <- prediction[ sort( names(prediction) )]
  
  all <- list(chains = chains, fit = fit, inputs = inputs, missing = missing,
              modelList = modelList, parameters = parameters,
              prediction = prediction)
  all$call <- match.call()
  all <- all[ sort(names(all)) ]
  class(all) <- "gjam"
  
  all
}

.contrastCoeff <- function(beta, sigma, sinv, notStand, stand, factorObject,
                           conditional=NULL){ 
  
  # if(!is.null(notStand)){
  #   beta[notStand,] <- beta[notStand,]*stand[notStand,]
  # }
  SO  <- ncol(beta)
  
  
  agg <- .sqrtRootMatrix(beta,sigma,DIVIDE=T)  #cor/stand scale
  
  if(factorObject$nfact > 0){          # center factors
    agg <- factorObject$lCont%*%agg    # standardized x, cor scale for w
    for(k in 1:factorObject$nfact){
      f2  <- factorObject$facList2[[k]]
      fk  <- paste(names(factorObject$facList2)[k],f2,sep='')
      
      amu <- colMeans(agg[drop=F,fk,])
      nl  <- length(fk)
      agg[fk,] <- agg[fk,] - matrix(amu,nl,SO,byrow=T)
      egg <- agg
    }
  } else {
    agg <- agg[drop=F,-1,]
    egg <- agg
  }
  if(is.null(conditional)){
    
    sens <- egg%*%sinv%*%t(egg)
    
  }else{
    con <- which(colnames(beta) %in% conditional)
    nc  <- c(1:SO)[-con]
    sg  <- sigma[con,con] - 
      sigma[con,nc]%*%solve(sigma[nc,nc])%*%sigma[nc,con]
    sens <- egg[,con]%*%solve(sg)%*%t(egg[,con])
  }
  
  
  list(ag = agg, eg = egg, sens = sens)
}


.chain2tab <- function(chain, snames = NULL, xnames = NULL){
  
  mu <- colMeans(chain)  
  SE <- apply(chain,2,sd)
  CI <- apply(chain,2,quantile,c(.025,.975))
  splus <- rep('', length=length(SE))
  splus[CI[1,] > 0 | CI[2,] < 0] <- '*'
  
  tab <- cbind( mu, SE, t(CI))
  tab <- signif(tab, 3)
  colnames(tab) <- c('Estimate','SE','CI_025','CI_975')
  tab <- as.data.frame(tab)
  tab$sig95 <- splus
  attr(tab, 'note') <- '* indicates that zero is outside the 95% CI'
  
  if(!is.null(snames)){
    Q <- length(xnames)
    S <- length(snames)
    
    mu <- matrix(mu,Q,S)
    colnames(mu) <- snames
    rownames(mu) <- xnames 
    mu <- signif(mu, 3)
  }
  
  list(mu = mu, tab = tab)
}


summary.gjam <- function(object,...){ 
  
  TRAITS <- F
  
  beta   <- object$parameters$betaMu # not standardized
  rb <- rownames(beta)
  cb <- colnames(beta)
  S  <- ncol(beta)
  Q  <- nrow(beta)
  n  <- nrow(object$inputs$y)
  notOther <- object$inputs$notOther
  other    <- object$inputs$other
  
  ng <- object$modelList$ng
  burnin <- object$modelList$burnin
  
  if("betaTraitTable" %in% names(object$parameters))TRAITS <- T
  
  sens <- .chain2tab(object$chains$fSensGibbs[burnin:ng,])$tab[,1:4]
  
  RMSPE       <- object$fit$rmspeBySpec
  imputed <- rep(0,S)
  missingx <- rep(0,Q)
  if(length(object$missing$xmiss) > 0){
    xtab <- table(object$missing$xmiss[,2])
    missingx[ as.numeric(names(xtab)) ] <- xtab
  }
  if(length(object$missing$ymiss) > 0){
    xtab <- table(object$missing$ymiss[,2])
    imputed[ as.numeric(names(xtab)) ] <- xtab
  }
  
  RMSPE <- RMSPE[notOther]
  #  imputedY <- imputed[notOther]
  bb   <- t( signif(rbind(beta[,notOther], RMSPE),3) )
  # imputedX <- c(missingx,NA,NA)
  # bb <- cbind(bb,imputedX)
  
  # cc <- as.vector( signif(beta[,notOther], 3) )
  # ss <- object$parameters$betaSe[,notOther]
  # ss <- as.vector( signif(ss, 3) )
  # rr <- as.vector( t(outer(cb[notOther],rb,paste,sep='_')) )
  # TAB <- data.frame(Estimate = cc, StdErr = ss)
  # rownames(TAB) <- rr
  
  # qb <- t( apply(object$chains$bgibbsUn,2,quantile,c(.025,.975)) )
  # mq <- match(rr,rownames(qb))
  # ci <- signif(qb[mq,],3)
  # TAB <- cbind(TAB,ci)
  
  cat("\nSensitivity by predictor variables f:\n")
  print( sens )
  
  cat("\nCoefficient matrix B:\n")
  print( t(bb) )
  
  cat("\nCoefficient matrix B:\n")
  print(object$parameters$betaTable)
  cat("\nLast column indicates if 95% posterior distribution contains zero.\n")
  
  cat("\nCoefficient matrix B, standardized for X:\n")
  print(object$parameters$betaStandXtable)
  cat("\nLast column indicates if 95% posterior distribution contains zero.\n")
  
  cat("\nCoefficient matrix B, standardized for X and W:\n")
  print(object$parameters$betaStandXWtable)
  cat("\nLast column indicates if 95% posterior distribution contains zero.\n")
  
  
  if(TRAITS){
    cat("\nCoefficient matrix for traits:\n")
    print(object$parameters$betaTraitTable)
    cat("\nLast column indicates if 95% posterior distribution contains zero.\n")
    
    cat("\nCoefficient matrix for traits, standardized for X and W:\n")
    print(object$parameters$betaTraitXWTable)
    cat("\nLast column indicates if 95% posterior distribution contains zero.\n")
    
  }
  
  
  if( length(object$modelSummary$missFacSpec) > 0 ){
    cat("\nMissing factor combinations:\n")
    print(object$modelSummary$missFacSpec)
  }
  
  dt <- object$input$designTable[-2,]
  cat("\n Design Table\n")
  print(dt)
  
  words <- .summaryWords(object)
  cat("\n",words)
  
  res <- list(DIC=object$fit$DIC, sensitivity = sens, 
              Coefficients=bb)
  class(res) <- "summary.gjam"
  invisible(res) 
}

.summaryWords <- function(object){
  
  Q  <- ncol(object$inputs$x)
  n  <- nrow(object$inputs$y)
  S  <- ncol(object$inputs$y)
  other    <- object$inputs$other
  notOther <- object$inputs$notOther
  
  nfact <- object$inputs$factorBeta$nfact
  nxmiss <- nrow( object$missing$xmiss )
  nymiss <- nrow( object$missing$ymiss )
  nholdout <- length(object$modelList$holdoutIndex)
  types    <- unique(object$modelList$typeNames)
  if(length(types) == 1)types <- rep(types,S)
  
  ef <- ""
  if( 'DA' %in% types ){
    wd  <- which(types == 'DA')
    rf  <- object$modelList$effort$values[,wd]
    gf <- signif( range(rf), 2)
    wr  <- signif( range(object$inputs$y[,wd]/rf), 3)
    if(gf[1] == gf[2]){
      ab <- paste(" DA effort is ",gf[1]," for all observations. ",sep="")
    }else{
      ab <- paste(" DA effort ranges from ", gf[1], " to ", gf[2],".",sep="")
    }
    ef <- paste(ab, " DA counts per effort (W) ranges from ", 
                wr[1], " to ", wr[2], ".",sep="")
  }
  if( 'CC' %in% types ){
    wd  <- which(types == 'CC')
    rr  <- round( range(object$inputs$y[,wd]) )
    ef <- paste(ef, " CC count range is (", rr[1],", ", rr[2], ").", sep="")
  }
  oc <- ""
  if(length(other) > 0){
    oc <- paste(" 'other' class detected in ydata, '", 
                colnames(object$inputs$y)[other], " column ",
                "', not fitted. ",sep='')
  }
  
  fc <- ""
  if(nfact > 0)fc <- paste(" There are",nfact,"factors in X. ")
  
  ty <- paste0( unique(types), collapse=", ")
  
  words <- paste("Sample contains n = ", n, " observations on S = ",
                 S, " response variables, and ", Q - 1, 
                 " predictors.  Data types (typeNames) include ", ty,
                 ".", fc, ef, oc, " There are ", nxmiss, 
                 " missing values in X and ", nymiss, 
                 " missing values in Y. The RMSPE is ",
                 signif(object$fit$rmspeAll,3),
                 ", and the DIC is ",round(object$fit$DIC),".", sep="")
  dr <- ""
  if(object$modelList$REDUCT){
    nr <- object$chains$kgibbs
    nd <- t( apply(nr,1,duplicated) )
    nr[!nd] <- 0
    nr[nr > 0] <- 1
    nk <- rowSums(1 - nr)
    nk <- max(nk[object$modelList$burnin:object$modelList$ng])
    
    dr <- paste("  Dimension reduction was implemented with N = ",nk,
                " and r = ",object$modelList$reductList$r,".",
                sep="")
  }
  ho <- ""
  if(nholdout > 0)ho <- paste(" Held out were",nholdout,"observations.")
  comp <- paste(" Computation involved ", object$modelList$ng,
                " Gibbs steps, with a burnin of ", object$modelList$burnin, 
                ".",dr,ho,sep='')
  paste(words, comp)
}

print.gjam <- function(x, ...){
  
  summary.gjam(x)
  
}

.getSigTable <- function(chain, SS, QQ, xn, sn){
  
  bci  <- apply(chain,2,quantile,c(.025,.975))
  tmp  <- .between(rep(0,SS*QQ),bci[1,],bci[2,],OUT=T)
  ii <- rep(' ',SS*QQ)
  ii[tmp[bci[1,tmp] < 0]] <- '-'
  ii[tmp[bci[2,tmp] > 0]] <- '+'
  
  bTab <- data.frame( matrix(ii,QQ,SS) )
  colnames(bTab) <- sn
  rownames(bTab) <- xn
  
  bTab <- data.frame( t(bTab) )
  
  bTab
}

.getPlotLayout <- function(np){
  
  # np - no. plots
  
  if(np == 1)return( c(1,1) )
  if(np == 2)return( c(1,2) )
  if(np == 3)return( c(1,3) )
  if(np <= 4)return( c(2,2) )
  if(np <= 6)return( c(2,3) )
  if(np <= 9)return( c(3,3) )
  if(np <= 12)return( c(3,4) )
  if(np <= 16)return( c(4,4) )
  if(np <= 20)return( c(4,5) )
  if(np <= 25)return( c(5,5) )
  if(np <= 25)return( c(5,6) )
  return( c(6,6) )
}


sqrtSeq <- function(maxval){ #labels for sqrt scale
  
  # maxval on sqrt scale
  
  by <- 2
  if(maxval >= 5)   by <- 10 
  if(maxval >= 10)  by <- 20 
  if(maxval >= 20)  by <- 100 
  if(maxval >= 30)  by <- 200 
  if(maxval >= 50)  by <- 500
  if(maxval >= 70)  by <- 1000
  if(maxval >= 100) by <- 2000
  if(maxval >= 200) by <- 10000
  if(maxval >= 500) by <- 50000
  if(maxval >= 700) by <- 100000
  if(maxval >= 1000)by <- 200000
  if(maxval >= 1500)by <- 400000
  
  labs <- seq(0, maxval^2, by = by)
  at   <- sqrt(labs)
  
  list(at = at, labs = labs)
  
}


.plotObsPred <- function(obs, yMean, ySE=NULL, opt = NULL){
  
  nbin <- nPerBin <- xlimit <- ylimit <- NULL
  add <- log <- SQRT <- F
  xlabel <- 'Observed'
  ylabel <- 'Predicted'
  trans <- .4
  col <- 'black'
  bins <- NULL
  atx <- aty <- labx <- laby <- NULL
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  if(!is.null(bins))nbin <- length(bins)
  
  if(log & SQRT)stop('cannot have both log and SQRT scale')
  
  yMean <- as.matrix(yMean)
  obs   <- as.matrix(obs)
  
  if(SQRT){
    xlim <- sqrt(xlimit)
    ylim <- sqrt(ylimit)
    obs   <- as.vector(sqrt(obs))
    yMean <- as.vector(sqrt(yMean))
    xlimit <- sqrt(range(obs,na.rm=T))
    xlimit[2] <- xlimit[2]*2
    ylimit <- sqrt(range(yMean,na.rm=T))
    ylimit[2] <- 2*ylimit[2]
    
    maxy <- max(yMean,na.rm=T)
    maxx   <- max(obs,na.rm=T)
    maxval <- max( c(maxx, maxy) )
    
    tt   <- sqrtSeq(1.2*maxx)
    if(is.null(atx))atx   <- tt$at
    if(is.null(labx))labx <- tt$labs
    
    tt   <- sqrtSeq(1.2*maxy)
    if(is.null(aty))aty   <- tt$at
    if(is.null(laby))laby <- tt$labs
    
    if(ylimit[2] < xlimit[2]) ylimit[2] <- xlimit[2]
    if(xlimit[2] < xlim[2])   xlimit[2] <- xlim[2]
    if(ylimit[2] < ylim[2])   ylimit[2] <- ylim[2]
  }
  
  if(is.null(xlimit))xlimit <- range(obs)
  if(is.null(ylimit) & !add){                      # can only happen if !SQRT
    if(!log){
      plot(obs,yMean,col=.getColor('black',.2),cex=.3, xlim=xlimit,
           xlab=xlabel,ylab=ylabel)
      if(log) suppressWarnings( plot(obs,yMean,col=.getColor('black',.2),cex=.3,
                                     xlim=xlimit,xlab=xlabel,ylab=ylabel,log='xy') )
    }
  }
  
  if(!is.null(ylimit)){
    if(!log & !add){
      if(!SQRT){
        plot(obs,yMean,col=.getColor('black',trans),cex=.3,
             xlim=xlimit,xlab=xlabel,ylab=ylabel,ylim=ylimit)
      }else{
        plot(obs,yMean,col=.getColor('black',trans),cex=.3,
             xlim=xlimit,xlab=xlabel,ylab=ylabel,ylim=ylimit,
             xaxt='n',yaxt='n')
        
        axis(1, at = atx, labels = labx)
        axis(2, at = aty, labels = laby, las=2)
      }
    }
    if(log & !add) plot(obs,yMean,col=.getColor('black',trans),cex=.3,
                        xlim=xlimit,xlab=xlabel,ylab=ylabel,log='xy',ylim=ylimit)
  }
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),
                                 col='grey',lwd=2)
  }
  
  if( !is.null(nbin) | !is.null(nPerBin) ){
    
    if(is.null(bins)){
      nbin <- 20
      bins <- seq(min(obs,na.rm=T),max(obs,na.rm=T),length=nbin)
    }else{
      nbin <- length(bins)
    }
    
    if(!is.null(nPerBin)){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm=T)
      bins <- bins[!duplicated(bins)]
      nbin <- length(bins)
    }
    
    yk <- findInterval(obs,bins)
    yk[yk == nbin] <- nbin - 1
    yk[yk == 1] <- 2
    
    wide <- diff(bins)/2
    db   <- 1
    for(k in 2:(nbin-1)){
      
      qk <- which(is.finite(yMean) & yk == k)
      q  <- quantile(yMean[qk],c(.5,.025,.158,.841,.975),na.rm=T)
      
      if(!is.finite(q[1]))next
      if(q[1] == q[2])next
      
      ym <- mean(yMean[qk])
      xx <- mean(bins[k:(k+1)])
      rwide <- wide[k]
      if(k == 2 & nbin < 5){
        xx <- mean(bins[1:2]) 
        rwide <- wide[1]
      }
      
      if(k > 1)db <- bins[k] - bins[k-1]
      
      if( xx > (bins[k] + db) ){
        xx <- bins[k] + db
        rwide <- wide[ max(c(1,k-1)) ]
      }
      
      suppressWarnings(
        arrows(xx, q[2], xx, q[5], lwd=2, angle=90, code=3, col=.getColor(col,.8),
               length=.02)
      )
      lines(c(xx-.5*rwide,xx+.5*rwide),q[c(1,1)],lwd=2, 
            col=.getColor(col,.8))
      rect(xx-.4*rwide,q[3],xx+.4*rwide,q[4], col=.getColor(col,.5))
    }
  }
  invisible( list(atx = atx, labx = labx, aty = aty, laby = laby) )
}


.gjamPlot <- function(output, plotPars){
  
  PLOTALLY <- TRAITS <- GRIDPLOTS <- SAVEPLOTS <- 
    REDUCT <- TV <- SPECLABS <- SMALLPLOTS <- F
  PREDICTX <- BETAGRID <- PLOTY <- PLOTX <- 
    CORLINES <- SIGONLY <- CHAINS <- RANDOM <- T
  omitSpec   <- trueValues  <- censor <- otherpar <- ng <- NULL
  traitList  <- specByTrait <- typeNames <- classBySpec <- 
    x <- y   <- burnin      <- richness <- betaTraitMu <-  
    corSpec  <- cutMu       <- ypredMu <- DIC <- yscore <- missingIndex <- 
    xpredMu  <- plotByTrait <- tMu <- tMuOrd <- traitTypes <- 
    isFactor <- betaMu <- betaMuUn <- corMu <- modelSummary <-
    randByGroup <- randGroupVarMu <- NULL
  unstandardX  <- NULL 
  ematAlpha     <- .5  
  ematrix      <- NULL
  ymiss <- eCont <- modelList <- timeList <- timeZero <- NULL
  random <- NULL
  kgibbs <- NULL
  chains <- inputs <- parameters <- prediction <- reductList  <- 
    bgibbs <- sgibbs <- sigErrGibbs <- factorBeta <- gsens <- 
    bFacGibbs <- alphaGibbs <- times <- alphaMu <- lambdaMu <-  
    factorLambda <- lambdaMuUn <- notOther <- other <- xtime <- NULL
  holdoutN <- 0
  
  TIME <- F
  
  cex <- 1
  holdoutIndex <- numeric(0)
  clusterIndex <- clusterOrder <- numeric(0)
  
  ncluster   <- min(c(4,ncol(y)))
  
  outFolder <- 'gjamOutput'
  outfile   <- character(0)
  
  width <- height <- 3
  
  oma <- c(1,1,0,0)
  mar <- c(1,1,1,0)
  tcl <- -0.1
  mgp <- c(0,0,0)
  
  specColor <- traitColor <- textCol <- 'black'
  
  for(k in 1:length(output))assign( names(output)[k], output[[k]] )
  
  for(k in 1:length(chains))assign( names(chains)[k], chains[[k]] )
  for(k in 1:length(fit))assign( names(fit)[k], fit[[k]] )
  for(k in 1:length(inputs))assign( names(inputs)[k], inputs[[k]] )
  for(k in 1:length(missing))assign( names(missing)[k], missing[[k]] )
  for(k in 1:length(modelList))assign( names(modelList)[k], modelList[[k]] )
  for(k in 1:length(parameters))assign( names(parameters)[k], parameters[[k]] )
  for(k in 1:length(prediction))assign( names(prediction)[k], prediction[[k]] )
  
  for(k in 1:length(reductList))assign( names(reductList)[k], reductList[[k]] )
  
  if(!is.null(plotPars))for(k in 1:length(plotPars))assign( names(plotPars)[k], plotPars[[k]] )
  
  if( !is.null(traitList) ){
    TRAITS <- T
    for(k in 1:length(traitList))assign( names(traitList)[k], traitList[[k]] )
  }
  if( 'trueValues' %in% names(plotPars) ){
    TV <- T
    for(k in 1:length(trueValues))assign( names(trueValues)[k], trueValues[[k]] )
    
    matchTrue <- match(colnames(betaMu),colnames(beta))
    beta      <- beta[,matchTrue]
    sigma     <- sigma[matchTrue,matchTrue]
    corSpec   <- corSpec[matchTrue,matchTrue]
  }
  if(!is.null(timeList)){
    for(k in 1:length(timeList))assign( names(timeList)[k], timeList[[k]] )
    ypredMu[timeZero,] <- NA
    TIME <- T
  }
  if(length(xpredMu) == 0)PREDICTX <- F
  if(!PREDICTX)PLOTX <- F
  
  if(!is.null(random)){
    RANDOM <- T
  }
  
  oma <- c(0,0,0,0)
  mar <- c(4,4,2,1)
  tcl <- -0.5
  mgp <- c(3,1,0)
  
  if(SAVEPLOTS){
    ff <- file.exists(outFolder)
    if(!ff)dir.create(outFolder)
  }
  
  chainNames <- names(chains)
  allTypes   <- unique(typeNames)
  
  ntypes   <- length(allTypes)
  typeCode <- match(typeNames,allTypes)
  specs    <- rownames(classBySpec)
  Q        <- ncol(x)
  nhold   <- length(holdoutIndex)
  ncut    <- ncol(classBySpec) + 1
  S       <- ncol(y)
  n       <- nrow(y)
  snames  <- colnames(y)
  xnames  <- colnames(x)
  # ng      <- nrow(chains$bgibbs)
  gindex  <- burnin:ng
  
  if(S < 20)SPECLABS <- T
  if(S > 10)CORLINES <- F
  if(S < 8){
    if(GRIDPLOTS)message('no GRIDPLOTS if S < 8')
    GRIDPLOTS <- F
  }
  
  if(length(specColor) == 1)specColor <- rep(specColor, S)
  boxCol    <- .getColor(specColor,.4)
  
  omit     <- c(which(colnames(y) %in% omitSpec),other)
  notOmit  <- 1:S
  SO       <- length(notOther)
  if(length(omit) > 0)notOmit <- notOmit[-omit]
  SM       <- length(notOmit)
  
  snames  <- colnames(y)
  xnames  <- colnames(x)
  
  cvec    <- c('black','brown','orange')
  if(ntypes > 4)cvec <- c(cvec,'green','blue')
  colF    <- colorRampPalette(cvec)
  
  ## richness prediction
  
  xSd <- sqrt( diag(cov(x)) )
  
  HOLD <- F
  if(holdoutN > 0)HOLD <- T
  
  if( !TRAITS & !is.null(richness) ){
    
    if(TIME)richness[timeZero,] <- NA
    
    w1 <- which(richness[,1] > 0)        # these are missing data
    if(HOLD)w1 <- w1[!w1 %in% holdoutIndex]
    
    xlimit <- range(richness[w1,1])
    
    if(diff(xlimit) > 0){
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'richness.pdf') )
      
      par(mfrow=c(1,2), bty='n', omi=c(.3,.3,0,0), mar=c(3,2,2,1), 
          tcl= tcl, mgp=mgp)
      
      xc <- c('obs','H_obs')
      yc <- c('predMu','H_pred')
      
      for(k in 1:2){
        
        kx <- richness[w1,xc[k]]
        ky <- richness[w1,yc[k]]
        if(k == 2){
          kx <- exp(kx)
          ky <- exp(ky)
        }
        
        ylimit <- range(ky)
        
        rr   <- range(kx,na.rm=T)
        bins <- seq(rr[1] - .5, ceiling(rr[2] + .5), by=1)
        nbin <- length(bins)
        
        rh <- hist(kx,bins,plot=F)
        xy <- rbind(c(bins[1],bins,bins[nbin]),c(0,rh$density,0,0))
        
        xy     <- .gjamBaselineHist(kx,bins=bins)
        xy[2,] <- ylimit[1] + .3*xy[2,]*diff(ylimit)/max(xy[2,])
        
        
        plot(xy[1,],xy[2,],col='tan',type='s',lwd=2, ylim=ylimit,
             xlab=' ',ylab='')
        #     axis(1,at=rr[1]:rr[2])
        
        polygon(xy[1,],xy[2,],border='brown',col='tan')
        
        if(HOLD){
          xhold <- richness[holdoutIndex,xc[k]]
          yhold <- richness[holdoutIndex,yc[k]]
          if(k == 2){
            xhold <- exp(xhold)
            yhold <- exp(yhold)
          }
          points(xhold,yhold, col='brown', cex=.3)
        }
        
        opt <- list(log=F, bins = bins,
                    nbin=nbin, xlabel='', ylabel='', col='darkblue', 
                    add=T)
        tmp <- .plotObsPred(kx, ky, opt = opt)
        abline(0,1,lty=2, lwd=2, col='grey')
        
        if(k == 1){
          .plotLabel('a) Richness (no. present)',cex=1.2,above=T)
        }else{
          .plotLabel('b) Diversity (H)',cex=1.2,above=T)
        }
      }
      mtext(side=1, 'Observed', outer=T, line=0)
      mtext(side=2, 'Predicted', outer=T, line=0)
      
      if(!SAVEPLOTS){
        readline('no. species, effective species -- return to continue ')
      } else {
        dev.off( )
      }
    } 
  }
  
  #######################################
  
  tmp <- .omitChainCol(bgibbs,'other')
  omitBC <- tmp$omit
  keepBC <- tmp$keep
  
  ns <- min( c(ng - burnin,1000) )
  simIndex <- sample(nrow(sgibbs),ns,replace=T)
  simIndex <- sort(simIndex)
  burn <- burnin/ng*1000
  
  tmp <- .expandSigmaChains(snames, sgibbs, otherpar, simIndex, 
                            sigErrGibbs, kgibbs, REDUCT)
  corMu <- tmp$rMu; corSe <- tmp$rSe; sigMu  <- tmp$sMu; sigSe  <- tmp$sSe
  
  if(REDUCT){
    sigmaerror <- mean(sigErrGibbs)
    sinv <- .invertSigma(sigMu,sigmaerror,otherpar,REDUCT)
  } else {
    sinv <- solveRcpp(sigMu[notOther,notOther])
  }
  
  bgibbsShort    <- bgibbs[simIndex,]
  sgibbsShort    <- tmp$chainList$schain      #lower.tri with diagonal
  rgibbsShort    <- tmp$chainList$cchain
  
  if(REDUCT){
    kgibbsShort  <- tmp$chainList$kchain
    otherpar     <- output$modelList$reductList$otherpar
  }
  
  SO <- length(notOther)
  
  fMat <- output$parameters$fmatrix
  
  betaLab   <- expression( paste('Coefficient matrix ',hat(bold(B))  ))
  corLab    <- expression( paste('Correlation matrix ',hat(bold(R))  ))
  cutLab    <- expression( paste('Partition matrix ',hat(bold(plain(P)))  ))
  
  AA <- F
  if(!SMALLPLOTS)AA <- T
  
  ################if true parameter values
  
  if(TV){
    
    mfcol <- c(1,2)
    if('OC' %in% typeNames)mfcol = c(2,2)
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'trueVsPars.pdf') )
    
    par(mfcol=mfcol,bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
    colF    <- colorRampPalette(c('darkblue','orange'))
    cols    <- colF(ntypes)
    
    if('beta' %in% names(trueValues)){
      beta <- trueValues$beta
      cols <- colF(ntypes)
      if(length(beta) < 100){
        .gjamTrueVest(chains$bgibbs[,keepBC],true=beta[keepBC],
                      typeCode,allTypes,colors=cols,label = betaLab)
      } else {
        opt <- list(xlabel='true',
                    ylabel='estimate', nPerBin=length(beta)/10,
                    fill='lightblue',box.col=cols,POINTS=T,MEDIAN=F,add=F)
        .plotObsPred(beta[,notOther],betaMu[,notOther],opt = opt)
        abline(0,1,lty=2)
      }
    }
    
    if( 'corSpec' %in% names(trueValues) ){
      
      cols <- colF(2^ntypes)
      
      corTrue <- corSpec
      diag(corTrue) <- NA
      if(length(other) > 0){
        corTrue[other,] <- NA
        corTrue[,other] <- NA
      }
      
      cindex <- which(lower.tri(corSpec,diag=T))            #location on matrix
      pindex <- which(lower.tri(corSpec,diag=T),arr.ind=T)
      if(!is.matrix(pindex)){
        pindex <- matrix(pindex,1)
      }
      
      rindex <- which(is.finite(corTrue[cindex]))     #location in sgibbs
      cindex <- cindex[rindex]
      pindex <- pindex[drop=F,rindex,]   
      
      cols <- colF(ntypes + ntypes*(ntypes-1)/2)
      
      rg <- rgibbsShort
      rg[rg == 1] <- NA
      
      xlim <- range(c(-.1,.1,corTrue[cindex]),na.rm=T)
      ylim <- range(c(-.1,.1,rg),na.rm=T)
      
      add <- F
      m   <- 1
      combNames <- character(0)
      combCols  <- numeric(0)
      
      box <- F
      
      for(k in 1:length(allTypes)){
        
        wk <- which(typeNames == allTypes[k])
        wk <- wk[wk %in% notOther]
        wp <- which(pindex[,1] %in% wk & pindex[,2] %in% wk)  
        
        if( length(wp) == 1 ){
          
          combNames <- c(combNames,allTypes[k])
          
          yci <- quantile( rgibbsShort[,rindex[wp]] ,c(.5,.025,.975))
          xmu <- corSpec[matrix(pindex[wp,],1)]
          if(!add){
            plot(xmu,yci[1],xlim=xlim,ylim=ylim,
                 pch=3,col=cols[m], xlab='true',ylab='')
            add <- T
          } else {
            points(xmu,yci[1],pch=3,col=cols[m])
          }
          lines( c(xmu,xmu),yci[2:3],col=cols[m],lwd=2)
        }
        
        if(length(wp) > 1){
          
          if(length(wp) < 100){
            .gjamTrueVest(rgibbsShort[,rindex[wp]],true=corSpec[cindex[wp]],
                          typeCode,allTypes,label=corLab,xlim=xlim,ylim=ylim,
                          colors=cols[m],legend=F,add=add)
          } else {
            box <- T
            opt <- list(xlabel='true',
                        ylabel='estimate', fill='lightblue',
                        nPerBin=length(wp)/10,box.col=cols[m], POINTS=T,
                        MEDIAN=F,add=add)
            .plotObsPred(corSpec[cindex[wp]],corMu[cindex[wp]],opt = opt)
            if(!add)abline(0,1,lty=2)
          }
          add <- T
          combNames <- c(combNames,allTypes[k])
          combCols  <- c(combCols,cols[m])
          m <- m + 1
        }
        
        if(k < length(allTypes)){
          
          for( j in (k+1):length(allTypes) ){
            
            wj <- which(typeNames == allTypes[j])
            wj <- wj[wj %in% notOther]
            wp <- which(pindex[,1] %in% wk & pindex[,2] %in% wj)
            
            if(length(wp) == 0){
              wp <- which(pindex[,2] %in% wk & pindex[,1] %in% wj)
            }
            
            if(length(wp) == 0)next
            
            if(length(wp) == 1){
              yci <- quantile( rgibbsShort[,rindex[wp]] ,c(.5,.025,.975))
              xmu <- corTrue[cindex[wp]]
              if(!add){
                plot(xmu,yci[1],xlim=xlim,ylim=ylim,
                     pch=3,col=cols[m])
              } else {
                points(xmu,yci[1],pch=3,col=cols[m])
              }
              lines( c(xmu,xmu),yci[2:3],col=cols[m],lwd=2)
              
            } else {
              if(!box){
                .gjamTrueVest(rgibbsShort[,rindex[wp]],
                              true=corTrue[cindex[wp]],
                              typeCode,allTypes,add=add,colors=cols[m],
                              legend=F, xlim=c(-.9,.9), ylim=c(-.9,.9))
              } else {
                opt <- list(nPerBin=length(wp)/10,
                            box.col=cols[m], fill='white',POINTS=T,
                            MEDIAN=F,add=add)
                .plotObsPred(corSpec[cindex[wp]],corTrue[cindex[wp]],
                             opt = opt)
              }
            }
            m <- m + 1
            
            mnames    <- paste(allTypes[k],allTypes[j],sep='-')
            
            combNames <- c(combNames,mnames)
            combCols  <- c(combCols,rep(cols[m],length(mnames)))
            add <- T
          }
          
        }
      }
      legend('topleft',combNames,text.col=cols,bty='n',ncol=3,cex=.7)
    }
    
    if('OC' %in% allTypes & 'cuts' %in% names(trueValues)){
      ctmp     <- cutMu #[,-1]
      wc       <- c(1:ncol( ctmp )) + 1
      ctrue    <- cuts[,wc]
      wf       <- which(is.finite(ctrue*ctmp)[,-1])
      cutTable <- .gjamTrueVest(chains$cgibbs[,wf],true=ctrue[,-1][wf],
                                typeCode,allTypes,colors='black',
                                label=cutLab,legend=F, add=F)
    }
    
    if(!SAVEPLOTS){
      readline('simulated beta, corSpec vs betaMu, corMu (95%) -- return to continue')
    } else {
      dev.off()
    }
  }
  
  ##################### partition for ordinal
  
  if('OC' %in% typeNames){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'partition.pdf') )
    
    par( mfrow=c(1,1), bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp )
    
    wk <- which(typeNames == 'OC')
    nk <- length(wk)
    
    cgibbs <- output$chains$cgibbs
    
    
    onames <- snames[wk]
    vnames <- sort(unique(.splitNames(colnames(cgibbs))$vnam[,1]))
    
    cgibbs[!is.finite(cgibbs)] <- NA
    cc <- colSums(abs(cgibbs),na.rm=T)
    cg <- cgibbs[,cc > 0]
    
    if('cuts' %in% names(trueValues))rownames(cuts) <- rownames(cutMu)
    
    c1 <- names(cc)[cc > 0]
    
    colc <- colF(ncol(cutMu))
    
    nk <- length(vnames)
    plot(0,0,xlim=c(0,max(cg,na.rm=T)),ylim=c(1,1+nk),cex=.1,
         xlab='Unit variance scale',
         ylab=' ',yaxt='n')
    .yaxisHorizLabs(vnames,at=c(1:nk))
    
    for(k in 1:length(vnames)){
      
      x1   <- 0
      ym   <- .5
      
      wcg <- grep(vnames[k],colnames(cg))
      if(length(wcg) == 0)next
      
      tmp <- .chains2density(cg,varName=vnames[k], cut=2.5)
      xt  <- tmp$x
      yt  <- tmp$y
      yt  <- .5*yt/max(yt)
      
      yt <- yt + k
      
      for(j in 1:nrow(xt)){
        
        if('cuts' %in% names(trueValues)){
          lines( rep(cuts[vnames[k],j+2],2),c(k,k+1),lty=2,col=colc[j],lwd=3)
        }
        
        xj <- c(xt[j,],xt[j,ncol(xt)],xt[j,1])
        yj <- c(yt[j,],k,k)
        
        x2 <- which.max(yj)
        xm <- .2*x1 + .8*xj[x2]
        
        polygon(xj, yj, border=colc[j], col=.getColor(colc[j], .4), lwd=2)
        if(k == length(vnames)) text(xm,ym+k,j,col=colc[j])
        x1 <- xj[x2]
      }
    }
    .plotLabel('Partition by species',above=T)
    
    if(!SAVEPLOTS){
      readline('cuts vs cutMu -- return to continue')
    } else {
      dev.off()
    }
  }  
  
  ############################
  
  rmspeAll <- sqrt( mean( (y[,notOther] - ypredMu[,notOther])^2,na.rm=T ) )
  
  eBySpec <- sqrt( colMeans( (y[,notOther]/rowSums(y[,notOther]) - 
                                ypredMu[,notOther]/rowSums(ypredMu[,notOther],
                                                           na.rm=T))^2 ) )
  ordFit  <- order(eBySpec)
  
  score <- mean(yscore)
  
  fit <- signif( c(DIC,score,rmspeAll), 5)
  names(fit) <- c('DIC','score','rmspe')
  
  ################## predict y
  
  if(PLOTY){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'yPred.pdf') )
    
    npp <- 0
    
    for(k in 1:length(allTypes)){
      wk    <- which(typeCode == k)
      if( length(censor) > 0 ){
        ncc <- 0
        if( typeNames[wk[1]] %in% names(censor) ){
          wm   <- which(names(censor) == typeNames[wk[1]])
          #    wall <- wm
          wnot <- wk
          for(m in wm){
            wnot <- wnot[!wnot %in% censor[[m]]$columns]
            npp  <- npp + 1
          }
          if(length(wnot) > 0)npp <- npp + 1
        } else {
          ncc <- ncc + 1
        }
      } else {
        ncc <- 1
      }
      npp <- npp + ncc
    }  
    
    mfrow <- .getPlotLayout(npp)
    par( mfrow=mfrow, bty='n', omi=c(.3,.3,0,0), mar=c(3,2,2,1), 
         tcl= tcl, mgp=mgp )
    
    ylab <- ' '
    mk   <- 0
    
    ypred <- ypredMu
    yobs  <- y
    ypred[ymiss] <- yobs[ymiss] <- NA
    if(TIME)ypred[timeZero,] <- NA
    
    for(k in 1:length(allTypes)){
      
      wk    <- which(typeCode == k)
      wk    <- wk[wk %in% notOther]
      wkm   <- wk
      nk    <- nkm <- length(wk)
      censm <- NULL
      wm    <- wall <- 1
      CENS  <- F
      add   <- F
      
      if( length(censor) > 0 ){
        if( typeNames[wk[1]] %in% names(censor) ){
          CENS <- T
          wm   <- which(names(censor) == typeNames[wk[1]])
          wall <- wm
          wnot <- wk
          for(m in wm){
            wnot <- wnot[!wnot %in% censor[[m]]$columns]
          }
          if(length(wnot) > 0)wall <- c(wall,max(wall) + 1)
        }
      }
      
      for(m in wall){
        
        if(CENS){
          if(m %in% wm){
            censm <- censor[[m]]
            wkm   <- censor[[m]]$columns
          } else {
            censm <- NULL
            wkm   <- wnot
          }
          nkm <- length(wkm)
        }
        
        mk <- mk + 1
        
        y1 <- yobs[,wkm,drop=F]
        yp <- ypred[,wkm,drop=F]
        
        tmp <- .gjamPlotPars(type=typeNames[wk[1]],y1,yp,censm)
        y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
        vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
        breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
        MEDIAN <- tmp$MEDIAN
        
        SQRT <- F
        if(LOG)SQRT <- T
        
        if(typeNames[wk[1]] == 'CA')nPerBin <- NULL
        
        tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
        breaks <- tmp$breaks
        bins   <- tmp$bins
        nbin   <- tmp$nbin
        
        if(length(bins) > 0){
          breaks <- bins
          nPerBin <- NULL
        }
        
        xy <- NULL
        if(typeNames[wk[1]] == 'PA'){
          atx <- labx <- c(0,1)
          aty <- laby <- c(0,1)
        }
        
        if( !typeNames[wk[1]] %in% c('PA','CAT') ){
          ncc <- max( c(100,max(y1, na.rm=T)/20) )
          if(min(y1, na.rm=T) < bins[1])bins[1] <- min(y1, na.rm=T)
          xy  <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
          xy[2,] <- ylimit[1] + .3*xy[2,]*diff(ylimit)/max(xy[2,])
          xy[1,xy[1,] < xlimit[1]] <- xlimit[1]
          xy[2,xy[2,] < ylimit[1]] <- ylimit[1]
          
          if(SQRT){
            y1     <- sqrt(y1)
            yp     <- sqrt(yp)
            ylimit <- 1.1*sqrt(ylimit)
            xlimit <- 1.1*sqrt(xlimit)
            xy     <- sqrt(xy)
            ss     <- sqrtSeq(ylimit[2])
            aty    <- ss$at
            laby   <- ss$labs
            ss     <- sqrtSeq(xlimit[2])
            atx    <- ss$at
            labx   <- ss$labs
            plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
                 xlab='',ylab='', xaxt='n',yaxt='n')
            axis(1, at = atx, labels = labx)
            axis(2, at = aty, labels = laby)
          }else{
            if(is.null(xy)){
              plot(NULL,xlim=xlimit,ylim=ylimit,
                   xlab='',ylab='')
              
            }else{
              plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,
                   ylim=ylimit, xlab='',ylab='')
              polygon(xy[1,],xy[2,],border='tan',col='wheat')
            }
          }
          
        } else {
          y11 <- mean(y1,na.rm=T)
          y00 <- 1 - y11
          x11 <- c(-.07,-.07,.07,.07,.93,.93,1.07,1.07,-.07)
          y11 <- c(0,y00,y00,0,0,y11,y11,0,0)
          
          if(SQRT){
            y1     <- sqrt(y1)
            yp     <- sqrt(yp)
            ylimit <- 1.1*sqrt(ylimit)
            xlimit <- 1.1*sqrt(xlimit)
            xy     <- sqrt(xy)
            ss     <- sqrtSeq(ylimit[2])
            aty    <- ss$at
            laby   <- ss$labs
            ss     <- sqrtSeq(xlimit[2])
            atx    <- ss$at
            labx   <- ss$labs
          }
          plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
               xlab='Observed',ylab='Predicted', xaxt='n',yaxt='n')
          axis(1, at = atx, labels = labx)
          axis(2, at = aty, labels = laby)
          polygon(x11,y11,border='tan',col='wheat')
        }
        abline(0,1,lty=2,lwd=3,col='grey')
        
        add <- T
        
        if(nhold > 0){
          y1h <- y[holdoutIndex,wkm,drop=F]
          yph <- ypredMu[holdoutIndex,wkm,drop=F]
          points(y1h,yph,col='brown',
                 pch=21, bg='green',cex=.3)
        } 
        
        if(xlimit[2] < max(bins, na.rm=T))xlimit[2] <- max(bins, na.rm=T) + 1
        
        opt <- list(log=F, xlabel='Observed', bins = bins,
                    nbin=nbin, ylabel='Predicted', col='blue', 
                    ylimit=ylimit, xlimit = xlimit, SQRT=F, add=T)
        tmp <- .plotObsPred(y1, yp, opt = opt)
        
        if(length(vlines) > 0)abline(v=vlines,lty=2)
        
        tf <- .gjamGetTypes(typeNames[wk[1]])$labels
        tf <- paste(letters[mk],tf, sep=') ')
        
        .plotLabel(tf,'topleft',above=AA)
      }
    }
    mtext('Observed', side=1, outer=T)
    mtext('Predicted', side=2, outer=T)
    
    
    if(!SAVEPLOTS){
      readline('obs y vs predicted y -- return to continue ')
    } else {
      dev.off()
    }  
  }##########################
  
  nfact <- factorBeta$nfact
  factorList <- factorBeta$factorList
  contrast <- factorBeta$contrast
  
  if(TIME){
    nfact <- nfact + factorLambda$nfact
    factorList <- append(factorList, factorLambda$factorList)
    contrast   <- append(contrast, factorLambda$contrast)
  }
  
  if( PLOTX & PREDICTX & length(xpredMu) > 0){
    
    noX <- character(0)
    colorGrad   <- colorRampPalette(c('white','brown','black'))
    
    iy <- c(1:n)
    if(!is.null(timeZero))iy <- iy[-timeZero]
    
    if(nfact > 0){
      
      nn <- length(unlist(factorList)) # + nfact
      mmat <- matrix(0,nn,nn)
      mnames <- rep('bogus',nn)
      samples <- rep(0,nn)
      
      ib <- 1
      
      par(mfrow=c(1,1),bty='n')
      
      mm <- max(nfact,2)
      useCols <- colorRampPalette(c('brown','orange','darkblue'))(mm)
      textCol <- character(0)
      
      for(kk in 1:nfact){
        
        gname <- names( factorList )[[kk]]
        fnames <- factorList[[kk]]
        nx     <- length(fnames)
        if(nx < 1)next
        
        ie <- ib + nx - 1
        
        noX <- c(noX,fnames)
        cont <- contrast[[kk]]
        
        refClass <- names(which( rowSums( cont ) == 0) )
        
        hnames <- substring(fnames, nchar(gname) + 1)
        
        #   ff <- strsplit(fnames,gname) 
        #   hnames <- matrix( unlist(ff ),nx,2,byrow=T)[,2]
        knames <- c(paste(gname,'Ref',sep=''),fnames)
        if(TIME){
          xtrue     <- xtime[iy,fnames,drop=F]
        }else{
          xtrue <- x[iy,fnames,drop=F]
        }
        nx    <- ncol(xtrue)
        
        xpred <- xpredMu[iy,fnames,drop=F]
        cmat  <- matrix(0,nx,nx)
        colnames(cmat) <- hnames
        rownames(cmat) <- rev(hnames)
        #    wt <- apply(xtrue,1,which.max)
        for(j in 1:nx){
          wj <- which(xtrue[,j] == 1)
          cmat[,j] <- rev( colSums(xpred[drop=F,wj,],na.rm=T)/length(wj) )
        }
        
        nb <- nn - ib + 1
        ne <- nn - ie + 1
        samples[ib:ie] <- colSums(xtrue)/n
        
        mmat[ne:nb,ib:ie] <- cmat
        mnames[ib:ie] <- hnames 
        
        textCol <- c(textCol,rep(useCols[kk],nx))
        
        ib <- ie + 1
      }
      
      colnames(mmat) <- mnames
      rownames(mmat) <- rev(mnames)
      
      if(length(mmat) == 1){
        
        mc <- c(mmat[1], 1 - mmat[1])
        mmat <- cbind(rev(mc),mc)
        rownames(mmat) <- colnames(mmat) <- factorBeta$facList2[[1]]
      }
      
      graphics.off()
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'xPredFactors.pdf' ) )
      par(mfrow=c(1,1),bty='n')
      
      slim <- 1.3*c(0,max(mmat))
      if(slim[2] > 1)slim[2] <- 1
      
      .corPlot(mmat,slim=slim,plotScale=.8, textCol = textCol,
               PDIAG=F,CORLINES=T, tri='both',
               SPECLABS = T, colorGrad = colorGrad,
               textSize=1, new=F)
      if(nx > 1){
        mloc <- par('usr')
        text(mean(mloc[1:2]),mloc[3] + .03*diff(mloc[3:4]),'Observed')
        mtext('Predicted',side=4)
      }
      
      if(!SAVEPLOTS){
        readline('x inverse prediction, factors -- return to continue ')
      } else {
        dev.off()
      }
    }
    
    noplot <- c(1,grep(':',xnames),grep('^2',xnames,fixed=T))
    vnames <- xnames[-noplot]
    vnames <- vnames[!vnames %in% noX]
    
    if(length(vnames) > 0){
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'xPred.pdf') )
      
      ylab <- ''
      xlab <- ''
      mfrow <- .getPlotLayout( length(vnames) )
      par( mfrow=mfrow, bty='n', omi=c(.3,.3,0,0), mar=c(3,2,2,1),
           tcl= tcl, mgp=mgp )
      
      missX <- missingIndex
      xmaxX <- apply(x,2,max,na.rm=T)
      
      k <- 0
      b <- 0
      
      for(j in 2:Q){
        
        if(!xnames[j] %in% vnames)next
        
        k <- k + 1
        b <- b + 1
        if(b == mfrow[2])b <- 0
        
        x1 <- x[iy,j]
        x2 <- xpredMu[iy,j]
        
        type <- 'CON'
        if(length(inputs$factorBeta$factorList) > 0){
          for(kk in 1:length(inputs$factorBeta$factorList)){
            if( xnames[j] %in% inputs$factorBeta$factorList[[kk]] )type <- 'PA'
            if(all(x[,j] %in% c(0,1)))type <- 'PA'
          }
        }
        
        tmp <- .gjamPlotPars(type=type,x1,x2)
        y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
        vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
        breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
        MEDIAN <- tmp$MEDIAN
        
        LOG <- add <- F
        
        if(nhold > 0){
          x1 <- x1[-holdoutIndex]
          x2 <- x2[-holdoutIndex]
          y1 <- y1[-holdoutIndex,,drop=F]
          yp <- yp[-holdoutIndex,,drop=F]
        }
        
        log <- ''
        if(LOG)log <- 'xy'
        
        SQRT <- F
        if(LOG)SQRT <- T
        
        
        tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG, POS=F)
        breaks <- tmp$breaks
        bins   <- tmp$bins
        nbin   <- tmp$nbin
        
        if(length(bins) > 0){
          breaks <- bins
          nPerBin <- NULL
        }
        
        if(nbin > 2){
          ncc   <- max( c(100,max(y1)/20) )
          xy <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
          xy[2,] <- ylimit[1] + .3*xy[2,]*diff(ylimit)/max(xy[2,])
          plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
               xlab=' ',ylab=ylab)
          polygon(xy[1,],xy[2,],border='tan',col='wheat')
          
          
          abline(0,1,lty=2,lwd=3,col='grey')
          
          add <- T
          
          if(nhold > 0){
            points(x[holdoutIndex,j],xpredMu[holdoutIndex,j],col='brown',
                   pch=21, bg='blue',cex=.4)
          } 
        }
        
        
        opt <- list(log=F, xlabel='Observed', bins = bins,
                    nbin=nbin, ylabel='Predicted', col='darkblue', 
                    ylimit=ylimit, xlimit = xlimit, SQRT=F, add=T)
        tmp <- .plotObsPred(y1, yp, opt = opt)
        
        
        
        if(nhold > 0)points(x[holdoutIndex,j],xpredMu[holdoutIndex,j],
                            col='brown',cex=.3)
        
        if(length(missX) > 0){
          ww <- which(missX[,2] == j)
          if(length(ww) > 0){
            wz <- missX[ww,]
            if(!is.matrix(wz))wz <- matrix(wz,1)
            points(jitter(ww*0+xmaxX[j]),xpredMu[wz],cex=.6,col='blue')
          }
        }
        
        .plotLabel(paste(letters[j-1],xnames[j],sep=') '), above=AA)
      }
      mtext('Observed',side=1, outer=T)
      mtext('Predicted',side=2,outer=T)
      
      if(!SAVEPLOTS){
        readline('x inverse prediction, covariates -- return to continue ')
      } else {
        dev.off()
      }
    }
  }   
  ######################
  
  if(PLOTALLY){
    
    np <- S <- ncol(y)
    npage <- 1
    o   <- 1:S
    if(S > 16){
      np    <- 16
      npage <- ceiling(S/16)
    }
    
    mfrow <- .getPlotLayout(np)
    
    k   <- 0
    add <- F
    
    o <- 1:S
    o <- o[o <= 16]
    
    for(p in 1:npage){
      
      file <- paste('yPredBySpec_',p,'.pdf',sep='')
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
      
      par(mfrow=mfrow, bty='n', omi=c(.3,.3,0,0), mar=c(3,2,2,1), 
          tcl= tcl, mgp=mgp)
      
      for(j in o){
        
        censm <- NULL
        if( length(censor) > 0 ){
          if( typeNames[j] %in% names(censor) ){
            wjc <- which(names(censor) == typeNames[j])
            if(j %in% censor[[wjc]]$columns)censm <- censor[[wjc]]
          }
        }
        
        y1 <- y[,j]
        if(min(y1) == max(y1))next
        y2 <- ypredMu[,j]
        
        tmp <- .gjamPlotPars(type=typeNames[j],y1,y2,censm)
        y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
        vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
        breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
        MEDIAN <- tmp$MEDIAN
        
        SQRT <- F
        if(LOG)SQRT <- T
        
        tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
        breaks <- tmp$breaks
        bins   <- tmp$bins
        nbin   <- tmp$nbin
        
        if(length(bins) > 0){
          breaks <- bins
          nPerBin <- NULL
        }
        
        if( !typeNames[wk[1]] %in% c('PA','CAT') ){
          ncc   <- max( c(100,max(y1)/20) )
          if(bins[1] > min(y1))bins <- c(min(y1),bins)
          ymm <- max(y1) + diff(range(y1,na.rm=T))*.01
          bins <- c(bins[bins < ymm], ymm) 
          
          xy <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
          xy[2,] <- ylimit[1] + .8*xy[2,]*diff(ylimit)/max(xy[2,])
          
          if(SQRT){
            y1     <- sqrt(y1)
            yp     <- sqrt(yp)
            ylimit <- 1.1*sqrt(ylimit)
            xlimit <- 1.1*sqrt(xlimit)
            xy     <- sqrt(xy)
            ss     <- sqrtSeq(ylimit[2])
            aty    <- ss$at
            laby   <- ss$labs
            ss     <- sqrtSeq(xlimit[2])
            atx    <- ss$at
            labx   <- ss$labs
          }
          plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
               xlab='Observed',ylab='Predicted', xaxt='n',yaxt='n')
          axis(1, at = atx, labels = labx)
          axis(2, at = aty, labels = laby)
          polygon(xy[1,],xy[2,],border='tan',col='wheat')
          
        } else {
          y11 <- mean(y1)
          y00 <- 1 - y11
          x11 <- c(-.07,-.07,.07,.07,.93,.93,1.07,1.07,-.07)
          y11 <- c(0,y00,y00,0,0,y11,y11,0,0)
          plot(x11,y11,col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
               xlab=' ',ylab=ylab)
          polygon(x11,y11,border='tan',col='wheat')
        }
        abline(0,1,lty=2,lwd=3,col='grey')
        add <- T
        
        if(nhold > 0){
          points(y1[holdoutIndex],yp[holdoutIndex],col='brown',
                 pch=21, bg='blue',cex=.4)
        } 
        
        fill <- .getColor('blue',.3)
        
        opt <- list(log=F, xlabel='Observed', bins = bins,
                    nbin=nbin, ylabel='Predicted', col='darkblue', 
                    add=T)
        
        tmp <- .plotObsPred(y1,yp,opt = opt)
        
        if(length(vlines) > 0)abline(v=vlines,lty=2)
        
        k <- k + 1
        if(k > 26)k <- 1
        
        lab <- paste(letters[k],') ',colnames(y)[j],' - ', 
                     typeNames[j], sep='')
        
        .plotLabel( lab,above=T )
        abline(0,1,lty=2)
        abline(h = mean(y2),lty=2)
      }
      mtext('Observed', 1, outer=T)
      mtext('Predicted', 2, outer=T)
      
      
      if(!SAVEPLOTS){
        readline('y prediction -- return to continue ')
      } else {
        dev.off()
      }
      o <- o + 16
      o <- o[o <= S]
    }
    
  }
  
  ############## traits
  
  if(TRAITS){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'traitPred.pdf') ) # start plot
    
    tt <- grep('other',colnames(plotByTrait))
    if(length(tt) > 0)colnames(plotByTrait)[tt] <- colnames(specByTrait)[tt]
    
    print(colnames(plotByTrait))
    
    yy <- plotByTrait
    o  <- 1:ncol(yy)
    
    if(ncol(yy) > 16){
      
      rmspe <- sqrt( colSums( (plotByTrait - tMu)^2 )/n )
      o <- order(rmspe)[1:16]
      yy <- plotByTrait[,o]
    }
    
    mfrow <- .getPlotLayout(length(o))
    par(mfrow=mfrow, bty='n', oma=oma, mar=c(3,3,1,1), tcl= tcl, mgp=mgp)
    k <- 0
    
    for(j in o){
      
      add   <- F
      jname <- colnames(tMu)[j]
      
      k <- k + 1
      
      td <- plotByTrait[,jname]
      
      tjj   <- tMu[,j]
      wj <- which(colnames(tMuOrd) == jname)
      
      tmp <- .gjamPlotPars(type=traitTypes[j],td,tjj)
      y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
      vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
      breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
      MEDIAN <- tmp$MEDIAN
      
      if(nhold > 0){
        add <- T
        log <- ''
        if(LOG)log <- 'xy'
        plot(td[holdoutIndex],tjj[holdoutIndex],xlab=' ',ylab=ylab,
             xlim=xlimit,ylim=ylimit,col='grey',pch=21,bg='brown',cex=.4,log=log)
      } 
      
      opt <- list( xlabel=' ',ylabel=ylab,nbin=nbin, 
                   nPerBin=nPerBin,
                   xlimit=xlimit,ylimit=ylimit,breaks=breaks,
                   wide=wide,LOG=LOG,
                   fill='grey',
                   POINTS=F,MEDIAN=MEDIAN,add=add )
      
      tmp <- .plotObsPred(td, tjj, opt = opt)
      if(length(vlines) > 0)abline(v=vlines,lty=2)
      abline(0,1,lty=2)
      abline(h=mean(td,na.rm=T),lty=2)
      
      .plotLabel( paste(letters[k],') ',.traitLabel(jname),sep=''),above=AA )
    }
    
    if(!SAVEPLOTS){
      readline('predictive trait distributions -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  ##############sensitivity 
  
  nfact <- factorBeta$nfact
  if(!is.matrix(fSensGibbs)){
    fSensGibbs <- matrix(fSensGibbs)
    colnames(fSensGibbs) <- xnames[-1]
  }
  
  wc <- c(1:ncol(fSensGibbs))
  wx <- grep(':',colnames(fSensGibbs))
  wx <- c(wx, grep('^2',colnames(fSensGibbs), fixed=T) )
  if(length(wx) > 0)wc <- wc[-wx]
  
  wx <- grep('intercept',colnames(fSensGibbs))
  if(length(wx) > 0)wc <- wc[-wx]
  wc <- c(1:ncol(fSensGibbs))
  
  tmp <- apply(fSensGibbs,2,range)
  wx <- which(tmp[1,] == tmp[2,])
  if(length(wx) > 0)wc <- wc[-wx]
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'sensitivity.pdf') ) # start plot
  
  xx   <- fSensGibbs[,wc,drop=F]
  tcol <- rep('black',ncol(xx))
  names(tcol) <- colnames(xx)
  
  if(nfact > 0){
    
    mm <- max(nfact,2)
    useCols <- colorRampPalette(c('brown','orange','darkblue'))(mm)
    
    for(i in 1:nfact){
      im <- which(colnames(xx) %in% rownames(factorBeta$contrast[[i]]))
      tcol[im] <- useCols[i]
    }
  }
  
  par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
  if(TIME)par(mfrow=c(1,2),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
  
  ord  <- order( colMeans(xx) )
  ylim <- c(min(xx),1.5*quantile(xx,.95))
  tmp <- .boxplotQuant( xx[,ord, drop=F], xaxt='n',outline=F, 
                        border=tcol[ord], whiskcol=tcol[ord],
                        boxfill=.getColor(tcol[ord],.4), 
                        pars = list(boxwex = 0.5, ylim=ylim), lty=1, log='y')
  mtext('Predictors in X',side=1,line=1)
  abline(h=0,lwd=2,col='grey')
  
  dy <- .05*diff(par()$yaxp[1:2])
  text(1:length(wc), dy + tmp$stats[5,],tmp$names,srt=90,pos=4,col=tcol[ord])
  sensLab   <- expression( paste('Sensitivity ',hat(bold(F))  ))
  .plotLabel(sensLab,'bottomleft',above=F, cex=1.1)   
  
  if(TIME){
    
    tiny <- 1e-6
    xg <- chains$gsens
    xg[is.na(xg)] <- 0
    w0 <- which(colSums(xg) == 0)
    if(length(w0) > 0)xg <- xg[,-w0,drop=F]
    
    if(length(w0) > 0){
      
      tcol <- rep('black',ncol(xg))
      names(tcol) <- colnames(xg)
      
      if(factorLambda$nfact > 0){
        
        mm <- max(nfact,2)
        useCols <- colorRampPalette(c('brown','orange','darkblue'))(mm)
        
        for(i in 1:factorLambda$nfact){
          im <- which(colnames(xg) %in% rownames(factorLambda$contrast[[i]]))
          tcol[im] <- useCols[i]
        }
      }
      xm <- colMeans(xg)
      ord  <- order( xm )
      
      ylim <- c(min(xg),2*quantile(xg,.9999,na.rm=T))
      if(ylim[1] < 1e-8)ylim[1] <- 1e-8
      tmp <- .boxplotQuant( xg[,ord, drop=F], xaxt='n',outline=F, 
                            border=tcol[ord],whiskcol=tcol[ord],
                            boxfill=.getColor(tcol[ord],.4), 
                            pars = list(boxwex = 0.5, ylim=ylim), lty=1, log='y')
      mtext('Predictors in V',side=1,line=1)
      abline(h=0,lwd=2,col='grey')
      dy <- .05*diff(par()$yaxp[1:2])
      text(1:length(ord), dy + tmp$stats[5,],tmp$names,srt=90,pos=4,col=tcol[ord])
      sensLab   <- expression( paste('Sensitivity ',hat(bold(lambda))  ))
      .plotLabel(sensLab,'bottomright',above=F, cex=1.1)  
      
    }
    
    if(!SAVEPLOTS){
      readline('sensitivity over full model -- return to continue ')
    } else { 
      dev.off()
    }
  }
  
  if(TIME){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'sensitivityAlpha.pdf') ) 
    par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
    
    tiny <- 1e-6
    xg <- chains$asens
    #   xg[xg < tiny] <- tiny
    xm <- colMeans(xg)
    ord  <- order( xm, decreasing=T )
    
    wo <- 50       # largest values
    if(wo > ncol(xg))wo <- ncol(xg)
    xc <- rev(ord[1:wo])
    
    ylim <- c(min(xg),2*quantile(xg,.9999))
    if(ylim[1] < 1e-8)ylim[1] <- 1e-8
    tmp <- .boxplotQuant( xg[,xc, drop=F], xaxt='n',outline=F,  
                          pars = list(boxwex = 0.5, ylim=NULL), lty=1, log='y')
    mtext('Predictors in U',side=1,line=1)
    abline(h=0,lwd=2,col='grey')
    dy <- .05*diff(par()$yaxp[1:2])
    text(1:wo, dy + tmp$stats[5,],tmp$names,srt=90,pos=4)
    sensLab   <- expression( paste('Sensitivity ',hat(bold(alpha))  ))
    .plotLabel(sensLab,'bottomright',above=F, cex=1.1)  
    
    if(!SAVEPLOTS){
      readline('sensitivity over species pairs -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  ######################  coefficient summary tables ############
  
  fnames <- rownames(factorBeta$eCont)
  
  # bTab   <- .getSigTable(bgibbs,S, Q, xnames, snames) 
  
  # q1    <- nrow(factorBeta$eCont)
  # 
  # bfTab <- .getSigTable(bFacGibbs,SO, q1, fnames, 
  #                       colnames(parameters$fBetaMu)) 
  
  # bfCoeffTable <- .processPars(bFacGibbs,sigOnly=SIGONLY)$summary
  # sigFbeta     <- rownames(bfCoeffTable)
  
  # bfSig <- bFacGibbs[,sigFbeta]
  
  # bCoeffTable <- .processPars(bgibbs[,keepBC],sigOnly=SIGONLY)$summary
  # sigBeta     <- rownames(bCoeffTable)
  # bCoeffTable <- .processPars(bgibbs[,keepBC],sigOnly=F)$summary
  
  # if(length(sigBeta) == 0)sigBeta <- c(1:ncol(bgibbs))
  
  # scaleNote <- 'W/X scale'
  
  # betaSig <- bgibbs[,sigBeta]
  
  # summaryCoeffs <- list(betaSig = bTab, fBetaSig = bfTab, 
  #                       betaCoeff = bCoeffTable, fBetaCoeff = bfCoeffTable)
  ##################################333333333
  
  
  
  tmp <- .splitNames(colnames(bgibbs),snames=colnames(y))
  vnames <- unique(tmp$vnam)
  xnam <- unique(tmp$xnam[tmp$xnam != 'intercept'])
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'betaChains.pdf') ) # start plot
  
  if(CHAINS){
    cseq <- 1:nrow(bgibbs)
    if(length(cseq) > 1000)cseq <- seq(1,length(cseq),length=1000)
    
    
    mfrow <- .getPlotLayout(length(xnam))
    par(mfrow=mfrow, bty='n', oma=oma, mar=c(2,2,1,1), tcl= tcl, mgp=mgp)
    
    flist <- factorBeta$factorList
    if(length(flist) > 0){
      flist <- sort(unique(unlist(flist)))
    }
    
    for(k in 1:length(xnam)){
      
      tname <- xnam[k]
      tmp   <- .chains2density(bgibbs[cseq,],varName=tname, cut=3)
      
      xt  <- tmp$x
      yt  <- tmp$y
      chainMat <- tmp$chainMat
      
      if(ncol(chainMat) > 20)chainMat <- chainMat[,sample(ncol(chainMat),20)]
      
      colF <- colorRampPalette(c('darkblue','orange'))
      cols <- colF(nrow(xt))
      
      snamek <- .splitNames(colnames(chainMat),colnames(y))$vnam
      
      nn <- nrow(chainMat)
      
      jk <- 1:ncol(chainMat)
      if(length(jk) > 20)jk <- sample(jk,20)
      plot(0,0,xlim=c(0,(1.4*nn)),ylim=range(chainMat[,jk]),
           xlab=' ',ylab=' ',cex=.01)
      
      for(j in jk){
        lines(chainMat[,j],col=cols[j])
        if(ncol(chainMat) < 15)text(nn,chainMat[nn,j],snamek[j],col=cols[j],pos=4)
        abline(v=burn,lty=2)
        
        if(k == 1 & j == 1).plotLabel( paste(burnin,":",ng),
                                       location='topright' )
      }
      .plotLabel(label=paste(letters[k],') ',tname,sep=''),
                 location='topleft',above=T)
      
      abline(h=0,lwd=4,col='white')
      abline(h=0,lty=2)
      
      if(ncol(chainMat) >= 15) text(nn,mean(par('usr')[3:4]),
                                    paste(ncol(chainMat),'spp'),pos=4)
    }
    
    if(!SAVEPLOTS){
      readline('beta coefficient thinned chains -- return to continue ')
    } else {
      dev.off()
    }
  }
  ######################### correlation chains, species at random
  
  if(CHAINS){
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'corChains.pdf') ) # start plot
    
    par(mfrow=c(2,2), bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
    
    w0 <- 1:ncol(sgibbs)
    if(REDUCT){
      same <- .sameByColumn(kgibbs)
      w0   <- order(same[lower.tri(same,diag=T)])
    } else {
      w0 <- sample(max(w0),80,replace=T)
    }
    ww <- 1:20
    
    for(jj in 1:4){
      
      ssj <- w0[ww]
      ssj <- ssj[is.finite(ssj)]
      if(length(ssj) == 0)break
      if(max(ssj) > max(w0))break
      ww  <- ww + 20
      
      tmp   <- .chains2density(rgibbsShort[,ssj])
      xt    <- tmp$x
      yt    <- tmp$y
      chainMat <- tmp$chainMat
      
      colF <- colorRampPalette(c('black','brown','orange'))
      cols <- colF(nrow(xt))
      stk  <- .splitNames(colnames(chainMat))$vnam
      
      ws <- which(stk[,1] == stk[,2])
      if(length(ws) > 0){
        stk <- stk[-ws,]
        chainMat <- chainMat[,-ws]
      }
      
      rr <- range(chainMat)
      if(!is.finite(rr[1]) | !is.finite(rr[2]))next
      
      if(is.matrix(chainMat)){
        
        snamek <- stk[,1]
        nn <- nrow(chainMat)
        plot(0,0,xlim=c(0,(1.4*nn)),ylim=range(chainMat),xlab=' ',ylab=' ',cex=.01)
        
        jk <- 1:ncol(chainMat)
        if(length(jk) > 20)jk <- sample(jk,20)
        
        for(j in jk){
          lines(chainMat[,j],col=cols[j])
          if(ncol(chainMat) < 15)text(nn,chainMat[nn,j],snamek[j],col=cols[j],pos=4)
          
        }
        
        if(jj == 1).plotLabel( paste(burnin,":",ng),location='topright' )
        abline(h=0,lwd=4,col='white')
        abline(h=0,lty=2)
        abline(v=burn,lty=2)
        
        if(ncol(chainMat) >= 15) text(nn,mean(par('usr')[3:4]),
                                      paste(ncol(chainMat),'spp'),pos=4)
      }
    }
    
    if(!SAVEPLOTS){
      readline('correlation thinned chains -- return to continue ')
    } else {
      dev.off()
    }
  }
  ##################### time chains
  
  if(TIME & CHAINS){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'lambdaChains.pdf') ) 
    
    tmp <- .splitNames(colnames(ggibbs),colnames(y))
    vnames <- unique(tmp$vnam)
    xnam <- unique(tmp$xnam)
    
    cseq <- 1:nrow(ggibbs)
    if(length(cseq) > 1000)cseq <- seq(1,length(cseq),length=1000)
    
    mfrow <- .getPlotLayout(length(xnam))
    par(mfrow=mfrow, bty='n', oma=oma, mar=c(2,2,1,1), tcl= tcl, mgp=mgp)
    
    for(k in 1:length(xnam)){
      
      tname <- xnam[k]
      
      tmp <- .chains2density(ggibbs[cseq,],varName=tname, cut=3)
      xt  <- tmp$x
      yt  <- tmp$y
      chainMat <- tmp$chainMat
      
      if(ncol(chainMat) > 20)chainMat <- chainMat[,sample(ncol(chainMat),20)]
      
      colF <- colorRampPalette(c('darkblue','orange'))
      cols <- colF(nrow(xt))
      
      snamek <- .splitNames(colnames(chainMat),colnames(y))$vnam
      
      nn <- nrow(chainMat)
      
      jk <- 1:ncol(chainMat)
      if(length(jk) > 20)jk <- sample(jk,20)
      plot(0,0,xlim=c(0,(1.4*nn)),ylim=range(chainMat[,jk]),
           xlab=' ',ylab=' ',cex=.01)
      
      for(j in jk){
        lines(chainMat[,j],col=cols[j])
        if(ncol(chainMat) < 15)text(nn,chainMat[nn,j],snamek[j],col=cols[j],pos=4)
        if(k == 1 & j == 1).plotLabel( paste('burn-in =',burnin),
                                       location='topright' )
      }
      if(k == 1)tname <- character(0)
      lab <- paste('lambda',tname)
      .plotLabel(label=paste(letters[k],') ',lab,sep=''),location='topleft',above=T)
      
      abline(h=0,lwd=4,col='white')
      abline(h=0,lty=2)
      
      if(ncol(chainMat) >= 15) text(nn,mean(par('usr')[3:4]),
                                    paste(ncol(chainMat),'spp'),pos=4)
    }
    
    if(!SAVEPLOTS){
      readline('lambda coefficient chains -- return to continue ')
    } else {
      dev.off()
    }
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'alphaChains.pdf') ) 
    
    cseq <- 1:nrow(alphaGibbs)
    if(length(cseq) > 1000)cseq <- seq(1,length(cseq),length=1000)
    
    np <- min(c(S,4))
    mfrow <- .getPlotLayout(np)
    par(mfrow=mfrow, bty='n', oma=oma, mar=c(2,2,1,1), tcl= tcl, mgp=mgp)
    
    kp <- min(c( 4, floor(S/4) ) )
    ka <- c(1:S)
    
    for(k in 1:np){
      
      wc <- sample(ka,kp)
      ka <- ka[!ka %in% wc]
      
      tmp <- .chains2density(alphaGibbs[cseq,wc], cut=3)
      xt  <- tmp$x
      yt  <- tmp$y
      chainMat <- tmp$chainMat
      
      colF <- colorRampPalette(c('darkblue','orange'))
      cols <- colF(nrow(xt))
      
      snamek <- .splitNames(colnames(chainMat),colnames(y))$vnam
      
      nn <- nrow(chainMat)
      
      jk <- 1:ncol(chainMat)
      if(length(jk) > 20)jk <- sample(jk,20)
      plot(0,0,xlim=c(0,(1.4*nn)),ylim=range(chainMat[,jk]),
           xlab=' ',ylab=' ',cex=.01)
      
      for(j in jk){
        lines(chainMat[,j],col=cols[j])
        if(ncol(chainMat) < 15)text(nn,chainMat[nn,j],snamek[j],
                                    col=cols[j],pos=4)
        if(k == 1 & j == 1).plotLabel( paste('burn-in =',burnin),
                                       location='topright' )
      }
      abline(h=0,lwd=4,col='white')
      abline(h=0,lty=2)
      
      if(ncol(chainMat) >= 15) text(nn,mean(par('usr')[3:4]),
                                    paste(ncol(chainMat),'spp'),pos=4)
    }
    
    if(!SAVEPLOTS){
      readline('alpha coefficient chains -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  ############################### beta posteriors as boxes
  
  fMu <- parameters$betaStandXWTable
  
  sigFbeta <- rownames(fMu)[fMu$sig95 == '*']
  bfSig <- bFacGibbs[,sigFbeta]
  
  if(length(bfSig) > 0){
    
    tmp <- .splitNames(colnames(bfSig), snames)
    vnam <- tmp$vnam
    xnam <- tmp$xnam
    
    xpNames <- .replaceString(fnames,':','X')
    xpNames <- .replaceString(xpNames,'I(','')
    xpNames <- .replaceString(xpNames,')','')
    xpNames <- .replaceString(xpNames,'^2','2')
    xpNames <- .replaceString(xpNames,'*','TIMES')
    
    fnames <- unique( xnam )
    
    brange <- apply(bfSig,2,range)
    
    for(j in 1:length(fnames)){
      
      wc <- which(xnam == fnames[j] & brange[2,] > brange[1,])
      if(length(wc) < 2)next
      
      plab <- paste('beta_',xpNames[j],'.pdf',sep='')
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,plab) ) # start plot
      
      par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
      
      .myBoxPlot( mat = bfSig[,wc], tnam = vnam[ wc ], snames = snames,
                  specColor, label=fnames[j], LEG=T)
      mtext(side=2,'Coefficient', line=2)
      
      
      if(!SAVEPLOTS){
        readline('standardized for W/X, 95% posterior -- return to continue ')
      } else {
        dev.off()
      }
    }
    
    #one figure
    
    if(length(fnames) > 1){
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'betaAll.pdf') )  
      
      npp <- length(which(table(match(xnam,fnames)) > 1))
      
      mfrow <- .getPlotLayout(npp)
      par( mfrow=mfrow, bty='n', omi=c(.3,.5,0,0), 
           mar=c(1,1,1,1), tcl= tcl )
      
      k <- 0
      for(j in 1:length(fnames)){
        
        wc <- which(xnam == fnames[j])
        if(length(wc) < 2)next
        
        k <- k + 1
        
        .myBoxPlot( mat = bfSig[,wc], tnam = vnam[ wc ], snames = snames,
                    specColor, label=' ', LEG=F)
        .plotLabel(fnames[j],'bottomleft')
      }
      mtext(side=2,'Coefficient value',outer=T, line=1)
      
      if(!SAVEPLOTS){
        readline('95% posterior -- return to continue ')
      } else {
        dev.off()
      }
    }
  }
  
  ############################## time #######################
  if(TIME){
    
    ggibbs <- chains$ggibbs  #lambda
    
    tmp  <- .splitNames(colnames(chains$ggibbs), snames)
    vnam <- tmp$vnam
    xnam <- tmp$xnam 
    gnames <- unique(xnam)
    
    k <- 0
    
    for(j in 1:length(gnames)){
      
      wc <- which(xnam == gnames[j])
      if(length(wc) < 2)next
      
      k <- k + 1
      
      plab <- paste('lambda_',gnames[j],'.pdf',sep='')
      if(j == 1){
        glab <- 'lambda'
      }else{
        glab <- paste('lambda:',gnames[j])
      }
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,plab) ) # start plot
      par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
      
      .myBoxPlot( mat = ggibbs[,wc], tnam = vnam[ wc ], snames = snames,
                  specColor, label=glab)
      if(j == 1)abline(h=1, col=.getColor('black',.3), lwd=2, lty=2)
      if(!SAVEPLOTS){
        readline('95% posterior -- return to continue ')
      } else {
        dev.off()
      }
    }
    
    # one plot
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'lambdaAll.pdf') )  
    
    npp <- length(which(table(match(xnam,gnames)) > 1))
    mfrow <- .getPlotLayout(npp)
    par( mfrow=mfrow, bty='n', oma=oma, mar=c(1,1,1,1), tcl= tcl, mgp=mgp )
    
    k <- 0
    for(j in 1:length(gnames)){
      
      wc <- which(xnam == gnames[j])
      if(length(wc) < 2)next
      
      k <- k + 1
      if(j == 1){
        glab <- 'lambda'
      }else{
        glab <- paste('lambda:',gnames[j])
      }
      .myBoxPlot( mat = ggibbs[,wc], tnam = vnam[ wc ], snames = snames,
                  specColor, label=glab)
      if(j == 1)abline(h=1, col=.getColor('black',.3), lwd=2, lty=2)
    }
    
    if(!SAVEPLOTS){
      readline('95% posterior -- return to continue ')
    } else {
      dev.off()
    }
  }  ### end time ##
  
  ############################### beta posteriors, traits
  
  if(TRAITS){
    
    M  <- nrow(specByTrait)
    nc     <- 0
    vnam   <- .splitNames(colnames(chains$bTraitFacGibbs))$vnam
    mnames <- colnames(specByTrait)
    
    if( length(is.finite(match(mnames,vnam[,1]))) > 0 )nc <- 2
    if( length(is.finite(match(mnames,vnam[,2]))) > 0 )nc <- 1
    
    ix <- 1
    if(nc == 1)ix <- 2
    xnam <- vnam[,ix]
    vnam <- vnam[,nc]
    
    if(length(traitColor) == 1)traitColor <- rep(traitColor, M)
    tboxCol <- .getColor(traitColor,.4)
    
    traitSd <- apply(plotByTrait,2,sd,na.rm=T)
    traitSd <- matrix(traitSd,nrow(chains$bTraitGibbs),length(traitSd),byrow=T)
    
    for(j in 2:length(xnames)){
      
      wc <- which(xnam == xnames[j])
      if(length(wc) < 2)next
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'traits.pdf') ) # start plot
      
      par(mfrow=c(1,1),bty='n', oma=oma, mar=mar, tcl= tcl, mgp=mgp)
      
      if(length(wc) > 100)wc <- sample(wc,100)
      
      mat <- chains$bTraitGibbs[,wc]*xSd[j]/traitSd
      vn  <- .splitNames(colnames(mat))$vnam[,1]
      
      .myBoxPlot( mat, tnam = vn, snames = mnames,
                  traitColor, label=' ', LEG=T)
      
      .plotLabel(xnames[j],location='bottomright')  
      
      if(!SAVEPLOTS){
        readline('traits, standardized for X/W, 95% posterior -- return to continue ')
      } else {
        dev.off()
      }
    }
  }
  
  ########### cluster analysis
  
  covx <- cov(x)
  covy <- cov(y[,notOmit])
  
  wo <- which(whichZero[,1] %in% other | whichZero[,2] %in% other)
  if(length(wo) > 0)whichZero <- whichZero[-wo,]
  wo <- which(whConZero[,1] %in% other | whConZero[,2] %in% other)
  if(length(wo) > 0)whConZero <- whConZero[-wo,]
  
  nsim <- 500
  if(S > 50)nsim  <- 100
  if(S > 100)nsim <- 20
  
  tmp <- eigen( ematrix[notOther,notOther] )
  
  eVecs   <- tmp$vectors
  eValues <- tmp$values
  rownames(eVecs) <- snames[notOther]
  
  if(!GRIDPLOTS){
    
    clusterIndex <- NULL
    clusterOrder <- NULL
    
    if(S >= 8){
      opt <- list( ncluster=ncluster, PLOT=F, DIST=F )
      clusterDat <- .clusterPlot( ematrix , opt)
      colCode    <- clusterDat$colCode
      cord       <- rev(clusterDat$corder)
      dord       <- notOther[!notOther %in% omit][cord]
      
      clusterIndex <- clusterDat$clusterIndex
      clusterOrder <- clusterDat$corder
    }
    
    invisible( return( list(fit = fit, 
                            ematrix = ematrix, clusterIndex = clusterIndex, 
                            clusterOrder = clusterOrder) ) )
  }
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'clusterDataE.pdf') ) # start plot
  
  mag <- mar
  mag[4] <- max(mar[4],6)
  par(mfrow=c(1,2), cex=.7, oma=oma, mar=mag, tcl= tcl, mgp=mgp)
  
  LABELS <- T
  if(S > 100 | !SPECLABS)LABELS <- F
  
  dcor <- .cov2Cor(covy)
  dcor[is.na(dcor)] <- 0
  
  opt <- list( main='',cex=.2,ncluster=ncluster,
               colCode=specColor[notOmit], textSize=.4, 
               LABELS = LABELS, DIST=F )
  tmp <- .clusterPlot( dcor, opt)
  colCode <- tmp$colCode
  
  clusterIndex <- tmp$clusterIndex
  clusterOrder <- tmp$corder
  
  .plotLabel('a) Data correlation',above=T, cex=1.7)
  
  
  tmp   <- .clustMat(ematrix[notOther,notOther], SYM = T)
  ecor <- tmp$cmat
  
  
  
  opt <- list( main='',cex=.2, ncluster=ncluster, 
               colCode=specColor[notOmit], textSize=.5, 
               LABELS = LABELS, DIST=F)
  tmp <- .clusterPlot( ecor , opt )
  .plotLabel('b) E correlation',above=T, cex=1.7)
  
  clusterIndex <- cbind( clusterIndex, tmp$clusterIndex )
  clusterOrder <- cbind( clusterOrder, tmp$corder )
  
  rownames(clusterIndex) <- rownames(clusterOrder) <- snames[notOmit]
  colnames(clusterIndex) <- colnames(clusterOrder) <- c('data','E')
  
  if(!SAVEPLOTS){
    readline('Data and E responses to X -- return to continue ')
  } else {
    dev.off()
  }
  
  ########### E communities
  
  imat <- output$inputs$y
  imat[imat > 0] <- 1
  iord <- colSums(imat)
  etab  <- table(clusterIndex[,'E'])
  eComs <- matrix(NA,ncluster, max(etab))
  ename <- rep( character(0), max(etab) )
  
  egroup <- clusterIndex[,'E']
  # bTab   <- cbind(egroup,bTab[notOther,])
  # summaryCoeffs$betaSig <- bTab
  
  # bfTab <- cbind(egroup, bfTab[notOther,])
  # summaryCoeffs$fBetaSig <- bfTab
  
  for(j in 1:ncluster){
    
    wj <- which(clusterIndex[,'E'] == j)
    jname <- rownames(clusterIndex)[wj]
    jname <- jname[order(iord[jname],decreasing=T)]
    eComs[j,1:length(jname)] <- jname
    mm    <- min( c(3,length(jname)) )
    jj    <- substr(jname[1:mm],1,6)
    ename[j] <- paste0(jj,collapse='_')
  } 
  rownames(eComs) <- ename
  eComs <- t(eComs)
  
  ########### ordination
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'ordination.pdf') ) # start plot
  
  clusNames <- eComs[1,]
  
  lambda <- eValues/sum(eValues)
  cl     <- cumsum(lambda)
  
  cbord <- .getColor(specColor[notOther],.4)
  cfill <- .getColor(specColor[notOther],.4)
  
  par(mfcol=c(2,2), bty='n', cex = cex, mar=c(4,4,1,1))
  
  p1 <- paste('Axis I (',round(100*lambda[1],0),'%)',sep='')
  p2 <- paste('Axis II (',round(100*lambda[2],0),'%)',sep='')
  p3 <- paste('Axis III (',round(100*lambda[3],0),'%)',sep='')
  
  xlim <- range(eVecs[,1])
  
  plot(eVecs[,1],eVecs[,2],cex=1,col=cbord, bg = cfill, pch=16,
       xlab=p1, ylab = p2) 
  abline(h=0,col=.getColor('black',.1),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.1),lwd=2,lty=2)
  
  text(eVecs[clusNames,1],eVecs[clusNames,2],substr(clusNames,1,7))
  
  plot(eVecs[,1],eVecs[,3],cex=1,col=cbord, bg = cfill, pch=16,
       xlab=p1, ylab = p3) 
  abline(h=0,col=.getColor('black',.1),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.1),lwd=2,lty=2)
  
  text(eVecs[clusNames,1],eVecs[clusNames,3],substr(clusNames,1,7))
  
  plot(eVecs[,2],eVecs[,3],cex=1,col=cbord, bg = cfill, pch=16,
       xlab=p2, ylab = p3)
  abline(h=0,col=.getColor('black',.1),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.1),lwd=2,lty=2)
  
  text(eVecs[clusNames,2],eVecs[clusNames,3],substr(clusNames,1,7))
  
  plot(cl,type='s',xlab='Rank',ylab='Proportion of variance',xlim=c(.9,S),
       ylim=c(0,1),log='x')
  lines(c(.9,1),c(0,cl[1]),lwd=2,type='s')
  for(j in 1:length(lambda))lines(c(j,j),c(0,cl[j]),col='grey')
  lines(cl,lwd=2,type='s')
  abline(h=1,lwd=2,col=.getColor('grey',.5),lty=2)
  
  if(!SAVEPLOTS){
    readline('ordination of E matrix -- return to continue ')
  } else {
    dev.off()
  }
  
  ########### dimension reduction ############
  
  if(REDUCT){ 
    
    graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'dimRed.pdf') ) # start plot
    
    mk <- .modalValuesInArray(kgibbs,2)[notOmit]
    NK <- table( table(mk) )
    mk <- length(NK)
    
    r <- otherpar$r
    
    par(bty='n')
    scale <- SO/3
    if(SMALLPLOTS)scale <- 10*scale
    .mapSetup(c(1,SO),c(1,SO),scale=scale)
    
    xl <- SO/15
    yl <- SO/8
    
    en <- SO*(SO+1)/2
    
    plot(0,0,xlim=c(0,SO+xl),ylim=c(0,SO+xl),cex=.01,xaxt='n',yaxt='n',
         xlab=' ',ylab=' ')
    
    rect(xl,yl,SO+xl,SO+yl,col='wheat',border='wheat',lty=2,lwd=2)
    polygon(c(xl,SO+xl,xl),c(yl,yl,SO+yl),col='blue',border='darkblue')
    rect(0,yl/10,r,mk+yl/10,col='blue',border='wheat', lwd=2)
    
    text(xl+SO/4,yl+SO/3,bquote(Sigma == .(en)), col='wheat', cex=1.4 )
    text(r, yl/20*(mk + 1),
         paste('Z (',mk,' x ',r,' = ',mk*r,')',sep=''),col='blue',
         cex=1.,pos=4)
    .plotLabel('Dimensions','topright')
    
    if(!SAVEPLOTS){
      readline('reduction from sigma to Z -- return to continue ')
    } else {
      dev.off()
    }
  } 
  
  ########### grid/correlation analysis
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'clusterGridR.pdf') )
  
  par(mfrow=c(1,1),bty='n',cex=1, oma=oma, mar=mag, tcl= tcl, mgp=mgp)
  
  colnames(corMu) <- rownames(corMu) <- colnames(y)
  
  psize <- .62
  if(SMALLPLOTS)psize <- psize/2
  
  par(plt=c(.03,.15,.1,.9), bty='n', new=F)
  opt <- list( main=' ',cex=.2, ncluster=ncluster, 
               colCode=specColor[notOmit], textSize=.5, 
               LABELS = F, DIST=F )
  tmp <- .clusterPlot( corMu[notOmit,notOmit] , opt)
  colCode   <- tmp$colCode
  corder    <- rev(tmp$corder)
  # specOrder <- snames[notOmit[corder]]
  rOrder <- snames[notOmit[corder]]
  
  clusterIndex <- cbind( clusterIndex, tmp$clusterIndex )
  clusterOrder <- cbind( clusterOrder, tmp$corder )
  
  ncc <- ncol(clusterIndex)
  colnames(clusterIndex)[ncc] <- colnames(clusterOrder)[ncc] <- 'R'
  
  if(LABELS){
    
    par(plt=c(.15,.33,.1,.9), bty='n', new=T)
    plot(c(0,0),c(0,0),col='white',xlim=range(c(0,1)),ylim=c(0,SO),
         xaxt='n',yaxt='n',xlab='',ylab='')
    xl <- rep(.5,SO)
    
    yl <- c(1:SO) + par('usr')[3] - .75
    cex <- .fitText2Fig(rOrder,fraction=1.2)
    text( xl,yl,rev(rOrder),pos=3,cex=cex, col=rev(colCode[corder]))
  }
  
  # knames <- snames[notOmit]
  
  tmp <- .invMatZero(sgibbs,nsim=nrow(sgibbs),snames=snames,
                     knames=rOrder,index=NULL, COMPRESS=T, 
                     REDUCT=REDUCT,
                     sigErrGibbs = output$chains$sigErrGibbs, 
                     kgibbs = output$chains$kgibbs,
                     otherpar = otherpar, alpha=ematAlpha)
  marIn <- tmp$inMarMat
  conIn <- tmp$inConMat
  wm    <- which(marIn[,1] %in% omit |  marIn[,2] %in% omit)
  if(length(wm) > 0)marIn <- marIn[-wm,]
  wm    <- which(conIn[,1] %in% omit |  conIn[,2] %in% omit)
  if(length(wm) > 0)conIn <- conIn[-wm,]
  
  sigCor <- c(nrow(marIn),nrow(conIn))/SM/(SM - 1)
  sigCor <- round(100*sigCor,0)
  names(sigCor) <- c('n_marIn','n_conIn')
  
  mor <- notOmit[corder]
  
  crr <- corMu[mor,mor]
  marIn[,1] <- match(marIn[,1],mor)
  marIn[,2] <- match(marIn[,2],mor)
  conIn[,1] <- match(conIn[,1],mor)
  conIn[,2] <- match(conIn[,2],mor)
  
  makeCR <- list('white' = conIn,'grey' = marIn)
  if(!is.null(specColor))textCol = colCode[mor]
  
  par(plt=c(.33, .33 + psize,.1,.9), bty='n', new=T)
  
  slim <- quantile(crr[lower.tri(crr)],c(.05,.95))
  
  SPECLABS <- F
  if(S < 30)SPECLABS <- T
  
  .corPlot(crr, slim=slim, makeColor=makeCR,plotScale=.99,
           PDIAG=T,CORLINES=CORLINES, textCol = colCode[corder], 
           SPECLABS = SPECLABS, squarePlot = F,
           textSize=1, widex = width, widey = height, new=T, add=F)
  
  ll <- paste(c('Cond Ind (white) = ', 'Cond & Marg Ind (grey) = '),
              sigCor,c('%','%'),sep='')
  legend('topright',ll,bty='n',cex=.8)
  
  .plotLabel(expression( paste(hat(bold(R)),'structure'  )),above=T, cex=.9)
  
  if(!SAVEPLOTS){
    readline('posterior correlation for model -- return to continue ')
  } else {
    dev.off()
  }
  
  ########################### cluster Fmat with beta
  
  fBetaMu <- output$parameters$betaStandXWmu
  
  if(Q > 4){
    
    graphics.off()
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'gridF_B.pdf') ) # start plot
    
    main1 <- expression( paste('Sensitivity ',hat(F)))
    main2 <- expression( paste('Responses ',hat(B)))
    
    ws <- which( rowSums(fMat) == 0)
    if(length(ws) > 0){
      not0 <- c(1:nrow(fMat))[-ws]
      fMat <- fMat[not0,not0] 
      fBetaMu <- fBetaMu[not0,]
    }
    
    
    mat1 <- fMat
    mat2 <- fBetaMu
    expand <- ncol(mat1)/ncol(mat2)
    expand <- max(c(1.5,expand))
    
    opt <- list(mainLeft=main1, main1=main1, main2 = main2,
                leftClus=T, topClus2=T, rightLab=F, topLab1=T, 
                topLab2 = T, leftLab=T, ncluster = ncluster,
                colCode2 = specColor[notOther], lower1 = T, diag1 = T,
                lower2 = F)
    .clusterWithGrid(mat1, mat2, expand=expand, opt)
    
    if(!SAVEPLOTS){
      readline('F & beta structure -- return to continue ')
    } else {
      dev.off()
    } 
  }
  #################################### cluster Emat
  
  graphics.off()
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'clusterGridE.pdf') ) # start plot
  
  mat1 <- ematrix[notOther,notOther]
  main1 <- expression(paste('Species ',hat(E)))
  opt <- list(mainLeft=main1, leftClus=T, leftLab=T, 
              colCode1 = specColor[notOther], rowCode = specColor[notOther],
              topLab1=T,ncluster = ncluster,
              lower1 = T, diag1 = F,horiz1=clusterIndex[,'E'])
  
  .clusterWithGrid(mat1, mat2=NULL, expand=1, opt)
  
  if(!SAVEPLOTS){
    readline('E: model-based response to X -- return to continue ')
  } else {
    dev.off()
  }
  
  ################# resid and Egrid
  
  graphics.off()
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'gridR_E.pdf') ) # start plot
  
  dcor <- .cov2Cor(covy)
  dcor[is.na(dcor)] <- 0
  
  mat1 <- dcor     
  mat2 <- ematrix[notOther,notOther]
  
  main1 <- expression(paste('Ordered by error ',hat(R)))
  main2 <- expression(paste('Response ',hat(E)))
  opt <- list(mainLeft='Species', main1=main1, main2 = main2,
              leftClus=T, leftLab=T, rowCode = specColor[notOther],
              topLab1 = T, topLab2 = T,rightLab=F,ncluster = ncluster,
              lower1 = T, diag1 = F,lower2 = T, diag2 = T)
  .clusterWithGrid(mat1, mat2, expand=1, opt)
  
  if(!SAVEPLOTS){
    readline('comparison R vs E -- return to continue ')
  } else {
    dev.off()
  }
  
  ################# data vs E grid
  
  graphics.off()
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'gridY_E.pdf') ) # start plot
  
  ytmp <- jitter(y[,mor],1e-10)
  cory <- cor(ytmp)
  
  mat1 <- cory
  mat2 <- ematrix[notOther,notOther]
  main1 <- 'Ordered by data, cor(Y)'
  main2 <- expression(paste('Response ',hat(E)))
  topLab1 <- topLab2 <- F
  if(S < 30)topLab1 <- topLab2 <- T
  
  
  opt <- list(mainLeft='Species', main1=main1, main2 = main2,
              leftClus=T, leftLab=T, lower1 = T, diag1 = F,
              topLab1 = topLab1, topLab2 = topLab2,ncluster = ncluster,
              lower2 = T, diag2 = T, sameOrder = T)
  .clusterWithGrid(mat1, mat2=mat2, expand=1, opt )
  
  if(!SAVEPLOTS){
    readline('raw data vs E -- return to continue ')
  } else {
    dev.off()
  }
  
  #################### beta grid
  if(BETAGRID & nrow(output$parameters$betaStandXWmu) > 2){
    
    graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'clusterGridB.pdf') ) # start plot
    
    mat1 <- output$parameters$ematrix[notOther,notOther]
    #   mat2 <- t(betaStandXWmu[,notOther])
    mat2 <- t(output$parameters$betaStandXWmu)
    main1 <- expression(paste('Species ',hat(E)))
    main2 <- expression(paste(hat(B),' by predictor'))
    topLab1 <- F
    if(S < 30)topLab1 <- T
    
    ee <- ncol(mat1)/ncol(mat2)
    ee <- max(c(ee,.8))
    ee <- min(c(ee, 1.2))
    opt <- list(mainLeft=main1, main1=main1, main2 = main2,
                topClus1=T, topClus2=T, topLab1 = topLab1, topLab2=T,
                leftLab=T,lower1 = T, diag1 = F, ncluster = ncluster,
                colCode1 = specColor[notOther],
                vert1=clusterIndex[,'E'], horiz2=clusterIndex[,'E'])
    .clusterWithGrid(mat1, mat2, expand=ee, opt)
    
    if(!SAVEPLOTS){
      readline('beta ordered by response to X -- return to continue ')
    } else {
      dev.off()
    }
    
    ################## random groups
    
    if(RANDOM){
      
      graphics.off()
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'randGroups.pdf') ) # start plot
      G <- ncol(randByGroup)
      
      mat1 <- randGroupVarMu[notOther,notOther]
      mat1 <- .cov2Cor(mat1)
      diag(mat1) <- 0
      mat2 <- randByGroup[notOther,] # + matrix(betaMu[1,notOther],length(notOther),G, byrow=T)
      main1 <- expression(paste('Species '))
      main2 <- expression('Group')
      topLab1 <- F
      if(S < 30)topLab1 <- T
      
      ee <- ncol(mat1)/ncol(mat2)
      ee <- max(c(ee,1))
      ee <- min(c(ee, 1.2))
      opt <- list(mainLeft=main1, main1=main1, main2 = main2,leftClus=T, 
                  topClus1=F, topClus2=T, topLab1 = topLab1, topLab2=T,
                  leftLab=T, lower1 = T, diag1 = F, 
                  colCode1 = specColor[notOther])
      .clusterWithGrid(mat1, mat2, expand=ee, opt)
      
      if(!SAVEPLOTS){
        readline('random groups correlation, coeffs -- return to continue ')
      } else {
        dev.off()
      }
    }
    
    
    ###################### Time grid
    if(TIME){
      
      graphics.off()
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'clusterTime.pdf') ) 
      
      mat1 <- alphaMu[notOther,notOther]
      lam  <- lambdaMuUn[,notOther]
      lam[1,] <- lam[1,] - 1
      mat2 <- t(lam)
      colnames(mat2)[1] <- 'lambda - 1'
      main1 <- expression(paste(hat(alpha),' from'))
      side1 <- expression(paste(hat(alpha),' to'))
      main2 <- expression(hat(lambda))
      mat1[is.na(mat1)] <- 0
      mat2[is.na(mat2)] <- 0
      
      topLab1 <- F
      if(S < 20)topLab1 <- T
      
      ee <- ncol(mat1)/(ncol(mat1) + ncol(mat2) )
      #  ee <- max(ee,.3)
      slim1 <- range(mat1)
      if(slim1[2] == 0)slim1[2] <- .0001
      opt <- list(mainLeft=side1, main1=main1, main2 = main2,
                  ncluster = ncluster,
                  topClus1=F, topClus2=F, topLab1 = topLab1, 
                  topLab2=T, rowOrder = c(1:S)[notOther], colOrder1 = c(1:S)[notOther], 
                  colOrder2 = 1:ncol(mat2), slim1 = slim1,
                  colCode1 = boxCol[notOther], lower1 = F, diag1 = F)
      .clusterWithGrid(mat1, mat2, expand=ee, opt)
      
      if(!SAVEPLOTS){
        readline('beta ordered by response to X -- return to continue ')
      } else {
        dev.off()
      }
      
      graphics.off()
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,'clusterGridLambda.pdf') ) # start plot
      
      mat1 <- ematrix[notOther,notOther]
      main1 <- expression(paste('Species ',hat(E)))
      main2 <- expression(paste(hat(Lambda),' by predictor'))
      topLab1 <- F
      if(S < 40)topLab1 <- T
      
      ee <- ncol(mat1)/(ncol(mat1) + ncol(mat2) )
      #  ee <- max(ee,.05)
      opt <- list(mainLeft=main1, main1=main1, main2 = main2,
                  colOrder2 = 1:ncol(mat2), ncluster = ncluster,
                  topClus1=T, topClus2=T, topLab1 = topLab1, topLab2=T,
                  colCode1 = boxCol[notOther], lower1 = T, diag1 = F)
      #                vert1=clusterIndex[,'E'], horiz2=clusterIndex[,'E'])
      .clusterWithGrid(mat1, mat2, expand=ee, opt)
      
      if(!SAVEPLOTS){
        readline('lambda ordered by response to X -- return to continue ')
      } else {
        dev.off()
      }
    }
    
    if(TRAITS){
      if(nrow(betaTraitMu) > 3){
        
        bb <- betaTraitMu[-1,]
        ord <- order(colSums(abs(bb)),decreasing=T)
        bb  <- bb[,ord]
        bl  <- bb[,ord]
        bh  <- bb[,ord]
        ror <- order(rowSums(abs(bb)),decreasing=T)
        bb  <- bb[ror,]
        bl  <- bl[ror,]
        bh  <- bh[ror,]
        
        white <- which(bl < 0 & bh > 0,arr.ind=T)
        
        makeColor <- list('white' = white )
        
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,'gridTraitB.pdf') ) 
        
        plotScale <- max(c(10,c(S,Q)/10))
        
        par(mfrow=c(1,1), bty='n', oma=c(1,1,1,1), 
            mar=c(5,4,4,2), tcl= tcl, mgp=mgp)
        
        ht <- nrow(bb)/ncol(bb)*width
        
        opt <- list(mainLeft='', main1='', main2 = '',
                    topClus1=T, topClus2=T, topLab1 = T, topLab2=F,
                    leftClus=T, 
                    leftLab=T, ncluster = ncluster,
                    colCode1 = traitColor)
        .clusterWithGrid(mat1=betaTraitMu[-1,], mat2=NULL, expand=1, opt)
        
        if(!SAVEPLOTS){
          readline('trait beta -- return to continue ')
        } else {
          dev.off()
        }
      }
    }
  }
  
  all <- list(fit = fit, ematrix = ematrix,
              eComs = eComs, ncluster = ncluster,
              clusterIndex = clusterIndex, clusterOrder = clusterOrder,
              eVecs = eVecs, eValues = eValues) 
  all <- all[ order(names(all)) ]
  invisible(all)
}

.gjamPrediction <- function(output, newdata, y2plot, PLOT, ylim, FULL){
  
  xnew <- ydataCond <- interBeta <- groupRandEff <- NULL
  tiny  <- 1e-10
  wHold <- phiHold <- ploHold <- sampleWhold <- NULL
  COND <- RANDOM <- F
  
  ng     <- output$modelList$ng
  burnin <- output$modelList$burnin
  
  nsim <- 500
  if('nsim' %in% names(newdata))nsim <- newdata$nsim
  
  if( is.null(newdata) ){
    
    if(PLOT){
      
      y1 <- output$inputs$y
      y2 <- output$prediction$ypredMu
      if(!is.null(y2plot)){
        y1 <- y1[,y2plot]
        y2 <- y2[,y2plot]
      }
      
      tmp <- .bins4data(y1)
      breaks <- tmp$breaks
      bins   <- tmp$bins
      nbin   <- tmp$nbin
      
      if(length(bins) > 0){
        breaks <- bins
        nPerBin <- NULL
      }
      
      opt <- list(nPerBin = NULL, breaks=breaks, ylimit = ylim,
                  fill='lightblue', box.col='darkblue', POINTS=F)
      
      .plotObsPred(y1, y2, opt = opt)
      abline(0,1,lwd=4,col='white')
      abline(0,1,lwd=2,col='grey',lty=2)
    }
    return(  list( ypredMu = output$modelSummary$ypredMu, 
                   ypredSe = output$modelSummary$ypredSd ) )
  }
  
  S <- SO <- S1 <- ncol(output$inputs$y)
  Q <- ncol(output$inputs$x)
  n <- nrow(output$inputs$x)
  y <- yp <- output$inputs$y
  x <- output$inputs$x
  
  xnames <- colnames(x)
  ynames <- colnames(y)
  
  cindex <- NULL
  
  notOther <- output$inputs$notOther
  other    <- output$inputs$other
  SO       <- length(notOther)
  
  otherpar <- output$modelList$reductList$otherpar
  censor   <- output$modelList$censor
  REDUCT   <- output$modelList$REDUCT
  
  notStandard <- output$modelList$notStandard
  
  NEWX <- F
  if('xdata' %in% names(newdata))NEWX <- T
  if('ydataCond' %in% names(newdata))COND <- T
  
  effort <- output$modelList$effort
  effMat <- effort$values
  inSamp <- 1:n
  
  REDUCT <- output$modelList$REDUCT
  sigmaerror <- NULL
  if(REDUCT){
    otherpar <- output$modelList$reductList$otherpar
    N  <- otherpar$N
    r  <- otherpar$r
    rndEff <- y*0
    sigmaerror <- otherpar$sigmaerror
  }
  cuts <- output$parameters$cutMu
  cuts <- cbind(-Inf,0,cuts,Inf)
  
  nfact      <- output$inputs$factorBeta$nfact
  isFactor   <- output$inputs$factorBeta$isFactor
  factorList <- output$inputs$factorBeta$factorList
  contrasts  <- output$inputs$factorBeta$contrast
  formula    <- output$modelList$formula
  xscale     <- output$inputs$standX
  if(is.matrix(xscale))  xscale <- t(xscale)
  facNames   <- names(factorList)
  
  typeNames <- output$modelList$typeNames
  tmp       <- .gjamGetTypes(typeNames)
  typeFull  <- tmp$typeFull
  typeCols  <- tmp$typeCols
  allTypes  <- unique(typeCols)
  typeCode  <- tmp$TYPES[typeCols]
  FCgroups  <- attr(typeNames,'FCgroups')
  CCgroups  <- attr(typeNames,'CCgroups')
  CATgroups <- attr(typeNames,'CATgroups')
  
  condCols <- numeric(0)
  
  standRows <- output$inputs$standRows
  standMat  <- output$inputs$standMat
  standX    <- output$inputs$standX
  xmu       <- standX[,1]
  xsd       <- standX[,2]
  intMat    <- interBeta$intMat
  
  notCorCols <- 1:S
  
  if( NEWX ){   ################ out-of-sample
    
    xnew <- newdata$xdata
    nx   <- n <- nrow(xnew)
    colnames(xnew) <- .cleanNames(colnames(xnew))
    
    wna <- which(is.na(xnew),arr.ind=T)
    if(length(wna) > 0)
      stop('cannot have NA in prediction grid newdata$xdata')
    
    
    effMat <- matrix(1, nx, S)
    holdoutN <- nx
    holdoutIndex <- 1:nx
    
    if( 'effort' %in% names(newdata) ){
      ev     <- newdata$effort$values
      effMat <-  matrix(1, nx, S)
      effMat[,newdata$effort$columns] <- ev
    }
    effort <- list(columns = c(1:S), values = effMat)
    
    ydataCond <- NULL
    
    if(nfact > 0){
      for(j in 1:nfact){
        
        nf <- names(factorList)[j]
        wf <- which(names(xnew) == nf)
        wo <- which(names(output$xnew) == nf)
        wc <- which(names(contrasts) == names(factorList)[j])
        cc <- contrasts[[wc]]
        
        xnew[[wf]] <- factor( xnew[[wf]], levels = rownames(cc) )
        
        attr(xnew[[wf]],'contrasts') <- cc
      }
    }
    
    y <- matrix(0,nx,S)
    colnames(y) <- ynames
    yp <- y
    
    wss <- names(standRows)[names(standRows) %in% names(xnew)]
    
    xnew[,wss] <- t( (t(xnew[,wss]) - standX[wss,'xmean'])/
                       standX[wss,'xsd'])
    
    tmp <- .gjamXY(formula, xnew, yp, typeNames, 
                   notStandard=names(xnew), checkX = F, xscale = xscale)
    x  <- tmp$x     
    
    beta <- output$parameters$betaMu
    
    w  <- x%*%beta
    yp <- w*effMat
    
    wca <- which(typeNames == 'CA')
    if(length(wca) > 0){
      yp[,wca][yp[,wca] < 0] <- 0
    }
    
    wda <- which(typeNames == 'DA')
    if(length(wda) > 0){
      yp[,wda] <- round(yp[,wda]*effMat[,wda],0)
      yp[,wda][yp[,wda] < 0] <- 0
    }
    
    ordCols <- which(typeNames == 'OC')
    if(length(ordCols) > 0){
      
      tmp   <- .gjamGetCuts(yp + 1,ordCols)
      cutLo <- tmp$cutLo
      cutHi <- tmp$cutHi
      
      for(k in ordCols){
        yp[,k] <- findInterval(yp[,k],cuts[k,]) - 1
      }
    }
    
    if(length(FCgroups) > 0){
      ntt <- max(FCgroups)
      for(i in 1:ntt){   
        wk      <- which( FCgroups == i )
        wo      <- which(wk %in% notOther)
        yp[,wk] <- .gjamCompW2Y(yp[,wk,drop=F], notOther=wo)$ww
      }
    }
    
    if(length(CCgroups) > 0){
      
      print( 'for CC data total effort (count) is taken as 1000' )
      ysum <- rep(1000,n)                   # CC use sum of 100
      ntt  <- max(CCgroups)
      if(ntt > 0){
        for(i in 1:ntt){  ## normalize y 
          wk      <- which( CCgroups == i )
          wo      <- which(wk %in% notOther)
          yp[,wk] <- .gjamCompW2Y(yp[,wk,drop=F], notOther=wo)$ww
          yp[,wk][yp[,wk] < 0] <- 0
          yp[,wk] <- round( sweep(yp[,wk],1,ysum,'*'), 0) 
        }
      }
    }
    
    
    tmp <- .gjamSetup(typeNames, x, yp, breakList=NULL, holdoutN=NULL, 
                      holdoutIndex=NULL, censor=NULL, effort=effort) 
    w <- tmp$w; z <- tmp$z; yp <- tmp$y; other <- tmp$other
    plo <- tmp$plo; phi <- tmp$phi
    ordCols  <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
    minOrd   <- tmp$minOrd;   maxOrd <- tmp$maxOrd;  censorCA <- tmp$censorCA
    censorDA <- tmp$censorDA;   censorCON <- tmp$censorCON;   ncut <- ncol(cuts)
    corCols <- tmp$corCols
    if(length(corCols) > 0)notCorCols <- notCorCols[-corCols]
    catCols  <- which(attr(typeNames,'CATgroups') > 0)
    sampleW  <- tmp$sampleW*0 + 1
    
    byCol <- byRow <- F
    if(attr(sampleW,'type') == 'cols')byCol <- T
    if(attr(sampleW,'type') == 'rows')byRow <- T
    indexW <- attr(sampleW,'index')
    inSamp <- 1:n
    
    byCol <- byRow <- F
    if(attr(sampleW,'type') == 'cols')byCol <- T
    if(attr(sampleW,'type') == 'rows')byRow <- T
    indexW <- attr(sampleW,'index')
    
    cdex <- c(1:S)
  }
  
  if(COND){
    
    ydataCond <- newdata$ydataCond
    colnames(ydataCond) <- .cleanNames(colnames(ydataCond))
    condNames <- colnames(ydataCond)
    
    if('other' %in% condNames){
      condNames <- condNames[condNames != 'other']
      ydataCond <- ydataCond[,condNames]
    }
    
    n  <- nrow(x)
    yp <- y
    
    condCols <- match(condNames, colnames(yp))
    
    yp[,condCols] <- as.matrix( ydataCond )
    
    tmp <- .gjamSetup(typeNames, x, yp, breakList=NULL, holdoutN=NULL, 
                      holdoutIndex=NULL,censor=NULL, effort=effort) 
    w <- tmp$w; z <- tmp$z; yp <- tmp$y; other <- tmp$other
    plo      <- tmp$plo; phi <- tmp$phi
    ordCols  <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
    minOrd   <- tmp$minOrd;   maxOrd <- tmp$maxOrd;  censorCA <- tmp$censorCA
    cuts     <- tmp$cuts
    censorDA <- tmp$censorDA;   censorCON <- tmp$censorCON; ncut <- ncol(cuts)
    corCols <- tmp$corCols
    if(length(corCols) > 0)notCorCols <- notCorCols[-corCols]
    effort   <- tmp$effort
    catCols  <- which(attr(typeNames,'CATgroups') > 0)
    sampleW  <- tmp$sampleW
    sampleW[,-condCols] <- 1
    
    standRows <- output$inputs$standRows
    standMat  <- output$inputs$standMat
    standMu <- output$inputs$standMu
    
    byCol <- byRow <- F
    if(attr(sampleW,'type') == 'cols')byCol <- T
    if(attr(sampleW,'type') == 'rows')byRow <- T
    indexW <- attr(sampleW,'index')
    
    cdex <- c(1:S)[-condCols]
    
    CCsums <- numeric(0)
    if(!is.null(CCgroups)){
      ncc    <- max(CCgroups)
      for(j in 1:ncc){
        wjk    <- which(CCgroups == j)
        CCsums <- append(CCsums,list( rowSums(y[,wjk]) ) )
      }
    }
  } ##############################
  
  if(length(other) > 0)cdex <- cdex[!cdex %in% other]
  S1   <- length(cdex)
  
  yg <- yp
  
  if(length(yp) < 10000 | FULL) FULL <- T
  
  if(FULL){
    ygibbs <- wgibbs <- matrix(0,nsim,length(yp))
  }
  
  #partition out-of-sample based max ever obs for species
  pmax <- apply(output$inputs$y/output$modelList$effort$values,2,max) 
  ptmp <- 10*matrix(pmax,n,S,byrow=T)                 
  
  ptmp[,ordCols]  <- length(ordCols) + 10
  ptmp[,compCols] <- 10
  ptmp[,catCols]  <- 10
  
  # note: all are holdouts for newdata, no holdouts for COND
  
  if(COND){
    holdoutN <- 0
    holdoutIndex <- NULL
    ploHold <- phiHold <- NULL
    plo[,-condCols] <- -ptmp[,-condCols]
    phi[,-condCols] <- ptmp[,-condCols]
  }else{
    holdoutN <- n
    holdoutIndex <- c(1:n)
    ploHold <- plo
    phiHold <- phi
    plo <- -ptmp
    phi <- ptmp
  }
  
  .updateW <- .wWrapper(REDUCT, RANDOM, S, effMat, corCols, notCorCols, typeNames, 
                        typeFull, typeCols, 
                        allTypes, holdoutN, holdoutIndex, censor, 
                        censorCA, censorDA, censorCON, notOther, sampleW, 
                        byRow, byCol,
                        indexW, ploHold, phiHold, sampleWhold, inSamp)
  
  ypred  <- matrix(0,n,S)
  colnames(ypred) <- ynames
  ypred2 <- wcred <- wcred2 <- ypred
  gvals  <- sample(burnin:ng,nsim,replace=T)
  
  pbar <- txtProgressBar(min=1,max=nsim,style=1)
  ig   <- 0
  
  corColC <- cdex[cdex %in% corCols]
  corColW <- which(cdex %in% corCols)
  
  ddex <- which(notOther %in% cdex)
  
  cutg <- cuts
  ncut <- ncol(cutg)
  
  ccols <- which(typeNames != 'CON')
  
  kg     <- 1
  
  rndEff <- 0
  
  prPresent <- w*0
  
  ############ E matrix
  emat <- matrix(0,S,S)
  colnames(emat) <- rownames(emat) <- ynames
  lo <- hi <- lm <- hm <- ess <- emat
  
  eCont <- output$inputs$factorBeta$eCont
  dCont <- output$inputs$factorBeta$dCont
  lCont <- output$inputs$factorBeta$lCont
  
  covE <- cov( x%*%dCont )  # note that x is standardized
  
  frow  <- NULL
  if(nfact > 0){
    frow <- rep(0,Q)
    for(j in 1:nfact){
      frow[ match(factorList[[j]], xnames) ] <- j
    }
  }
  
  q1 <- nrow(eCont)
  fnames   <- rownames(eCont)
  facList2 <- factorList
  if(nfact > 0){
    for(j in 1:nfact){
      wj <- which(names(xnew) == names(factorList)[j])
      facList2[[j]] <- levels(xnew[[wj]])
    }
  }
  notPA <- which(!typeNames == 'PA' & !typeNames == 'CON')
  
  for(g in gvals){
    
    bg  <- matrix( output$chains$bgibbs[g,], Q, S)
    muw <- x%*%bg
    
    if(REDUCT){
      Z  <- matrix(output$chains$sgibbs[g,],N,r)
      sigmaerror <- output$chains$sigErrGibbs[g]
      K    <- output$chains$kgibbs[g,]
      sg   <- .expandSigma(sigmaerror, S, Z = Z, K, REDUCT = T)
    } else {
      sg <- .expandSigma(output$chains$sgibbs[g,], S = S, REDUCT = F)
    }
    
    alpha <- .sqrtRootMatrix(bg,sg,DIVIDE=T)
    
    bgg <- bg[,notOther]
    agg <- .sqrtRootMatrix(bgg,sg[notOther,notOther],DIVIDE=T)  
    
    if(nfact > 0){
      agg <- lCont%*%agg    #standardized for x and cor scale for y
      for(k in 1:nfact){
        fk  <- factorList[[k]]
        mua <- colMeans(agg[drop=F,fk,])
        nl  <- length(fk)
        agg[fk,] <- agg[fk,] - matrix(mua,nl,SO,byrow=T)
      }
    } else {
      agg <- agg[drop=F,-1,]
    }
    
    egg         <- lCont%*%bgg          #standardized for x, not cor for y
    
    if( 'OC' %in% typeCode ){
      cutg[,3:(ncut-1)] <- matrix( output$chains$cgibbs[g,], length(ordCols))
      tmp   <- .gjamGetCuts(yg + 1,ordCols)
      cutLo <- tmp$cutLo
      cutHi <- tmp$cutHi
      
      plo[,ordCols] <- cutg[cutLo]
      phi[,ordCols] <- cutg[cutHi]
    }
    
    tmp <- .updateW( rows=1:n, x, w, yg, bg, sg, alpha, cutg, plo, phi, 
                     rndEff=rndEff, groupRandEff, sigmaerror, wHold )
    w   <- tmp$w
    
    if(!COND){
      
      yg  <- tmp$yp   
      
    }else{
      
      tmp <- .conditionalMVN(w, muw, sg, cdex = ddex, S)  
      muc    <- tmp$mu
      sgp    <- tmp$vr
      if(S1 == 1){
        w[,ddex] <- matrix(rnorm(n,muc,sqrt(sgp[1])))
      } else {
        w[,ddex] <- .rMVN(n,muc,sgp)
      }
      muw[,ddex] <- muc
      
      if( length(corColC) > 0 ){    #expanded w on this scale
        sgs  <- .cov2Cor(sg)
        mus  <- x%*%alpha
        muw[,corColC] <- mus[,corColC]
        
        tmp <- .conditionalMVN(w, mus, sgs, cdex = cdex, S)
        mus    <- tmp$mu
        sgs    <- tmp$vr
        muw[,cdex] <- mus
        
        if(S1 == 1){
          w[,ddex] <- matrix(rnorm(n,mus,sqrt(sgs[1])))
        } else {
          w[,ddex] <- .rMVN(n,mus,sgs)
        }
      } 
      yg[,-condCols] <- (w*effMat)[,-condCols]
      if(length(notPA) > 0){
        mmm <- yg[,notPA]
        mmm[mmm < 0] <- 0
        yg[,notPA]   <- mmm
      } 
      
      for(k in allTypes){    # predicting from w (not from yg)
        
        wk  <- which(typeCols == k)
        nk  <- length(wk)
        wo  <- which(wk %in% notOther)
        wu  <- which(typeCols[notOther] == k)
        wp  <- w[, wk, drop=F]
        
        groups <- NULL
        
        if( typeFull[wk[1]] == 'countComp' ){
          
          groups <- CCgroups[wk]
          nkk    <- max(groups)
          
          for(j in 1:nkk){
            
            wjk <- which(typeCols[wk] == k & CCgroups[wk] == j)
            wno <- which(wk %in% notOther)
            woo <- which(wk %in% other)
            www <- w[,wk]
            www[www < 0] <- 0
            
            www <- .gjamCompW2Y(www,notOther=wno)$ww
            
            if(COND){
              www <- sweep(www,1,CCsums[[j]],'*')
            } else {
              www <- sweep(www,1,ysum,'*')
            }
            yg[,wk] <- www
          }
          
        } else {
          
          if(typeFull[wk[1]] == 'fracComp') groups <- FCgroups[wk]
          
          glist <- list(wo = wo, type = typeFull[wk[1]], yy = yg[,wk,drop=F], 
                        wq = wp, yq = yg[,wk,drop=F], cutg = cutg, 
                        censor = censor, censorCA = censorCA, 
                        censorDA = censorDA, censorCON = censorCON, 
                        eff = effMat[,wk,drop=F], groups = groups, 
                        k = k, typeCols = typeCols, notOther = notOther,
                        wk = wk, sampW = sampleW[,wk])
          
          tmp <- .gjamWLoopTypes( glist )
          yg[,wk] <- tmp[[2]] #[,wk]
          yg[,wk] <- .censorValues(censor,yg,yg)[,wk]
        }
      }
      
      yg[,condCols] <- as.matrix( ydataCond )
    }
    ####################
    
    if(length(ccols) > 0){
      mmm <- muw[,ccols]
      mmm[mmm < 0] <- 0
      muw[,ccols] <- mmm
    }
    
    yy <- yg
    
    if('PA' %in% typeNames){
      wpa <- which(typeNames == 'PA')
      yy[,wpa] <- round(yg[,wpa])
    }
    
    if(length(notPA) > 0){
      w0 <- which(yy[,notPA] <= 0)
      w1 <- which(yy[,notPA] > 0)
      yy[,notPA][w0] <- 0
      yy[,notPA][w1] <- 1
    }
    
    prPresent <- prPresent + yy
    
    ig <- ig + 1
    setTxtProgressBar(pbar,ig)
    
    ypred  <- ypred + yg
    ypred2 <- ypred2 + yg^2
    wcred  <- wcred + muw
    wcred2 <- wcred2 + muw^2
    
    ess[notOther,notOther]  <- .cov2Cor( t(agg)%*%covE%*%agg ) 
    emat[notOther,notOther] <- emat[notOther,notOther] + ess[notOther,notOther]
    
    if(FULL){
      ygibbs[kg,] <- as.vector(yg)
      wgibbs[kg,] <- as.vector(muw)
    }
    kg <- kg + 1
  } ###################
  
  prPresent <- prPresent/nsim
  
  ematrix  <- emat/nsim
  
  xunstand <- .getUnstandX(x,standRows,xmu,xsd,intMat)$xu
  
  yMu  <- ypred/nsim
  res  <- ypred2/(nsim - 1) - yMu^2
  res[res < tiny] <- tiny
  yPe <- sqrt(res) 
  
  wMu  <- wcred/nsim
  res  <- wcred2/(nsim - 1) - wMu^2
  res[res < tiny] <- tiny
  wSe <- sqrt(res)
  
  colnames(yMu) <- colnames(yPe) <- colnames(wMu) <- 
    colnames(wSe) <- ynames
  
  sdList <- list( yMu = yMu, yPe = yPe, wMu = wMu, wSe = wSe )
  
  piList <- NULL
  if(FULL){
    wLo <- matrix( apply(wgibbs,2,quantile,.05), n, S )
    wHi <- matrix( apply(wgibbs,2,quantile,.95), n, S )
    yLo <- matrix( apply(ygibbs,2,quantile,.05), n, S )
    yHi <- matrix( apply(ygibbs,2,quantile,.95), n, S )
    
    colnames(wLo) <- colnames(wHi) <- colnames(yLo) <- 
      colnames(yHi) <-  ynames
    piList <- list( wLo = wLo, wHi = wHi, yLo = yLo, yHi = yHi )
  }
  
  if(PLOT){
    
    oma <- c(0,0,0,0)
    mar <- c(4,4,2,1)
    tcl <- -0.5
    mgp <- c(3,1,0)
    
    par(oma = oma, mar = mar, tcl = tcl, mgp = mgp, bty='n')
    
    wy <- which(colnames(y) %in% y2plot & c(1:S) %in% notOther)
    t2plot <- typeNames[wy]
    allTypes <- unique(t2plot)
    
    mfrow <- .getPlotLayout(length(allTypes) + 1)
    par(mfrow=mfrow, bty='n', mar=c(1,2,3,1) )
    
    k   <- 0
    add <- F
    
    for(j in 1:length(allTypes)){
      
      wk <- which(typeNames == allTypes[j] & c(1:S) %in% notOther)
      ws <- colnames(y)[wk]
      wm <- which(colnames(yMu) %in% colnames(y)[wk])
      wk <- match(colnames(yMu)[wm],colnames(y))
      
      y1 <- y[,wk]
      if(min(y1) == max(y1))next
      y2 <- yMu[,wm]
      
      tmp <- .gjamPlotPars(type=allTypes[j],y1,y2)
      y1 <- tmp$y1; yp <- tmp$yp; nbin <- tmp$nbin; nPerBin <- tmp$nPerBin
      vlines <- tmp$vlines; xlimit <- tmp$xlimit; ylimit <- tmp$ylimit
      breaks <- tmp$breaks; wide <- tmp$wide; LOG <- tmp$LOG; POINTS <- F
      MEDIAN <- tmp$MEDIAN
      
      log <- ''
      if(LOG)log <- 'xy'
      
      if(LOG){
        wn <- which(y1 > 0 & yp > 0)
        y1 <- y1[wn]
        yp <- yp[wn]
      }
      
      tmp <- .bins4data(y1,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
      breaks <- tmp$breaks
      bins   <- tmp$bins
      nbin   <- tmp$nbin
      
      if( !allTypes[j] %in% c('PA','CAT') ){
        ncc   <- max( c(100,max(y1)/20) )
        xy <- .gjamBaselineHist(y1,bins=bins,nclass=ncc)
        xy[2,] <- ylimit[1] + .8*xy[2,]*diff(ylimit)/max(xy[2,])
        plot(xy[1,],xy[2,],col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
             xlab='Observed',ylab='Predicted')
        polygon(xy[1,],xy[2,],border='tan',col='wheat')
        
      } else {
        y11 <- mean(y1)
        y00 <- 1 - y11
        x11 <- c(-.07,-.07,.07,.07,.93,.93,1.07,1.07,-.07)
        y11 <- c(0,y00,y00,0,0,y11,y11,0,0)
        plot(x11,y11,col='tan',type='s',lwd=2,xlim=xlimit,ylim=ylimit,
             xlab='Observed',ylab='Predicted')
        polygon(x11,y11,border='tan',col='wheat')
      }
      abline(0,1,lty=2,lwd=3,col='brown')
      abline(h = mean(y1),lty=2,lwd=3,col='tan')
      
      add <- T
      
      
      opt <- list(xlabel='Observed',ylabel='Predicted',nbin=nbin, 
                  nPerBin=nPerBin, xlimit=xlimit,ylimit=ylimit,
                  breaks=breaks, wide=wide, LOG=LOG, fill='lightblue', 
                  box.col='darkblue',POINTS=F, MEDIAN=MEDIAN, add=add)
      
      .plotObsPred(y1, y2, opt = opt)
      
      if(length(vlines) > 0)abline(v=vlines,lty=2)
      
      tt <- allTypes[j]
      if(length(ws) == 1)tt <- paste(ws,tt,sep='-')
      
      lab <- paste(letters[j],') ',tt, sep='')
      .plotLabel( lab,above=T )
    }
    
    yp  <- colMeans(yMu)
    wy  <- match(colnames(yMu),colnames(y))
    
    
    opt <- list(xlabel='Observed', xlimit=NULL, ylimit=NULL,
                breaks=breaks, wide=wide, LOG=LOG, fill='lightblue', 
                box.col='darkblue', POINTS=T, ptcol='darkblue')
    .plotObsPred( colMeans(y[,wy]),yp, opt = opt)
    
    abline(0, 1,lty=2,lwd=3,col='brown')
    abline(h = mean(y1),lty=2,lwd=3,col='tan')
    .plotLabel( paste(letters[j+1],') By Species',sep=''),above=T )
    
  }
  
  bk <- list( x = xunstand, sdList = sdList, piList = piList, prPresent = prPresent,
              ematrix = ematrix)
  if(FULL)bk <- append( bk, list(ychains = ygibbs) )
  bk
}

.updateBetaTime <- function(X, Y, sig, rows, pattern, lo=NULL, hi=NULL){
  
  SS <- ncol(Y)
  B  <- t(lo)*0
  tiny <- 1e-5
  
  QX  <- ncol(X)
  XX <- crossprod(X)
  
  omega <- sig*solveRcpp(XX)
  muB   <- t(omega%*%crossprod((1/sig)*X, Y))
  
  for(k in 1:nrow(rows)){
    krow <- rows[k,]
    krow <- krow[is.finite(krow)]
    notk <- c(1:QX)[-krow]
    if(length(notk) == 1){
      M1 <- omega[krow,notk, drop=F]/omega[notk,notk]
    }else{
      OI <- try( solveRcpp(omega[notk,notk]), T)
      if( inherits(OI,'try-error') ){
        OI <- diag(1/diag(omega[notk,notk]))
      }
      M1 <- omega[krow,notk, drop=F]%*%OI
    }
    pk  <- pattern[k,]
    pk  <- pk[is.finite(pk)]
    muk <- muB[pk, krow, drop=F] - muB[pk,notk]%*%t(M1)
    Mk  <- omega[krow,krow] - M1%*%omega[notk,krow]
    
    if(is.null(lo)){
      if(length(Mk) == 1){
        B[pk,krow] <- rnorm(length(pk),muk,sqrt(Mk))
      }else{
        B[pk,krow] <- .rMVN( length(pk), rep(0,length(krow)), Mk) + muk
      }
    } else { 
      if(length(Mk) == 1){
        B[pk,krow] <- .tnorm(length(pk),lo[krow,pk],hi[krow,pk],muk,sqrt(Mk))
      } else {
        ll <- t(lo)[pk,krow,drop=F]
        hh <- t(hi)[pk,krow,drop=F]
        test <- try( .tnormMVNmatrix( avec=muk, muvec=muk, smat=Mk,
                                      lo=ll, hi=hh), T)
        if( inherits(test,'try-error') ){
          mm <- diag(Mk)
          mm[mm < tiny] <- tiny
          test <- .tnorm(length(ll),ll,hh,muk,sqrt(mm))
        }
        B[pk,krow] <- test
      }
    }
  }
  t(B)
}

.updateTheta <- function(w,tg,cutLo,cutHi,ordCols,holdoutN,
                         holdoutIndex,minOrd,maxOrd){
  
  word <- w[,ordCols,drop=F]
  ncut <- ncol(tg)
  nc   <- ncut - 1
  n    <- nrow(w)
  nk   <- length(ordCols)
  
  c1 <- cutLo[,1]
  c2 <- cutLo[,2]
  c3 <- cutHi[,1]
  c4 <- cutHi[,2]
  
  if(holdoutN > 0){
    word <- word[-holdoutIndex,]
    ss   <- seq(0,(nk-1)*n,by=n)
    wh <- as.vector( outer(holdoutIndex,ss,'+') )
    c1 <- c1[-wh]
    c2 <- c2[-wh]
    c3 <- c3[-wh]
    c4 <- c4[-wh]
  }
  
  cmin <- .byGJAM(as.vector(word),c1,c2,fun='min')
  cmax <- .byGJAM(as.vector(word),c1,c2,fun='max')
  
  cmin[!is.finite(cmin[,1]),1] <- -10
  cmin[,2] <- 0
  cmax[,1] <- 0
  cmax[cmax == -Inf] <- Inf
  
  tmp <- .interpRows(cmax,startIndex=minOrd+1,endIndex=maxOrd-1,
                     INCREASING=T,minVal=0,maxVal=Inf,
                     defaultValue=NULL,tinySlope=.001)
  
  cmax[!is.finite(cmax)] <- tmp[!is.finite(cmax)]
  
  ww <- which(!is.finite(cmin) & is.finite(cmax),arr.ind=T)
  if(length(ww) > 0){
    w0 <- ww
    w0[,2] <- w0[,2] - 1
    cmin[ww] <- runif(nrow(ww),cmax[w0],cmax[ww])
  }
  
  clo <- cmax[drop=F,,-nc]
  chi <- cmin[drop=F,,-1]
  clo[,1] <- -1
  
  ww <- which(is.finite(clo))
  cl <- clo[ww]
  ch <- chi[ww]
  wc <- which(cl > ch,arr.ind=T)
  cl[cl > ch] <- ch[cl > ch]
  
  chi[ww] <- .tnorm(length(ww),cl,ch,cl,3)
  chi[,1] <- 0
  cmax <- cbind(-Inf,chi,Inf)
  
  
  cmax[,ncut] <- Inf
  if( ncol(cmax) > max(maxOrd) )cmax[ cbind(1:nk,maxOrd+1) ] <- Inf
  
  wmin <- which(minOrd > 1)
  if(length(wmin) > 0){
    for(j in wmin)cmax[j,2:c(minOrd[j]+1)] <- 0:(minOrd[j] - 1)
  }
  cmax
}

.censorValues <- function(censor,y,yp){
  
  mm  <- length(censor)
  if(mm == 0)return(yp)
  
  if(mm > 0){
    for(m in 1:mm){
      wc  <- censor[[m]]$columns
      nc  <- ncol( censor[[m]]$partition )
      ym  <- yp[,wc,drop=F]
      cp  <- censor[[m]]$partition
      for(k in 1:nc){
        wlo <- which( ym >  cp[2,k] & ym <  cp[3,k])
        ym[wlo] <- cp[1,k]
      }
      yp[,wc] <- ym
    }
  }
  yp
}

.gjamWLoopTypes <- function( glist ){
  
  wo <- type <- yy <- wq <- yq <- cutg <- censor <- 
    censorCA <- censorDA <- censorCON <- eff <- groups <- k  <- 
    typeCols <- notOther <- wk <- sampW <- NULL
  
  for(k in 1:length(glist))assign( names(glist)[k], glist[[k]] )
  
  #returns [[1]] in-sample w for x prediction, and 
  #        [[2]] out-of-sample y prediction
  
  if( type == 'continuous' ){
    yy[sampW == 1] <- wq[sampW == 1]
    return( list(yy,yq) )   # w = y
  }
  
  nk  <- ncol(wq)
  wkk <- c(1:nk)
  n  <- nrow(wq)
  
  if( type == 'ordinal' ){
    for(s in 1:nk)yq[,s] <- findInterval(yq[,s],cutg[s,]) - 1
    return( list(wq,yq) )
  }
  
  if( type == 'presenceAbsence' ){
    yq <- pnorm(yq)          # probit
    return( list(wq,yq) )
  }
  
  if( type == 'contAbun' ){
    yq[yq < 0]    <- 0
    return( list(wq,yq) )
  }
  
  if( type == 'discAbun' ){
    yq[yq < 0] <- 0
    if(length(censorDA) > 0) wq[-censorDA] <- yy[-censorDA]
    
    yq <- yq*eff
    
    return( list(wq,yq) )
  }
  
  if( type == 'categorical' ){ ## only prediction
    
    ntt <- max( groups )
    
    for(i in 1:ntt){  
      
      if(ntt == 1){
        wki <- wkk
      } else {
        wki <- which( groups == i )
      }
      nki  <- length(wki)
      wko  <- wki
      wmax <- apply( yq[,wko],1, which.max) 
      
      yq[,wki] <- 0
      yq[,wki][ cbind(1:n,wmax) ] <- 1
    }
    return( list(wq,yq) )
  }
  
  if( type == 'countComp' ){  ##  w and y 
    
    ntt <- max( groups )
    ww  <- wq                 
    ww[ww < 0] <- 0     
    yq[yq < 0] <- 0
    
    for(i in 1:ntt){  ## normalize w and y 
      
      if(ntt == 1){
        wki <- wkk
      } else {
        wki <- which( groups == i )
      }
      
      io <- which(wki %in% wo)
      wc <- .gjamCompW2Y(ww[,wki,drop=F], notOther=io)$ww
      wq[,wki][wq[,wki] > 0] <- wc[wq[,wki] > 0]
      
      yq[,wki] <- .gjamCompW2Y(yq[,wki,drop=F],notOther=io)$ww
      ysum     <- rowSums(yy[,wki,drop=F])
      yq[,wki] <- round( sweep(yq[,wki,drop=F],1,ysum,'*'), 0) 
    }
    return( list(wq,yq) )
  }
  
  ## fracComp: w and y 
  
  ntt <- max( groups )
  
  wy     <- which(yy > 0)  
  wq[wy] <- yy[wy]     
  yq[yq < 0] <- 0
  
  for(i in 1:ntt){  ## normalize w and y 
    
    if(ntt == 1){
      wki <- wkk
    } else {
      wki <- which(groups == i)
    }
    
    io <- which(wki %in% wo)
    yq[,wki] <- .gjamCompW2Y(yq[,wki,drop=F],notOther=io)$ww
  }
  return( list(wq,yq) )
}

.gjamWcatLoop <- function(y, ws, mus, sgs, notOther, plo, phi, groups, 
                          REDUCT = F){
  
  # if REDUCT, sgs is length-S sigvec
  # if !REDUCT, sgs is css[notOther,notOther]
  
  ntt <- max( groups )
  n <- nrow(y)
  
  for(i in 1:ntt){  
    
    wki <- which(groups == i)
    nki <- length(wki)
    wko <- wki[wki %in% notOther]
    
    w0 <- apply( ws[,wko]*(1 - y[,wko]),1, max) # max(w, 0) for y = 0
    w1 <- apply( ws[,wko]*y[,wko],1, max)       # w for y = 1
    w0[w0 < 0] <- 0                             # when y is reference
    
    si <- sample(wko)
    
    for(s in si){
      
      y1         <- which(y[,s] == 1)
      plo[-y1,s] <- -500
      phi[y1,s]  <- 500
      plo[y1,s]  <- w0[y1]
      phi[-y1,s] <- w1[-y1]
      
      if(REDUCT){
        
        ws[,s] <- .tnorm(n,plo[,s],phi[,s],mus[,s],sqrt(sgs[s]))
        
      } else {
        sm  <- which(notOther == s)
        tmp <- .conditionalMVN(ws[,notOther], mus[,notOther], 
                               sgs, sm)
        mue <- tmp$mu
        vr  <- max(tmp$vr,1e-8)
        ws[,s] <- .tnorm(n,plo[,s],phi[,s],mue,sqrt(vr))
      }
      
      w1[y1]  <- ws[y1,s]    #new w for y = 1
      w0[-y1] <- apply( ws[-y1,wki]*(1 - y[-y1,wki]),1, max)
    }
  }
  list(w = ws, plo = plo, phi = phi)
}

.gjamWcatLoop2 <- function(y, ws, mus, sgs, notOther, plo, phi, groups, 
                           REDUCT = F){
  
  # if REDUCT, sgs is length-S sigvec
  # if !REDUCT, sgs is css[notOther,notOther]
  
  ntt <- max( groups )
  n <- nrow(y)
  
  for(i in 1:ntt){  
    
    wki <- which(groups == i)
    nki <- length(wki)
    wko <- wki[wki %in% notOther]
    
    w1 <- apply( ws[,wko]*y[,wko],1, max)        # w for y = 1
    
    so <- match(wko,notOther)
    
    for(s in wko){   
      
      y1 <- which(y[,s] == 1)
      
      #   if(length(y1) == 0)next
      sm  <- which(notOther == s)   #index in sgs = sg[notOther,notOther]
      sn  <- so[so != sm]           #index in sgs for so
      qs  <- wko[wko != s]          
      
      if(REDUCT){
        ws[y1,s] <- .tnorm(length(y1),plo[y1,s],phi[y1,s],
                           mus[y1,s],sqrt(sgs[s]))
      } else {
        tmp <- .conditionalMVN(ws[y1,notOther], mus[y1,notOther], sgs, sm)
        mue <- tmp$mu
        vr  <- max(tmp$vr,1e-8)
        ws[y1,s] <- .tnorm(length(y1),plo[y1,s],phi[y1,s],mue,sqrt(vr))
      }
      
      w1[y1] <- ws[y1,s]        # w for y = 1
      phi[y1,wki] <- w1[y1]
      phi[y1,s]   <- 500
      
      if(REDUCT){     # the zeros
        tmp <- .tnorm(length(y1)*length(qs),plo[y1,qs],phi[y1,qs],
                      mus[y1,qs],sqrt(sgs[s]))
      } else {
        tmp <- .tnormMVNmatrix(ws[y1,notOther],mus[y1,notOther],
                               smat=sgs, plo[y1,notOther], 
                               hi=phi[y1,notOther],
                               whichSample=so)[,sn,drop=F]
      }
      ws[y1,qs] <- tmp
      
      ###########
      if(length(sn) > 0)tmp <- apply( tmp, 1, max ) #########
      tmp[tmp < 0] <- 0
      plo[y1,s]  <- tmp
    }
    ##############
    s <- wki[!wki %in% wko]   #  y = 1 is ref class
    y1 <- which(y[,s] == 1)
    tmp <- .tnormMVNmatrix(ws[y1,notOther],mus[y1,notOther],
                           smat=sgs, plo[y1,notOther], 
                           hi=phi[y1,notOther],
                           whichSample=so)
    ws[y1,wko] <- tmp[,so]
    #############
  }
  list(w = ws, plo = plo, phi = phi)
}

.gjamWLoop <- function( llist ){
  
  ws <- mus <- sgs <- wkk <- lo <- hi <- sampW <- indexW <- NULL
  byCol <- T
  byRow <- F
  
  llist <- for(k in 1:length(llist))assign( names(llist)[k], llist[[k]] )
  
  n <- nrow(lo)
  tiny <- .00001
  
  if(byCol){
    
    iss <- wkk[wkk %in% indexW]
    
    for(s in iss){
      
      rs <- which(sampW[,s] > 0)
      ls <- lo[drop=F,rs,s]
      hs <- hi[drop=F,rs,s]
      
      tmp <- .conditionalMVN(ws[drop=F,rs,],mus[drop=F,rs,],sgs,s)
      mu  <- tmp$mu
      vr  <- max(tmp$vr,tiny)
      tmp <- .tnorm(length(rs),ls,hs,mu,sqrt(vr))
      
      wl  <- which(tmp == ls)
      if(length(wl) > 0) tmp[wl] <- ls[wl] + tiny*(ls[wl])
      
      wl  <- which(tmp == hs)
      if(length(wl) > 0) tmp[wl] <- hs[wl] - (1 - tiny)*hs[wl]
      
      ws[rs,s] <- tmp
    }
    return(ws)
  }
  
  for(i in indexW){
    
    rs  <- which(sampW[i,] > 0)
    rs  <- rs[rs %in% wkk]
    ws[i,rs] <- .tnormMVNmatrix(ws[drop=F,i,], mus[drop=F,i,],
                                smat=sgs, lo[drop=F,i,], hi[drop=F,i,],
                                whichSample=rs)[,rs]
  }
  ws
}

.setContrasts <- function(xx){
  
  # contrasts where each level is compared to the reference level
  # data must have an attribute for 'reference' class assigned as, e.g.,
  # attr(xdata$soil,'reference') <- 'reference'
  # where xx is xdata$soil and 'reference' is the name that appears in xx
  
  levs  <- attr(xx,'levels') 
  nl    <- length(levs)
  ml    <- nl - 1
  ref <- levs[1]
  
  intType <- attr(xx,'intType')
  
  if(is.null(intType))intType <- 'ref'
  
  wr  <- which(levs == ref)
  
  cj <- matrix(-1/nl,ml,ml)
  diag(cj) <- ml/nl
  rownames(cj) <- levs[-wr]
  colnames(cj) <- levs[-wr]
  
  rj <- rep(-1/nl,ml)
  cj <- rbind(rj,cj)
  rownames(cj)[1] <- ref
  
  levs <- as.character(levs)
  
  cj <- cj[drop=F,levs,]
  if(intType == 'ref'){
    cj[cj > 0] <- 1
    cj[cj < 0] <- 0
  }
  list(levs = levs, cont = cj)
}

.gjamXY <- function(formula, xdata, y, typeNames, notStandard, 
                    checkX = T, xscale = NULL){
  
  n        <- nrow(xdata)
  S        <- ncol(y)
  snames   <- colnames(y)
  facNames <- character(0)
  factorList <- contrast <- NULL
  colnames(xdata) <- .cleanNames(colnames(xdata))
  NOX <- T
  xmean <- 1
  
  original <- colnames(xdata)
  
  xdataNames <- original
  
  if(!is.null(notStandard))notStandard <- .cleanNames(notStandard)
  
  form <- attr( terms(formula), 'term.labels' )
  
  if(length(form) > 0){       # not done if formula = ~ 1
    NOX  <- F
    form <- .cleanNames(form)
    form <- paste0(form,collapse=' + ')
    formula <- as.formula( paste('~',form) )
    
    # no transformation
    t1 <- attr( terms(formula), 'term.labels' )
    wi <- grep('I(',t1,fixed=T)
    if(length(wi) > 0)t1 <- t1[-wi]      # linear terms
    
    wi <- grep(':',t1,fixed=T)
    if(length(wi) > 0)t1 <- t1[-wi]
    
    
    xdata0 <- xdata[,t1, drop=F]
    xnames <- colnames(xdata0)
    
    standX <- !sapply(xdata0,is.factor)
    facNames <- names(standX)[!standX]
    standX   <- names(standX)[standX]
    standX   <- standX[!standX %in% notStandard]
    
    tmp <- .getStandX(xdata0,standX)
    xdata0 <- tmp$xstand
    xmean  <- tmp$xmu
    xsd    <- tmp$xsd
    xscale <- rbind(xmean,xsd)
    
    factorList <- contrast <- vector('list',length = length(facNames))
    names(factorList) <- facNames
    
    if(length(facNames) > 0){
      
      for(j in 1:length(facNames)){
        
        wj <- which(names(xdata0) == facNames[j])
        xf <- as.character(xdata0[[wj]])
        
        cj <- attr(xdata0[[wj]],'contrasts')
        contrast[[j]] <- cj
        tt <- .setContrasts(xdata0[[wj]])$cont
        factorList[[j]] <- paste(facNames[j],colnames(tt),sep='')
        
        if(!is.null(cj))next                       # contrasts previously set
        
        contrast[[j]] <- tt
        attr(xdata0[[wj]],'contrasts') <- tt
      }
      names(contrast) <- facNames
    }  
    
    www <- match(colnames(xdata0),colnames(xdata))
    if(length(www) > 0)xdata[,www] <- xdata0
  }
  
  tmp <- model.frame(formula,data=xdata,na.action=NULL)
  x   <- model.matrix(formula, data=tmp)
  
  colnames(x)[1] <- 'intercept'
  
  xnames <- colnames(x)
  snames <- colnames(y)
  Q      <- ncol(x)
  predXcols <- 2:Q
  isFactor <- character(0)
  
  facBySpec <- missFacSpec <- NULL
  
  VIF <- isNonLinX <- designTable <- NULL
  isInt <- intMat <- numeric(0)
  
  if(!NOX){
    
    if(length(facNames) > 0){
      
      iy <- y*0
      iy[y > 0] <- 1
      facBySpec <- numeric(0)
      missFacSpec <- character(0)
      
      for(j in 1:length(facNames)){
        
        #   ij <- grep(facNames[j],colnames(x))
        
        ij <- which( colnames(x) %in% factorList[[j]] )
        ij <- xnames[ij]
        #  ix <- grep(':',ij)
        #  if(length(ix) > 0)ij <- ij[-ix]
        isFactor <- c(isFactor,ij)
        
        print(paste('observations in factor',facNames[j]))
        print(colSums(x, na.rm=T)[ij])
        
        fs <- matrix(NA,S,length(factorList[[j]]))
        colnames(fs) <- factorList[[j]]
        rownames(fs) <- snames
        
        for(k in 1:length(ij)){
          xi   <- ij[k]
          fs[,k] <- colSums( matrix(x[,xi],n,S)*iy, na.rm=T )
        }
        ms <- 'none missing'
        missFS <- which(fs == 0,arr.ind=T)
        if(length(missFS) > 0){
          ms <- paste(rownames(missFS),ij[missFS[,2]],sep='_')
        }
        
        facBySpec <- append(facBySpec,list(fs))
        missFacSpec <- append(missFacSpec,list(ms))
        
      }
      names(facBySpec) <- names(missFacSpec) <- facNames
    }
    
    
    # check design
    
    if(checkX & length(standX) > 0){
      checkInt <- range(x[,1])
      if(checkInt[1] != 1 | checkInt[2] != 1)
        stop( paste('x[,1] must be intercept (ones)') )
      
      tmp <- .checkDesign(x[,c('intercept',standX)])
      if(tmp$rank < tmp$p)stop( 'x not full rank' )
      VIF         <- tmp$VIF
      designTable <- tmp$designTable$table
    }
    
    if(Q > 2 & length(standX) > 0){
      
      wx <- grep('^2',colnames(x),fixed=T)
      if(length(wx) > 0){
        mm <- unique(unlist(strsplit(colnames(x)[wx],'^2)',fixed=T)))
        mm <- .replaceString(mm,'I(','')
        mm <- match(mm,colnames(x))
        mat <- cbind(wx,mm,mm)
        colnames(mat) <- c('int','main1','main2')
        intMat <- mat
        isInt <- wx
        isNonLinX <- sort(unique( c(isNonLinX,mm,isInt)))
      }
      
      wx <- grep(':',colnames(x))
      if(length(wx) > 0){
        mm  <- matrix(unlist(strsplit(colnames(x)[wx],':')),ncol=2,byrow=T)
        mat <- matrix( match(mm,colnames(x)), ncol=2)
        mat <- cbind(wx,mat)
        colnames(mat) <- c('int','main1','main2')
        wx <- c( which(colnames(x) %in% mm),wx )
        isInt <- sort(c(isInt,wx))
        intMat <- rbind(intMat,mat)
      }
      if(!is.null(isInt))isNonLinX <- sort(unique( c(isNonLinX,isInt)))
    }
    
  }
  
  
  standMat <- matrix(1,Q,1)
  rownames(standMat) <- xnames
  standMu <- standMat - 1
  
  xss <- colnames(xscale)
  
  if(length(xss) > 0){
    standMu[xss,]  <-  xscale['xmean',xss]
    standMat[xss,] <-  xscale['xsd',xss]
  }
  
  # standardize in interactions
  
  if(length(intMat) > 0){
    
    for(j in 1:nrow(intMat)){
      im <- intMat[j,]
      s1 <- s2 <- 1
      if( xnames[im[2]] %in% colnames(xscale) )s1 <- xscale['xsd',xnames[im[2]]]
      if( xnames[im[3]] %in% colnames(xscale) )s2 <- xscale['xsd',xnames[im[3]]]
      
      standMat[im[1],] <- s1*s2
    }
  }
  
  standRows <- which(standMat[,1] != 1 | standMu[,1] != 0)
  standRows <- standRows[!names(standRows) %in% notStandard]
  
  colnames(y) <- .cleanNames(colnames(y))  
  
  # check composition
  
  tiny <- 1 + 1e-10
  
  if('FC' %in% typeNames){
    
    groups <- attr(typeNames,'FCgroups')
    
    if(is.null(groups)){
      groups <- rep(0,S)
      groups[typeNames == 'FC'] <- 1
      attr(typeNames,'FCgroups') <- groups
    }
    
    ngg    <- max(groups)
    for(kk in 1:ngg){
      wf <- which(groups == kk)
      if(length(wf) == 0)stop( 'FC data must have > 1 column' )
      ww <- which(y[,wf] < 0)
      if(length(ww) > 0)stop( 'FC values cannot be < 0' )
      wr <- rowSums(y[,wf],na.rm=T)
      vv <- unique(wr)
      ww <- which(vv != 0 & vv > 1.01)
      if(length(ww) > 0){
        wx <- which(wr %in% vv)
        ii <- paste0(wx, collapse=', ')
        stop( paste('FC data must sum to zero (all absent) or one, check obs:',ii))
      }
    }
  }
  
  if('CC' %in% typeNames){
    wf <- which(typeNames == 'CC')
    if(length(wf) < 2)stop( 'CC data must have > 1 column' )
  }
  
  if(is.null(snames))snames <- paste('S',1:S,sep='-')
  if(is.null(xnames))xnames <- paste('x',1:Q,sep='-')
  
  snames <- sub('_','-',snames)
  xnames <- sub('_','-',xnames)
  
  colnames(y) <- snames
  colnames(x) <- xnames
  
  if(length(isNonLinX) == 0)isNonLinX <- NULL
  if(length(notStandard) == 0)notStandard <- NULL
  
  if( !is.null(notStandard) ){
    ns <- notStandard
    for(k in 1:length(ns)){
      wk <- grep(ns[k],colnames(x))
      ns <- c(ns,colnames(x)[wk])
    }
    notStandard <- unique(ns)
  }
  
  factorAll <- list(nfact = length(factorList), factorList = factorList, 
                    isFactor = isFactor, contrast = contrast,
                    facBySpec = facBySpec, missFacSpec = missFacSpec,
                    facNames  = facNames)
  interaction <- list(isInt = isInt, intMat = intMat, isNonLinX = isNonLinX)
  
  list(x = x, y = y, snames = snames, xnames = xnames, predXcols = predXcols,
       interaction = interaction,factorAll = factorAll,
       xdata = xdata, designTable = designTable, xmean = xmean, xscale = xscale,
       standMu = standMu, standMat = standMat, standRows = standRows,
       notStandard = notStandard, xdataNames = xdataNames, formula = formula)
}

.gjamCompW2Y <- function(ww,notOther=c(1:(ncol(ww)-1))){
  
  pg <- .995
  
  n  <- nrow(ww)
  W  <- rowSums(ww[,notOther,drop=F])
  wh <- which(W > pg)
  other <- c(1:ncol(ww))[-notOther]
  
  if(length(wh) > 0){
    contract <- (1 - (1 - pg)^(W[wh]/pg))/W[wh]
    ww[wh,]  <- ww[wh,]*contract        
  }
  
  ww[,other] <- 1 - rowSums(ww[,notOther,drop=F])
  
  list(pg = pg, ww = ww )
}

.imputX_MVN <- function(xx,yy,beta,xmiss,sinv,xprior=0,xbound=NULL,priorWT=1){
  
  # priorWT is inverse of variance
  
  wcol <- unique(xmiss[,2])
  S    <- nrow(sinv)
  Q    <- nrow(beta)
  
  if(is.null(xbound))xbound <- apply(xx,2,range,na.rm=T)
  
  for(j in wcol){
    
    wx <- which(xmiss[,2] == j)
    
    wj <- xmiss[drop=F,wx,]     # rows, col, missing col j
    wr <- wj[,1]                            # rows missing col j
    xp <- xprior[wx]                    # prior mean
    
    bj <- matrix(beta[j,],1)                # row for missing x
    bn <- matrix(beta[-j,],Q - 1)           # other rows
    
    xn <- xx[drop=F,wr,-j]                  # other cols
    z  <- yy[drop=F,wr,] - xn%*%bn          # y - not missing xb
    datwt <- bj%*%sinv%*%t(bj)              # conditional var
    V     <- 1/( datwt + priorWT*datwt )     
    v     <- z %*%sinv%*%t(bj) + xp*priorWT # conditional 
    xx[wj] <- .tnorm(length(wr),xbound[1,j],xbound[2,j],v%*%V,sqrt(V))
  }
  xx
}

.interp <- function(y,INCREASING=F,minVal=-Inf,maxVal=Inf,defaultValue=NULL,
                    tinySlope=NULL){  #interpolate vector x
  
  if(is.null(defaultValue))defaultValue <- NA
  
  tiny <- .00001
  if(!is.null(tinySlope))tiny <- tinySlope
  
  y[y < minVal] <- minVal
  y[y > maxVal] <- maxVal
  
  n  <- length(y)
  wi <- which(is.finite(y))
  
  if(length(wi) == 0)return(rep(defaultValue,n))
  if(length(wi) == 1)ss <- tiny
  
  xx  <- c(1:n)
  z  <- y
  
  if(wi[1] != 1) wi <- c(1,wi)
  if(max(wi) < n)wi <- c(wi,n)
  
  ss <- diff(z[wi])/diff(xx[wi])
  
  ss[is.na(ss)] <- 0
  
  if(length(ss) > 1){
    if(length(ss) > 2)ss[1] <- ss[2]
    ss[length(ss)] <- ss[length(ss)-1]
  }
  if(INCREASING)ss[ss < tiny] <- tiny
  
  if(is.na(y[1]))  z[1] <- z[wi[2]] - xx[wi[2]]*ss[1]
  if(z[1] < minVal)z[1] <- minVal
  if(z[1] > maxVal)z[1] <- maxVal
  
  for(k in 2:length(wi)){
    
    ki <- c(wi[k-1]:wi[k])
    yk <- z[wi[k-1]] + (xx[ki] - xx[wi[k-1]])*ss[k-1]
    yk[yk < minVal] <- minVal
    yk[yk > maxVal] <- maxVal
    z[ki] <- yk
  }
  z
}

.interpRows <- function(x,startIndex=rep(1,nrow(x)),endIndex=rep(ncol(x),nrow(x)),
                        INCREASING=F,minVal=-Inf,maxVal=Inf,
                        defaultValue=NULL,tinySlope=.001){  
  #interpolate rows of x subject to increasing
  
  nn  <- nrow(x)
  p  <- ncol(x)
  xx <- c(1:p)
  
  if(length(minVal) == 1)minVal <- rep(minVal,nn)
  if(length(maxVal) == 1)maxVal <- rep(maxVal,nn)
  
  ni   <- rep(NA,nn)
  flag <- numeric(0)
  
  z <- x
  
  for(i in 1:nn){
    if(startIndex[i] == endIndex[i]){
      z[i,-startIndex[i]] <- NA
      next
    }
    z[i,startIndex[i]:endIndex[i]] <- .interp(x[i,startIndex[i]:endIndex[i]],
                                              INCREASING,minVal[i],maxVal[i],
                                              defaultValue,tinySlope)
  }
  
  z
}

.invertSigma <- function(sigma,sigmaerror=NULL,otherpar=NULL, REDUCT){
  
  if(REDUCT){
    sinv <- invWbyRcpp(sigmaerror, otherpar$Z[otherpar$K,])
  } else {
    testv <- try( chol(sigma) ,T)
    if( inherits(testv,'try-error') ){
      tiny  <- .1*diag(sigma)
      sigma  <- sigma + diag(diag(sigma + tiny))
      testv <- try(chol(sigma),T)
    }
    sinv    <- chol2inv(testv)
  }
  sinv
}

.invMatZero <- function(sgibbs,nsim=2000,snames,knames,index=NULL,
                        COMPRESS = F, REDUCT=F,
                        sigErrGibbs = NULL, kgibbs = NULL,
                        otherpar = NULL, alpha = .95){   
  # return conditional independence
  # if COMPRESS, sgibbs is as.vector(lower.tri(Sigma,diag=T) )
  # alpha: prob that covariance/inverse is not zero 
  
  S <- length(snames)
  
  if(is.null(index))index <- c(1:nrow(sgibbs))
  simIndex   <- sample(index,nsim,replace=T)
  
  if(!REDUCT){
    
    if(COMPRESS){
      tmp <- .expandSigmaChains(snames, sgibbs, otherpar, 
                                simIndex, sigErrGibbs, kgibbs, 
                                REDUCT=REDUCT, CHAINSONLY=T)$chainList$schain
      sgibbs <- tmp
    }
    S1 <- sqrt(ncol(sgibbs))
  } else {
    N  <- otherpar$N
    r  <- otherpar$r
    S1 <- S
    SS <- matrix(0,S1,S1)
  }
  
  SK     <- length(knames)
  sindex <- match(knames,snames)
  mm     <- matrix(0,SK,SK)
  rownames(mm) <- colnames(mm) <- knames
  hiSS <- loSS <- hiSI <- loSI <- mm
  
  for(j in simIndex){
    if(!REDUCT){
      ss <- matrix(sgibbs[j,],S1,S1) 
      si <- chol2inv(chol( ss ) ) 
      
    } else {
      Z  <- matrix(sgibbs[j,],N,r)
      ss <- .expandSigma(sigErrGibbs[j], S1, Z = Z, kgibbs[j,], REDUCT = T)
      si <- invWbyRcpp(sigErrGibbs[j], Z[kgibbs[j,],])
    }
    
    ss <- ss[sindex,sindex]
    si <- si[sindex,sindex]
    
    hiSS[ss > 0] <- hiSS[ss > 0] + 1/nsim
    loSS[ss < 0] <- loSS[ss < 0] + 1/nsim
    hiSI[si > 0] <- hiSI[si > 0] + 1/nsim
    loSI[si < 0] <- loSI[si < 0] + 1/nsim
  }
  
  loMar <- which(loSS > alpha)
  hiMar <- which(hiSS > alpha)
  inMar <- which(loSS < alpha & hiSS < alpha)   # not different from zero
  
  loCon <- which(loSI > alpha)
  hiCon <- which(hiSI > alpha)
  inCon <- which(loSI < alpha & hiSI < alpha)
  
  inMarMat <- which(loSS < alpha & hiSS < alpha,arr.ind=T)
  inConMat <- which(loSI < alpha & hiSI < alpha,arr.ind=T)
  
  list( inMarMat = inMarMat, inConMat = inConMat )
}

.mapSetup <- function(xlim,ylim,scale=NULL,widex=10.5,widey=6.5){  
  
  #scale is x per inch
  #new means not a new plot
  
  if(is.null(scale))scale <- 1
  
  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  
  if(px > widex){
    dx <- widex/px
    px <- widex
    py <- py*dx
  }
  if(py > widey){
    dx <- widey/py
    py <- widey
    px <- px*dx
  }
  
  par(pin=c(px,py))
  invisible( c(px,py) )
}

.sameByColumn <- function(mat,fraction=F){
  
  nc <- ncol(mat)
  
  sameMat <- matrix(0,nc,nc)
  
  for(j in 2:nc){
    for(k in 1:(j - 1)){
      wj <- which(mat[,j] == mat[,k])
      sameMat[j,k] <- length(wj)
    }
  }
  fraction <- sameMat/nrow(mat)
  fraction[upper.tri(fraction, diag=T)] <- NA
  fraction
}

.modalValuesInArray <- function(mat,idim = 1){
  
  # modal values for each row (idim = 1) or column (idim = 2)
  
  as.numeric( apply(mat,idim,
                    function(x) names(which.max(table(x)))) )
}

.multivarChainNames <- function(rowNames,colNames){
  as.vector( t(outer(colNames,rowNames,paste,sep='_')) )
}

.rMVN <- function (nn, mu, sigma){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.matrix(mu)) mu <- matrix(mu,1)
  if(length(mu) == 1)mu <- matrix(mu,1,nrow(sigma))
  if(ncol(mu) == 1)  mu <- t(mu)
  
  m <- ncol(sigma)
  
  if(ncol(mu) != m)stop('dimension mismatch mu, sigma')
  if(nn > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
  if(nn != nrow(mu))stop('sample size does not match mu')
  
  si <- try(svd(sigma),T)
  
  if( inherits(si,'try-error') ){
    ev <- eigen(sigma, symmetric = TRUE)
    si <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    si <- t(si$v %*% (t(si$u) * sqrt(si$d)))
  }
  p <- matrix(rnorm(nn * m), nn) %*% si
  p + mu
}

.omitChainCol <- function(cmat,omitCols){
  
  #omitCols - characterVector
  
  keep <- c(1:ncol(cmat))
  ocol <- numeric(0)
  for(j in 1:length(omitCols)){
    ocol <- c(ocol,grep(omitCols[j],colnames(cmat)))
  }
  if(length(ocol) > 0)keep <- keep[-ocol]
  list(keep = keep, omit = ocol)
}

.outFile <- function(outFolder=character(0),file){
  paste(outFolder,file,sep='/')
}

.plotLabel <- function(label,location='topleft',cex=1.3,font=1,
                       above=F,below=F,bg=NULL){
  
  if(above){
    adj <- 0
    if(location == 'topright')adj=1
    title(label,adj=adj, font.main = font, font.lab =font)
    return()
  }
  if(below){
    adj <- 0
    if(location == 'bottomright')adj=1
    mtext(label,side=1,adj=adj, outer=F,font.main = font, font.lab =font,cex=cex)
    return()
  }
  
  if(is.null(bg)){
    tmp <- legend(location,legend=' ',bty='n')
  } else {
    tmp <- legend(location,legend=label,bg=bg,border=bg,text.col=bg,bty='o')
  }
  
  xt <- tmp$rect$left # + tmp$rect$w
  yt <- tmp$text$y
  
  pos <- 4
  tmp <- grep('right',location)
  if(length(tmp) > 0)pos <- 2
  
  XX <- par()$xlog
  YY <- par()$ylog
  
  if(XX)xt <- 10^xt
  if(YY)yt <- 10^yt
  
  text(xt,yt,label,cex=cex,font=font,pos=pos)
}

.bins4data <- function(obs, nPerBin=NULL, breaks=NULL, nbin=NULL, LOG=F, POS=T){
  
  if(!is.null(nPerBin)){
    mb <- 20
    if(length(obs)/nPerBin > mb)nperBin <- length(obs)/mb
  }
  
  if( is.null(breaks) ){
    
    if( is.null(nbin) )nbin <- 20
    
    br   <- range(obs[is.finite(obs)],na.rm=T)
    bins <- seq(br[1],br[2],length=nbin)
    if(LOG){
      yy <- obs[obs > 0]
      oo <- min( yy,na.rm=T )
      
      ybin <- seq(log10(oo),log10(max(yy, na.rm=T)),length=20)
      bins <- 10^c(log10(.1*oo),ybin)
      bins <- unique(bins)
      nbin <- length(bins)
      
      nPerBin <- NULL
    }
    
    if( !is.null(nPerBin) ){
      nbb <- nPerBin/length(obs)
      if(nbb < .05)nbb <- .05
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      oo   <- obs
      if(POS)oo <- obs[obs > 0]
      bins <- quantile(oo,nbb,na.rm=T)
      bins <- c(min(oo,na.rm=T),bins)
      bins <- sort(unique(bins))
      
      db <- diff(bins)
      qo <- quantile(obs,c(.1,.9),na.rm=T)
      wb <- which( db/diff(range(qo)) < .02)
      wb <- wb[wb != 1]
      if(length(wb) > 0)bins <- bins[-wb]
      
      nbin <- length(bins)
    }
  } else {
    bins <- breaks
    nbin <- length(bins)
  }
  
  list(breaks = breaks, bins = bins, nbin = nbin)
}

.plotObsPredOld <- function(obs,yMean,ySE=NULL, add=F, box.col='black', opt=NULL){
  
  boxPerc <- .6826895; whiskerPerc <- .95
  nbin <- nPerBin <- breaks <- xlimit <- ylimit <- ptcol <-
    fill <- wide <- NULL
  LOG <- F
  POINTS <- MEDIAN <- T
  xlabel <- 'Observed'; ylabel <- 'Predicted'
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  aa <- (1 - boxPerc)/2
  boxQuant <- c(aa, 1 - aa )
  aa <- (1 - whiskerPerc)/2
  whiskerQuant <- c(aa, 1 - aa )
  
  if(is.null(ptcol)){
    ptcol <- 'black'
  }
  if(length(ptcol) == 1)ptcol <- rep(ptcol,length(obs))
  
  if(is.null(xlimit))xlimit <- quantile(obs[is.finite(obs)],c(.01,.99),na.rm=T)
  if(is.null(ylimit))ylimit <- range(yMean[is.finite(yMean)],na.rm=T)
  
  xxx <- obs
  yyy <- yMean
  
  if(LOG){
    if(is.null(xlimit))xlimit <- range( obs[obs > 0],na.rm=T )
    if(is.null(ylimit))ylimit <- range( yMean[yMean > 0],na.rm=T )
    if(xlimit[1] <= 0)xlimit[1] <- .001
  }
  
  if(!POINTS){
    xxx <- xlimit[1]
    yyy <- ylimit[1]
  }
  
  if(!add){
    if(is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel,log='xy')
    }
    if(!is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel,
                   xlim=xlimit,ylim=ylimit)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.03,xlab=xlabel,ylab=ylabel,
                   xlim=xlimit,log='xy',ylim=ylimit)
    }
  }
  
  if(POINTS)points(xxx,yyy,pch=16,col=.getColor(ptcol,.5), cex=.5)
  
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),
                                 col='grey',lwd=2)
  }
  
  tmp    <- .bins4data(obs,nPerBin=nPerBin,breaks=breaks,LOG=LOG)
  breaks <- tmp$breaks
  bins   <- tmp$bins
  nbin   <- tmp$nbin
  
  if(is.null(wide))wide <- diff(bins)/2.1
  if(length(wide) == 1)wide <- rep(wide,nbin)
  minmax <- par('usr')[1:2]
  dff    <- diff(minmax)
  if(!LOG)wide[wide > dff/5] <- dff/5
  
  maxx <- 0
  last <- F
  
  for(k in 1:(nbin-1)){
    
    mb <- bins[k+1]
    if(mb >= xlimit[2]){
      last <- T
      mb   <- Inf
    }
    ok <- which(obs >= bins[k] & obs < mb)
    if(length(ok) == 0)next
    qk <- which(is.finite(yMean) & obs >= bins[k] & obs <= mb)
    q  <- quantile(yMean[qk],c(.5,whiskerQuant[1],boxQuant[1],
                               boxQuant[2],whiskerQuant[2]),na.rm=T)
    if(LOG)q[q <= 0] <- ylimit[1]
    
    ym <- q[1]
    xx <- mean(bins[k:(k+1)])      # bounded by bins
    if(!LOG){
      if(MEDIAN)xx <- median(obs[ok],na.rm=T)
    } else {
      xx <-  sqrt( prod(bins[k:(k+1)]) )
    }
    points(xx,q[1],pch=3,col=box.col)
    yy    <- q[c(2,5)]
    yy[1] <- max(c(yy[1],ylimit[1]),na.rm=T) + .0000001
    yy[2] <- max(yy)
    
    yy1 <- q[3]
    yy1 <- max(yy1,ylimit[1],na.rm=T) + .00000001
    yy2 <- max(yy1,q[4])
    
    minx <- xx - .3*(xx - bins[k])
    maxx <- xx + .3*(mb - xx)
    
    dx1 <- xx - minx
    dx2 <- maxx - xx
    if(dx1 > dx2)minx <- xx - dx2
    if(dx1 < dx2)maxx <- xx + dx1
    
    figRange <- par('usr')
    totalwide <- (maxx - minx)/diff(figRange[1:2])
    
    if(is.null(nPerBin)){
      
      if(maxx >= xlimit[2])maxx <- xlimit[2]
      
      if(LOG & k == 1){
        
        if(xx == 0)xx <- .5*bins[k+1]
        
        dx   <- log10(bins[k+1]) - log10(xx)
        maxx <- 10^(log10(xx) + .2*dx)
        if(k == 1){
          dx <- -log10(xlimit[1]) + log10(xx)
        } else {
          dx <- -log10(bins[k-1]) + log10(xx)
        }
        minx <- 10^(log10(xx) - .2*dx)
        if(minx < xlimit[1])minx <- xlimit[1]
        totalwide <- (log10(maxx) - log10(minx))/diff(figRange[1:2])
      }
      
      
      rect(minx,yy1,maxx,yy2,col=fill,border=box.col)
      lines(c(minx,maxx),c(ym,ym),lwd=2,col=box.col)
    }
    
    if(!is.null(nPerBin)){
      qo <- quantile(obs[ok],c(.3,.7,.25,.75),na.rm=T)
      if(qo[1] == qo[2] | !MEDIAN)qo <- c(xx-.2*wide[k],
                                          xx+.2*wide[k],xx-.3*wide[k],
                                          xx+.3*wide[k])
      rect(qo[1],yy1,qo[2],yy2,col=fill,border=box.col)
      lines(c(qo[3],qo[4]),c(ym,ym),lwd=2,col=box.col)
      lines(rep(mean(qo[1:2]),2),yy,lwd=2,col=box.col)
    } else {
      lines(c(xx,xx),yy,lwd=2,col=box.col)
    }
    if(last)break
  }
  
  invisible( bins )
}

.predictY2X_linear <- function(xpred,yy,bb,ss,sinv=NULL, 
                               priorIV = diag(1e-10,ncol(xpred)), 
                               priorX = matrix(0,ncol(xpred)), 
                               predCols = c(2:ncol(xpred)),REDUCT, lox, hix){
  
  #inverse prediction for multivariate linear in x
  
  prX <- priorX[predCols]
  if(!is.matrix(prX))prX <- matrix(prX)
  
  nn <- nrow(yy)
  notPred <- c(1:ncol(xpred))[-predCols]
  
  bp <- matrix(bb[drop=F,predCols,],length(predCols))
  
  if(length(notPred) > 0){
    bn <- matrix(bb[notPred,],length(notPred))
    yy <- yy - xpred[,notPred]%*%bn
  }
  pp <- length(predCols)
  
  if(is.null(sinv))sinv <- chol2inv(chol(ss))
  
  bs <- bp%*%sinv
  
  V <- chol2inv(chol( bs%*%t(bp) + priorIV[predCols,predCols] ) )
  v <- yy%*%t(bs) + matrix( priorIV[predCols,predCols] %*% prX,nn,pp,byrow=T)
  mu <- v%*%V
  
  qq <- ncol(mu)
  
  if(qq > 1){
    xpred[,predCols] <- .tnormMVNmatrix(avec=xpred[,predCols],muvec=mu,smat=V,
                                        lo=matrix(lox[predCols],nn,qq,byrow=T),
                                        hi=matrix(hix[predCols],nn,qq,byrow=T))
  } else {
    xpred[,predCols] <- .tnorm(nn,lox[predCols],hix[predCols], mu,sqrt(V))
  }
  xpred
}

.predictY2X_nonLinear <- function(xx,yy,bb,ss,priorIV = diag(1e-10,ncol(xx)), 
                                  priorX=matrix(0,ncol(xx)),
                                  factorObject, interObject, lox, hix){
  
  #inverse prediction for multivariate nonlinear in x and factors, metropolis
  
  predCols <- interObject$isNonLinX
  isInt    <- interObject$isInt
  intMat   <- interObject$intMat
  isFactor <- factorObject$isFactor
  factorList <- factorObject$factorList
  contrast  <- factorObject$contrast
  
  iFcol  <- NULL
  priorX <- priorX[predCols]
  if(!is.matrix(priorX))priorX <- matrix(priorX)
  
  nn <- nrow(yy)
  intercept <- xx[,1]
  
  xnew <- xx
  
  xv <- as.vector(xx[,predCols])
  nv <- length(xv)
  lo <- rep(lox[predCols],each=nn)
  hi <- rep(hix[predCols],each=nn)
  xnew[,predCols] <- .tnorm(nv,lo,hi,xv,.01)
  
  if(length(isFactor) > 0){          # all factors, main effects
    np <- length(factorList)
    for(k in 1:np){
      nf <- length(factorList[[k]]) + 1
      tm <- contrast[[k]][sample(nf,nn,replace=T),]
      xnew[,factorList[[k]]] <- tm
    }
    iFcol <- match(isFactor,colnames(xx))
  }
  
  if(length(intMat) > 0){     # some of the nlin terms interactions?
    xnew[,intMat[,1]] <- xnew[,intMat[,2]]*xnew[,intMat[,3]]
  }
  
  pnow <- .dMVN(yy,xx%*%bb,ss,log=T)
  pnew <- .dMVN(yy,xnew%*%bb,smat=ss,log=T)
  
  a  <- exp(pnew - pnow)
  z  <- runif(nn,0,1)
  wa <- which(z < a)
  xx[wa,] <- xnew[wa,]
  
  list(x = xx, accept = length(wa))
}

.predVsObs <- function(true,p,xlim=range(true),ylim=range(p,na.rm=T),xlab=' ',
                       ylab=' ', colors=rep(1,length(true)),lwd=2,add=F){ 
  
  #true  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  if(!is.matrix(p))p <- matrix(p,ncol=1)
  
  nn <- length(true)
  y  <- apply(p,2,quantile,c(.5,.025,.975))
  
  if(!add)plot(true,y[1,],xlim=xlim,ylim=ylim,xlab=xlab,
               ylab=ylab,col=colors,pch=3,lwd=lwd)
  points(true,y[1,],col=colors,pch=3,lwd=lwd)
  
  for(j in 1:nn)lines(c(true[j],true[j]),y[2:3,j],col=colors[j],lwd=lwd)
  abline(0,1,lty=2)
  
  invisible(y)
}

.processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                         sigOnly = F,burnin=1,xlimits = NULL){  
  
  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  NOPARS <- F
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
      btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }
    
    wq   <- apply(btmp,2,quantile,c(.025,.975),na.rm=T)  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    
    if(length(wq) == ncol(btmp))NOPARS <- T
    if(NOPARS) warning('no significant pars to plot')
    if(length(wq) > 0 & !NOPARS){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
  }
  
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
    if(burnin > (nrow(xgb) + 100))stop("burnin too large")
    xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)
  
  out <- t(rbind(apply(xgb,2,mean,na.rm=T),apply(xgb,2,sd,na.rm=T),
                 apply(xgb,2,quantile,c(.025,.975),na.rm=T)))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('estimate','se','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('estimate','se','0.025','0.975','true value')
  }
  
  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf),mar=c(4,2,2,2))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc),mar=c(4,2,2,2))
  
  if(CPLOT & !NOPARS){
    for(j in 1:nc){
      plot(xgb[,j],type='l')
      abline(h=out[j,],lty=2)
      if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
      abline(h = 0, col='grey',lwd=2)
      title(colnames(xgb)[j])
    }
  }
  xlims <- xlimits
  if(DPLOT & !NOPARS){
    for(j in 1:nc){
      xj <- density(xgb[,j])
      if(is.null(xlimits))xlims <- range(xj$x)
      plot(xj$x,xj$y,type='l',xlim=xlims)
      abline(v=out[j,],lty=2)
      if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
      title(colnames(xgb)[j])
    }
  }
  list(summary = signif(out,4))
}

.replaceString <- function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx,fixed=T)
  if(length(ww) == 0)return(xx)
  
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now,fixed=T) )
    ss <- s[1]
    if(length(s) == 1)ss <- paste( ss,new,sep='')
    if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}

.cleanNames <- function(xx){
  
  xx <- .replaceString(xx,'-','')
  xx <- .replaceString(xx,'_','')
  xx <- .replaceString(xx,' ','')
  xx <- .replaceString(xx,"'",'')
  
  xx
}

.buildYdata <- function(ydata, ytypes){
  
  # when y has factors, data.frame to matrix
  
  S <- ncol(ydata)
  
  wd <- which(duplicated(colnames(ydata)))
  if(length(wd) > 0){
    warning('duplicated colummn names in ydata')
    for(k in 1:length(wd)){
      dname <- colnames(ydata)[wd[k]]
      wk    <- which(colnames(ydata) == dname)
      colnames(ydata)[wk] <- paste(dname,1:length(wk),sep='')
    }
  }
  
  original <- colnames(ydata)
  colnames(ydata) <- .cleanNames(colnames(ydata))
  
  new <- colnames(ydata)
  ydataNames <- rbind(original,new)
  
  CCgroups  <- attr(ytypes,'CCgroups')
  FCgroups  <- attr(ytypes,'FCgroups')
  CATgroups <- attr(ytypes,'CATgroups')
  
  ngroup  <- 0
  ccg     <- CCgroups
  fcg     <- FCgroups
  y       <- numeric(0)
  
  snames <- colnames(ydata)
  nc     <- ncol(ydata)
  wfact  <- .whichFactor(ydata)
  nfact  <- length(wfact)
  wnot   <- c(1:nc)
  if(nfact > 0)wnot <- wnot[-wfact]
  
  ntypes <- character(0)
  
  if(length(wnot) > 0){
    
    if(is.null(ccg)) ccg <- rep(0,length(wnot)) # if not assigned, assume same 
    if(is.null(fcg)) fcg <- rep(0,length(wnot))
    
    snames <- snames[wnot]
    ntypes <- ytypes[wnot]
    y      <- ydata[,wnot,drop=F]
    
    wcomp <- grep('CC',ytypes[wnot])
    ncomp <- length(wcomp)
    if(ncomp > 0){
      if( max(ccg[wnot[wcomp]]) == 0 )ccg[wnot[wcomp]] <- 1  #assume same group
      goo <- grep('other',snames[wcomp])
      if( length(goo) == 0 )snames[wcomp[ncomp]] <- 
        paste(snames[wcomp[ncomp]],'other',sep='')
    }
    
    wcomp <- grep('FC',ytypes[wnot])
    ncomp <- length(wcomp)
    if(ncomp > 0){
      if( max(fcg[wnot[wcomp]]) == 0)fcg[wnot[wcomp]] <- 1  #assume same group
      goo <- grep('other',snames[wcomp])
      if(length(goo) == 0)snames[wcomp[ncomp]] <- 
        paste(snames[wcomp[ncomp]],'other',sep='')
    }
  }
  
  if(nfact > 0){   # categorical
    
    ngroup <- 0
    ycat   <- cag <- numeric(0)
    if(length(ccg) > 0)cag    <- ccg*0
    
    for(j in 1:nfact){
      
      ngroup <- ngroup + 1
      
      conj   <- contrasts(ydata[,wfact[j]],contrasts=F)
      cj     <- colnames(conj)
      yj     <- conj[ydata[,wfact[j]],]
      colnames(yj) <- paste(colnames(ydata)[wfact[j]],cj,sep='')
      
      w11    <- which(colSums(yj) > 0)  #drop empty levels
      yj     <- yj[,w11]
      cj     <- cj[w11]
      
      goo <- grep('other',colnames(yj))
      if(length(goo) == 0){
        colnames(yj)[ncol(yj)] <- paste(colnames(ydata)[wfact[j]],'other',sep='')
        cj[ncol(yj)] <- colnames(yj)[ncol(yj)]
      }
      
      ycat   <- cbind(ycat, yj)
      cag    <- c(cag,rep(ngroup,length(cj)))
      fcg    <- c(fcg,rep(0,length(cj)))
      ccg    <- c(ccg,rep(0,length(cj)))
      ntypes <- c(ntypes,rep('CAT',length(cj)))
    }
    
    rownames(ycat) <- NULL
    n1 <- ncol(y) + 1
    n2 <- ncol(ycat)
    y <- cbind(y,ycat)
    attr(ntypes,'CATgroups')  <- cag
  }
  
  if(max(ccg) > 0)attr(ntypes,'CCgroups') <- ccg
  if(max(fcg) > 0)attr(ntypes,'FCgroups') <- fcg
  
  list(y = as.matrix(y), CCgroups = ccg, FCgroups = fcg, 
       CATgroups = attr(ntypes,'CATgroups'), typeNames = ntypes,
       ydataNames = ydataNames)
}

.setUpSim <- function(n, S, Q, x, typeNames){
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  
  notOther <- c(1:S)
  snames   <- character(0)
  tnames   <- character(0)
  sN       <- S
  catCols  <- NULL
  
  ngroup <- fgroup <- cgroup <- 1
  GROUPS <- F
  CCgroups <- FCgroups <- CATgroups <- numeric(0)
  s      <- 0
  
  wcc <- which(!typeNames %in% c('CC','FC','CAT'))
  ncc <- length(wcc)
  if(ncc > 0){
    snames   <- paste('S',c(1:ncc),sep='')
    CCgroups <- FCgroups <- CATgroups <- rep(0,ncc)
    tnames   <- typeNames[wcc]
    s <- ncc
  }
  wcc <- which(typeNames == 'CC')
  ncc <- length(wcc)
  if(ncc > 0){
    ss <- c( (s+1):(s+ncc))
    CCgroups <- c(CCgroups,rep(1,ncc))
    FCgroups <- c(FCgroups,rep(0,ncc))
    tnames   <- c(tnames,rep('CC',ncc))
    snn      <- paste('S',ss,sep='')
    snn[ncc] <- paste(snn[ncc],'other',sep='')
    snames   <- c(snames, snn)
    ngroup <- 1
    s <- max(ss)
  }
  wcc      <- which(typeNames == 'FC')
  ncc <- length(wcc)
  if(ncc > 0){
    ss <- c( (s+1):(s+ncc))
    FCgroups <- c(FCgroups,rep(1,ncc))
    CCgroups <- c(CCgroups,rep(0,ncc))
    tnames   <- c(tnames,rep('FC',ncc))
    snn      <- paste('S',ss,sep='')
    snn[ncc] <- paste(snn[ncc],'other',sep='')
    snames   <- c(snames, snn)
    fgroup <- 1
    s <- max(ss)
  }
  
  CATgroups <- CCgroups*0
  
  if( 'CAT' %in% typeNames ){
    
    wk     <- which(typeNames == 'CAT')
    ncomp  <- length(wk)
    ncat   <- sample(3:4,ncomp,replace=T)
    nall   <- sum(ncat)
    ntot   <- s + nall
    CATgroups <- rep(0,s)
    
    js <- s
    for(j in 1:ncomp){
      js     <- js + 1
      sseq   <- (s+1):(s + ncat[j])
      cj <- paste('S',js,letters[1:ncat[j]],sep='')
      cj[ncat[j]] <- paste('S',js,'other',sep='')
      snames    <- c(snames,cj)
      CATgroups <- c(CATgroups,rep(j,ncat[j]))
      tnames    <- c(tnames,rep('CAT',ncat[j]))
      s <- max(sseq)
    }
    CCgroups <- c(CCgroups,rep(0,sum(ncat)))
    FCgroups <- c(FCgroups,rep(0,sum(ncat)))     
    catCols  <- which(CATgroups > 0)
    cgroup   <- ncomp
  }
  sN     <- length(tnames)
  oo     <- grep('other',snames)
  notOther <- c(1:sN)[-oo]
  
  tmp <- .gjamGetTypes(tnames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  typeNames <- tmp$typeNames
  
  if(is.null(x)){
    x <- matrix( rnorm(n*Q,.1), n, Q)  
    x[,1] <- 1
  }
  
  beta <- matrix(0, Q, sN)
  ss   <- diag(.01,sN)    
  
  colnames(beta) <- colnames(ss) <- rownames(ss) <- snames
  wkeep <- numeric(0)
  cnames <- tnames <- character(0)
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    
    if( typeFull[wk[1]] == 'presenceAbsence' ){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(Q*nk,-1.5,1.5)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    if(typeFull[wk[1]] %in% c('continuous','contAbun')){
      diag(ss)[wk] <- .4
      beta[,wk]    <- runif(Q*nk,-.5,2)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    if(typeFull[wk[1]] == 'discAbun'){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(Q*nk,-.1,2)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    if(typeFull[wk[1]] == 'ordinal'){
      diag(ss)[wk] <- 1
      beta[,wk]    <- runif(Q*nk,-.4,2)
      wkeep <- c(wkeep,wk)
      tnames <- c(tnames,typeNames[wk])
      cnames <- c(cnames,colnames(beta)[wk])
    }
    
    if( typeFull[wk[1]] %in% c('fracComp','countComp','categorical') ){
      
      if(length(wk) < 2)stop('composition data must have at least 2 columns')
      
      ntt <- cgroup
      if( typeFull[wk[1]] == 'fracComp' ){
        ntt <- fgroup
        attr(tnames,'FCgroups') <- FCgroups
      }
      if( typeFull[wk[1]] == 'countComp' ){
        ntt <- ngroup
        attr(tnames,'CCgroups') <- CCgroups
      }
      if( typeFull[wk[1]] == 'categorical' ){
        attr(tnames,'CATgroups') <- CATgroups
      }
      
      for(i in 1:ntt){
        
        if(ntt == 1){   
          wki <- wk
        } else {
          if( typeFull[wk[1]] == 'countComp' )wki <- 
              which(typeCols == k & CCgroups == i)
          if( typeFull[wk[1]] == 'fracComp' )wki <- 
              which(typeCols == k & FCgroups == i)
          if( typeFull[wk[1]] == 'categorical' )wki <- 
              which(typeCols == k & CATgroups == i)
        }
        
        nki    <- length(wki)
        
        if( typeFull[wk[1]] == 'categorical' ){
          
          bb <- matrix( rnorm(Q*nki,0,.5), Q,nki)
          bb[1,] <- bb[1,]*0
          
          for(kk in 1:5){
            
            mu   <- x%*%bb
            w    <- mu
            cols <- apply(w,1,which.max)
            mindex <- cbind( c(1:n),cols )
            
            wmax <- w[mindex]
            ww   <- which(wmax < 0)
            nw   <- length(ww)
            if(nw > 0) w[mindex[ww,]] <- .tnorm(nw,0,10,mu[mindex[ww,]],1)
            
            bb <- solveRcpp(crossprod(x))%*%crossprod(x,w)
          }
          
          keep <- as.numeric( names(table(cols))  )
          wkeep <- c(wkeep,wki[keep])
          tnames <- c(tnames,rep('CAT',length(keep)))
          
          bbb  <- colnames(beta)[wki[keep]]
          if(length(keep) < nki){
            bbb  <- substr(bbb,1,2)
            labs <- c(letters[1:(length(bbb) - 1)],'other')
            bbb  <- paste(bbb,labs,sep='')
          }
          cnames <- c(cnames,bbb)
          beta[,wki] <- bb
          diag(ss)[wk] <- 1
          
        } else {
          
          bb <- matrix( rnorm(Q*nki,0,1/nki), Q, nki)
          bb[1,] <- bb[1,]*0
          
          w  <- x%*%bb
          
          for(m in 1:3){
            w1 <- w
            w1[w < 0] <- 0
            w2 <- sweep(w1,1,rowSums(w1),'/')
            w[w >= 0] <- w2[w >= 0]
            bb <- solveRcpp(crossprod(x))%*%crossprod(x,w)
            w  <- x%*%bb
          }
          
          wkeep <- c(wkeep,wki)
          tnames <- c(tnames,typeNames[wki])
          cnames <- c(cnames,colnames(beta)[wki])
          diag(ss)[wk] <- .1/nk^2.5
          beta[,wki] <- bb
        }
        
      }
    }
    
  }
  
  S <- length(wkeep)
  beta      <- beta[,wkeep]
  sigma     <- ss[wkeep,wkeep]
  colnames(beta) <- colnames(sigma) <- rownames(sigma) <- cnames
  CCgroups  <-  CCgroups[wkeep] 
  FCgroups  <-  FCgroups[wkeep]
  CATgroups <-  CATgroups[wkeep]
  snames    <- cnames
  other     <- numeric(0)
  notOther  <- c(1:S)
  other     <- grep('other',snames)
  if(length(other) > 0)notOther <- notOther[-other]
  
  list(beta = beta, x = x, sigma = sigma, CCgroups = CCgroups, 
       FCgroups = FCgroups, CATgroups = CATgroups, typeNames = tnames,
       other = other, notOther = notOther, snames = snames)
}

.between <- function(x,lo,hi,ILO = T, IHI = T, OUT=F){
  
  if(length(x) == 0) return( numeric(0) )
  
  if(OUT)return( which(x < lo | x > hi) )
  if(!ILO & !IHI ) return( which(x > lo & x < hi) )
  if(!ILO &  IHI ) return( which(x > lo & x <= hi) )
  if( ILO & !IHI ) return( which(x >= lo & x < hi) )
  if( ILO &  IHI ) return( which(x >= lo & x <= hi) )
  
}        

.simData <- function( n, S, Q, x, typeNames, nmiss, effort ){
  
  #  pg <- .95
  
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  
  typeNotCat <- typeNames
  
  cgrep <- grep('CAT',typeNames)
  if(length(cgrep) > 0){
    ycat <- vector( mode = 'list', length=length(cgrep) )
    names(ycat) <- paste('CAT',1:length(cgrep),sep='_')
  }
  
  cuts <- numeric(0)
  
  tmp    <- .setUpSim(n, S, Q, x, typeNames)
  beta   <- tmp$beta
  x      <- tmp$x
  sig    <- tmp$sigma
  snames <- colnames(beta)
  typeNames <- tmp$typeNames
  other     <- tmp$other
  notOther  <- tmp$notOther
  CCgroups  <- tmp$CCgroups
  FCgroups  <- tmp$FCgroups
  CATgroups <- tmp$CATgroups
  
  tmp <- .gjamGetTypes(typeNames)
  typeCols  <- tmp$typeCols
  typeFull  <- tmp$typeFull
  typeCode  <- tmp$TYPES[typeCols]
  allTypes  <- sort(unique(typeCols))
  
  S      <- length(typeNames)
  xnames <- paste('x',1:Q,sep='')
  
  SS    <- matrix(1,S,S)
  SS[lower.tri(SS)] <- runif(S*(S - 1)/2,-.98,.98)
  SS[upper.tri(SS)] <- SS[lower.tri(SS)]
  
  SS    <- cor( .rMVN(S+5,0,SS) )
  SS    <- .cor2Cov(diag(sig),SS)
  
  sigma <- .rwish(S+2,SS)/(S + 2)
  
  corCols <- which(typeNames %in% c('PA','OC','CAT'))
  if(length(corCols) > 0){
    corSpec <- .cov2Cor(sigma)
    sigma[corCols,corCols] <- corSpec[corCols,corCols]
  }
  
  beta[,other] <- 0
  mu <- w <- matrix(0,n,S)
  
  mu[,notOther] <- x%*%beta[,notOther]
  w[,notOther]  <- mu[,notOther] + .rMVN(n,0,sigma[notOther,notOther]) 
  colnames(w) <- snames
  
  y  <- w
  z  <- w*0
  z[w <= 0]   <- 1
  z[w > 0]    <- 2
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk) 
    
    if( typeFull[wk[1]] %in% c('fracComp','countComp','categorical') ){
      
      if( typeFull[wk[1]] == 'fracComp' )
        groups <- attr(typeNames,'FCgroups') <- FCgroups
      if( typeFull[wk[1]] == 'countComp' )
        groups <- attr(typeNames,'CCgroups') <- CCgroups
      if( typeFull[wk[1]] == 'categorical' )
        groups <- attr(typeNames,'CATgroups') <- CATgroups
      
      ntt <- max(c(1,groups))
      
      for(i in 1:ntt){
        
        if(ntt == 1){
          wki <- wk
        } else {
          wki <- which(typeCols == k & groups == i)
        }
        nki <- length(wki)
        
        if( typeFull[wk[1]] == 'categorical' ){
          
          wko  <- wki[1:(nki-1)]                  
          wcol <- apply(w[,wko],1,which.max)
          w0   <- which( w[,wko][ cbind( c(1:n),wcol ) ] < 0 )
          if(length(w0) > 0)wcol[w0] <- nki
          
          wtab <- tabulate(wcol)
          if(length(wtab) < nki){
            ww <- rep(0,nki)
            ww[1:length(wtab)] <- wtab
            wtab <- ww
          }
          
          if(min(wtab) < 5){
            wlo <- which(wtab < 5)
            for(s in 1:length(wlo)){
              wro <- sample(n,5)
              wcol[wro] <- wlo[s]
              tmp <- w[wro,wki]
              if(wlo[s] == nki){
                tmp[tmp > -.01] <- -.01   # all values neg
                tmp[,nki] <- .1
              } else {
                mm <- pmax(0,apply(tmp,1,max))
                tmp[,wlo[s]] <- mm + .1
              }
              w[wro,wki] <- tmp
            }
          }
          
          mindex <- cbind(1:n,wcol)
          
          vv <- colnames(w)[wki[wcol]]
          mm <- nchar(vv)
          vv <- substr(vv,3,mm)                  
          
          ycat[[i]] <- vv
          
          yk   <- w[,wki]*0
          yk[ mindex ] <- 1
          y[,wki] <- yk
          z[,wki] <- yk + 1
          
        } else {
          
          noto <- c(1:nki)[-grep('other',snames[wki])]
          
          ww     <- w[,wki]
          
          for(j in 1:5){
            
            w0     <- which(ww < 0)
            ww[w0] <- 0  
            
            yk  <- .gjamCompW2Y(ww,notOther=noto)$ww
            
            yplus <- which(yk > 0)
            yminu <- which(yk < 0)
            
            ww[yplus] <- yk[yplus]
            
            bb <- solveRcpp(crossprod(x))%*%crossprod(x,ww)
            mu <- x%*%bb
            ww <- mu + .rMVN(n,0,sigma)[,wki]
          }
          zk     <- ww*0 + 1
          zk[w0] <- 0
          w[,wki] <- ww  
          beta[,wki] <- bb
          
          if(typeFull[wk[1]] == 'fracComp'){
            y[,wki] <- yk
            z[,wki] <- zk
          }
          if( typeFull[wk[1]] == 'countComp' ){
            
            mm <- S*20
            a  <- 4
            b  <- mm/a
            ee <- rpois(n,rgamma(n,shape=a,scale=b))
            yy <- sweep(yk,1,ee,'*')
            
            ww <- ceiling(yy)
            ww[ww < 0] <- 0
            y[,wki] <- ww
            z[,wki] <- ww + 1
          }
        }
      }
    }
    
    if( typeFull[wk[1]] != 'continuous' ) y[,wk][y[,wk] < 0] <- 0  # not cens
    
    if( typeFull[wk[1]] == 'presenceAbsence' )y[,wk] <- z[,wk] - 1       
    
    if( typeFull[wk[1]] == 'discAbun' ){
      
      if(!is.null(effort)){
        we     <- wk[wk %in% effort$columns]
        y[,we] <- round( w[,we]*effort$values,0 )
      } else {
        w0 <- round(w[,wk,drop=F],0)
        y[,wk] <- w0
      }
      y[,wk][y[,wk] < 0] <- 0
      z[,wk]        <- y[,wk] + 1
    }
    
    if( typeFull[wk[1]] == 'ordinal' ){
      
      yy   <- w[,wk,drop=F]
      ncut <- 8
      maxw <- floor(max(yy))
      
      cuts  <- t( matrix( c(-Inf, seq(0,(maxw-1),length=(ncut-2)) ,Inf),
                          ncut,nk) )
      rownames(cuts) <- snames[wk]
      
      for(j in 1:nk){
        z[,wk[j]]   <- findInterval(yy[,j],cuts[j,])
      }
      
      y[,wk] <- z[,wk] - 1
    }
    
  }
  
  #####################################
  
  noMore <- F
  if( 'categorical' %in% typeFull & noMore){
    
    wss  <- w*0
    wss[,notOther] <- .sqrtRootMatrix(w[,notOther],sigma[notOther,notOther],
                                      DIVIDE=T)
    css  <- .cov2Cor(sigma[notOther,notOther])
    alpha <- .sqrtRootMatrix(beta,sigma, DIVIDE=T)
    muss <- x%*%alpha
    
    wk <- which(typeNames == 'CAT')
    wo <- which(wk %in% notOther)
    
    plo  <- w*0 - 500
    phi  <- w*0 + 500
    
    phi[y == 0] <- 0
    plo[y == 1] <- w[y == 1]
    IXX <- solveRcpp(crossprod(x))
    
    for(k in 1:25){
      
      tmp <- .gjamWcatLoop2(y, ws = wss, mus = muss, sgs = css, 
                            notOther, plo, phi, groups = CATgroups)
      wss[,wk] <- tmp$w[,wk]
      plo     <- tmp$plo
      phi     <- tmp$phi
      beta[,wo] <- IXX%*%crossprod(x,wss[,wo])
      muss[,wo] <- x%*%beta[,wo]
    }
    w[,wo] <- wss[,wo]
  }
  
  
  beta <- solveRcpp(crossprod(x))%*%crossprod(x,w)
  sigma[notOther,notOther] <- var(w[,notOther] - x%*%beta[,notOther]) ### NO
  sigma[other,] <- sigma[,other] <- 0
  diag(sigma)[other] <- diag(sig)[other]
  
  ydata <- data.frame(y)
  typeFrame <- typeNames
  
  if('CAT' %in% typeNames){
    wcat <- grep('CAT',typeNames)
    wnot <- c(1:S)[-wcat] 
    nss  <- length(wnot) + 1
    ncc  <- length(wnot) + length(ycat)
    names(ycat) <- paste('S',nss:ncc,sep='')
    ydata <- as.data.frame(ycat)
    if(length(wnot) > 0)ydata <- cbind(y[,wnot,drop=F],ydata)
    typeFrame <- c(typeNames[wnot], rep('CAT',length(ycat)))
  }
  
  if(nmiss > 0){
    x[ sample(length(x),nmiss) ] <- NA
    x[,1] <- 1
    wmiss <- which(is.na(x),arr.ind=T)
    nmiss <- nrow(wmiss)
  }
  
  xnames[1]      <- 'intercept'
  colnames(y)    <- snames
  colnames(beta) <- rownames(sigma) <- colnames(sigma) <- snames
  colnames(x)    <- rownames(beta) <- xnames
  
  form <- as.formula( paste('~ ',paste(colnames(x)[-1],collapse='+' )) )
  
  list(formula = form, xdata = data.frame(x), ydata = ydata,
       y = y, w = w,  typeNames = typeFrame, typeY = typeNames, effort = effort,
       trueValues = list(beta = beta, sigma = sigma, 
                         corSpec = .cov2Cor(sigma), cuts = cuts))
}

.tnorm <- function(n,lo,hi,mu,sig){   
  
  #normal truncated lo and hi
  
  tiny <- 10e-6
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] - tiny
  z
}

.traitLabel <- function(tname){
  
  tname <- .replaceString(tname,now='soilFactor',new='')
  tname[tname == 'gmPerSeed'] <- 'Seed mass'
  tname[tname == 'gmPerCm']   <- 'Wood dens'
  tname[tname == 'woodSG']    <- 'Wood dens (green)'
  tname[tname == 'maxHt']     <- 'Max ht'
  tname[tname == 'leafN']     <- 'leaf [N]'
  tname[tname == 'leafP']     <- 'leaf [P]'
  tname[tname == "other"]  <- 'Deciduous'
  tname[tname == "broaddeciduous"]  <- 'Deciduous'
  tname[tname == "broadevergreen"]  <- 'BL evergrn'
  tname[tname == "needleevergreen"] <- 'NL evergrn'
  tname[tname == "dioecious"] <- 'Dioecious'
  tname[tname == "u1"] <- 'Slope'
  tname[tname == "u2"] <- 'Aspect 1'
  tname[tname == "u3"] <- 'Aspect 2'
  tname[tname == "ringPorous"] <- 'RP xylem'
  tname[tname == "temp"] <- 'Winter temperature'
  tname[tname == "stdage"] <- 'Stand age'
  for(j in length(tname)){
    tname[j] <- paste(toupper(substring(tname[j], 1, 1)), substring(tname[j], 2),sep = "", collapse = " ")
  }
  tname
}

.updateWishartNoPrior <- function(xx,yy,df,beta=NULL,IXX=NULL,WX=NULL,WIX=NULL,
                                  TRYPRIOR=F){
  #df  <- n - Q + S - 1
  S     <- ncol(yy)
  index <- 0
  XX  <- crossprod(xx)
  IXX <- solveRcpp(XX)
  
  D  <- diag(1,nrow(xx)) - xx%*%IXX%*%t(xx)
  SS  <-  t(yy)%*%D%*%yy
  testv <- try(chol(SS),T)
  
  if( inherits(testv,'try-error') ){
    tiny <- 1e-8
    SS[SS < tiny] <- tiny
    message('warning: updateWishartNoPrior')
    SS <- crossprod(yy - xx%*%beta) +  diag(diag(SS)*.001)#*nrow(SS)
    SS <- SS + diag(diag(SS)*.1)
    testv <- try(chol(SS),T)
    index <- 1
  }
  SI <- chol2inv(testv)
  
  z  <- matrix(rnorm(df*S),df,S)%*%chol(SI)
  
  sinv  <- crossprod(z)
  sigma <- solveRcpp(sinv)
  list( sigma = sigma, sinv = sinv, indicator = index )
}

.sqrtRootMatrix <- function(xmat,sigma,DIVIDE=F){
  
  # xmat is n by p
  # sigma is p by p
  
  if(DIVIDE){
    if(length(sigma) == 1)return(xmat/sqrt(sigma))
    return( xmat%*%diag(1/sqrt(diag(sigma))) )
  }
  
  if(length(sigma) == 1)return(xmat*sqrt(sigma))
  xmat%*%diag(sqrt(diag(sigma)) )
}

.yaxisHorizLabs <- function( labels, at=c(1:length(labels)), xshift=.05,
                             col = 'black', pos=NULL){
  
  #add horizontal y axis labels to existing plot
  #pos should be either NULL, 2 (left)
  
  text(par('usr')[3] - xshift*par('usr')[4] - par('usr')[3], y=at,
       labels, xpd=T, pos = pos, col=col)
}

.sampleP <- function(N, avec, bvec, K){
  
  a    <- avec + vapply(1:(N-1), function(k)sum(K == k), 0)
  b    <- bvec + vapply(1:(N-1), function(k)sum(K > k), 0)
  V    <- rbeta((N - 1), a, b)
  p    <- vector("numeric",length=N)
  p[1] <- V[1]
  for(l in 2:(N - 1))p[l] <- prod(1 - V[1:(l - 1)])*V[l]
  p[N] <- prod(1 - V)   
  p
}

### BNP functions for sampling the weights in the PYM process

lt.temp_st_pdf <- function(s, c, sigma, k) {
  exp( - c*( (s+k)^(sigma)  - k^(sigma) ))
}


mult_PY <- function(alpha,sigma, H) {
  Uv<- rgamma(1,alpha/sigma,alpha/sigma)
  # Uv<- rgamma(1,alpha/sigma,1)
  U<- (Uv)^(1/sigma)
  
  x.rlap <- rlaptrans(H, lt.temp_st_pdf, c=alpha/(sigma*H), sigma, k=U)
  #x.rlap <- rlaptrans(H, lt.temp_st_pdf, c=1/H, sigma, k=U)
  pk_vec <- x.rlap /sum(x.rlap)
  return(pk_vec)
}





pdf_lk_mat<- function(l,v, n_k, sigma,H, mat){
  return( (v^l)*exp(mat[n_k,l]))
}

sample_lk_mat<- function(nk_vec,v,sigma,H,M){
  l_post<-c()
  k<- length(nk_vec)
  for (i in 1:k){
    l_vec<- 1:nk_vec[i]
    if (length(l_vec)==1){
      l_post[i]=l_vec
    }
    else{
      p_v<- sapply(l_vec, function(x) pdf_lk_mat(x,v,nk_vec[i],sigma,H,mat=M))
      pv_norm<- p_v/sum(p_v)
      l_post[i]<- sample(1:(nk_vec[i]),size=1, replace=TRUE, prob=pv_norm)
    }
  }
  return(l_post)
}



.sampleP_PYM <- function(N, alpha_val, sigma_val, K, Mat, func){
  n_k<- table(K)
  lh<-  rep(0,N)
  alpha= alpha_val
  sigma=sigma_val
  #sample v
  ptr_logv_comp_mat <- create_xptr("log_v_pdf_comp_mat")
  v_s = ru_rcpp(logf = func,alpha=alpha, sigma=sigma,H=N,k = length(n_k), nk_vec=n_k,Cnk_mat=Mat, n=1,  d=1, init=1)
  #sample lk
  lk <- sample_lk_mat(n_k,v_s$sim_vals[1],sigma,N,Mat)
  lh[c(as.numeric(c(names(n_k))))]= lk
  vh <- rep(0,N)  # initialize
  W_h <- rep(0,N)
  P_h<-  rep(0,N)
  p_vec<- n_k - lk*sigma
  W_h<- rdirichlet(1,c(p_vec, sum(lk)*sigma + alpha))
  ### R
  alpha_post<- alpha + sum(lk)*sigma 
  Uv<- rgamma(1,alpha_post/sigma,alpha_post/sigma)
  U<- (Uv)^(1/sigma)
  x.rlap <- rlaptrans(N, lt.temp_st_pdf, c=alpha_post/(sigma*N), sigma, k=U)
  R_h<- x.rlap /sum(x.rlap)
  P_h[c(as.numeric(c(names(n_k))))]<- W_h[1:length(n_k)] + W_h[length(n_k)]* R_h[1:length(n_k)]
  P_h[-c(as.numeric(c(names(n_k))))] <-  W_h[length(n_k)]* R_h[(length(n_k)+1):N]
  return(P_h)
}









.getPars <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                     alpha.DP, inSamples,...){      
  
  # Y includes all terms but x%*%beta
  
  nn   <- length(inSamples)
  p    <- ncol(x)
  S    <- ncol(Y)
  ntot <- nrow(Y)
  nn   <- length(inSamples)
  
  covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
  rndEff <- RR%*%t(Z[K,])
  
  res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
  
  if(CLUST){   #only until convergence
    avec <- 1/rgamma(r, shape = (2 + r )/2, 
                     rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
    
    D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
    Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
                    Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
    
    pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                         Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                         sigmasqk = sigmaerror)
    K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
    pvec <- .sampleP(N = N, avec = rep(alpha.DP/N,(N-1)),
                     bvec = ((N-1):1)*alpha.DP/N, K = K)  
  }
  
  list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
       sigmaerror = sigmaerror, rndEff = rndEff)
} 


.wWrapperTime <- function(sampleW, y, timeZero, i1, i2, tindex, gindex, uindex,
                          notOther, n, S, REDUCT, RANDOM){
  
  function(w,plo,phi,wpropTime,xl,yp,Lmat,Amat,mub,rndEff, groupRandEff,sdg,muw,
           Umat,Vmat,sinv){
    
    if(RANDOM)muw <- muw + groupRandEff
    
    W <- matrix(.tnorm(n*S,plo,phi,w,wpropTime),n,S)
    W[sampleW == 0] <- y[sampleW == 0]
    
    ii  <- i1
    ni  <- length(ii)
    yp  <- yp*0
    
    for(im in 1:2){
      
      w[sampleW == 0] <- y[sampleW == 0]
      
      if(im == 2)ii <- i2
      
      i00 <- tindex[ii,1]
      i11 <- tindex[ii,2]
      
      ww <- W[i00,]
      ww[ww < 0] <- 0
      
      mugStar <- (ww[,gindex[,'colW']]*xl[i11,gindex[,'rowG']])%*%Lmat
      muaStar <- (ww[,uindex[,1]]*ww[,uindex[,2]] )%*%Amat
      muStar  <- mub[i11,] + mugStar + muaStar + rndEff[i11,]
      
      if(REDUCT){
        pnow <- dnorm(w[i00,notOther],muw[i00,notOther],sdg,log=T) +
          dnorm(w[i11,notOther],muw[i11,notOther],sdg,log=T)
        pnew <- dnorm(ww[,notOther],muw[i00,notOther],sdg,log=T) + 
          dnorm(w[i11,notOther],muStar[,notOther],sdg,log=T)
        za <- which( runif(length(pnow),0,1) < exp(pnew - pnow) )
        if(length(za) > 0){
          w[i00,][za] <- W[i00,][za]
          ww <- w[i00,]
          ww[ww < 0]  <- 0
          muw[i11,][za] <- muStar[za]
          Umat[i11,]    <- ww[,uindex[,1]]*ww[,uindex[,2]]
          Vmat[i11,]    <- ww[,gindex[,'colW']]*xl[i11,gindex[,'rowG']]
        }
      }else{
        #      pnow <- .dMVN(w[i00,notOther],muw[i00,notOther],sinv=sinv,log=T) + 
        #         .dMVN(w[i11,notOther],muw[i11,notOther],sinv=sinv,log=T)
        #      pnew <- .dMVN(ww[,notOther],muw[i00,notOther],sinv=sinv,log=T) + 
        #        .dMVN(w[i11,notOther],muStar[,notOther],sinv=sinv,log=T)
        
        pnow <- .dMVN(w[i00,notOther],muw[i00,notOther],sinv=sinv,log=T) + 
          .dMVN(w[i11,notOther],muw[i11,notOther],sinv=sinv,log=T)
        pnew <- .dMVN(ww[,notOther],muw[i00,notOther],sinv=sinv,log=T) + 
          .dMVN(w[i11,notOther],muStar[,notOther],sinv=sinv,log=T)
        
        za <- which( runif(length(pnow),0,1) < exp(pnew - pnow) )
        if(length(za) > 0){
          w[i00[za],] <- W[i00[za],]
          ww          <- w[i00,]
          ww[ww < 0]  <- 0
          muw[i11[za],] <- muStar[za,]
          Umat[i11,]    <- ww[,uindex[,1]]*ww[,uindex[,2]]
          Vmat[i11,]    <- ww[,gindex[,'colW']]*xl[i11,gindex[,'rowG']]
        }
      }
      W[i00,] <- w[i00,]
    }
    if(REDUCT){
      yp <- matrix(rnorm(n*S,muw,sdg),n,S)
    }else{
      yp[,notOther] <- .rMVN(n,muw[,notOther],sdg[notOther,notOther])
    }
    
    nz <- length(timeZero)
    ww <- w[timeZero,]
    ww[ww < 0] <- 0
    
    mugStar <- (ww[,gindex[,'colW']]*xl[timeZero+1,gindex[,'rowG']])%*%Lmat
    muaStar <- (ww[,uindex[,1]]*ww[,uindex[,2]] )%*%Amat
    muStar  <- mub[timeZero+1,] + mugStar + muaStar + rndEff[timeZero+1,]
    if(RANDOM)muStar <- muStar + groupRandEff[timeZero+1,] ###################
    
    if(REDUCT){
      pnow <- dnorm(w[timeZero+1,notOther],muw[timeZero+1,notOther],sdg,log=T)
      pnew <- dnorm(w[timeZero+1,notOther],muStar[,notOther],sdg,log=T)
      za   <- which( runif(length(pnow),0,1) < exp(pnew - pnow) )
      if(length(za) > 0){
        w[timeZero,][za]   <- W[timeZero,][za]
        ww <- w[timeZero,]
        ww[ww < 0]  <- 0
        muw[timeZero+1,][za] <- muStar[za]
        Umat[timeZero+1,]    <- ww[,uindex[,1]]*ww[,uindex[,2]]
        Vmat[timeZero+1,]    <- ww[,gindex[,'colW']]*xl[timeZero+1,gindex[,'colX']]
      }
    }else{
      pnow <- .dMVN(w[timeZero+1,notOther],muw[timeZero+1,notOther],
                    sinv=sinv,log=T)
      pnew <- .dMVN(w[timeZero+1,notOther],muStar[,notOther],
                    sinv=sinv,log=T)
      
      za   <- which( runif(length(pnow),0,1) < exp(pnew - pnow) )
      if(length(za) > 0){
        w[timeZero[za],]   <- W[timeZero[za],]
        ww <- w[timeZero,]
        ww[ww < 0]  <- 0
        muw[timeZero[za]+1,] <- muStar[za,]
        Umat[timeZero+1,]    <- ww[,uindex[,1]]*ww[,uindex[,2]]
        Vmat[timeZero+1,]  <- ww[,gindex[,'colW']]*xl[timeZero+1,gindex[,'rowG']]
      }
    }
    
    list(Umat = Umat, Vmat = Vmat, w = w, muw = muw, yp = yp)
  }
}

.wWrapper <- function(REDUCT, RANDOM, S, effMat, corCols, notCorCols, typeNames, 
                      typeFull, typeCols, 
                      allTypes, holdoutN, holdoutIndex, censor, 
                      censorCA, censorDA, censorCON, notOther, sampleW, 
                      byRow, byCol,
                      indexW, ploHold, phiHold, sampleWhold, inSamp){
  if(REDUCT){
    
    function(rows=1:nrow(x), x, w, y, bg, sg, alpha, cutg, plo, phi, 
             rndEff, groupRandEff, sigmaerror, wHold){
      
      n    <- nrow(y)
      w0   <- which(sampleW == 1)
      SC   <- ncol(y)
      scol <- c(1:S)
      sigvec <- rep(sigmaerror,S)
      
      if(holdoutN > 0){ # in-sample to predict X out-of-sample
        wHold <- w[drop=F,holdoutIndex,] 
      }
      
      yPredict  <-  w*0
      SN <- length(notCorCols)
      
      if(length(notCorCols) > 0){ ###### covariance scale
        mue    <- x%*%bg
        if(RANDOM)mue <- mue + groupRandEff
        muf    <- mue + rndEff
        
        w[w0]  <- .tnorm(length(w0), plo[w0], phi[w0], muf[w0], sqrt(sigmaerror))
        w[-w0] <- y[-w0]          
        
        if(holdoutN < n){   # in-sample prediction, known RE
          yPredict[,notCorCols] <- rnorm(n*SN,muf[,notCorCols],sqrt(sigmaerror))
        }
        
        if(holdoutN > 0){ # in-sample for holdouts to predict X out-of-sample
          
          # in-sample with RE
          if(holdoutN < n){
            wHold[,notCorCols] <- 
              matrix( .tnorm(holdoutN*SN, as.vector(ploHold[,notCorCols]), 
                             as.vector(phiHold[,notCorCols]),
                             as.vector(muf[drop=F,holdoutIndex,notCorCols] ),
                             sqrt(sigmaerror)),holdoutN,SN) 
          }
          # marginalized RE out-of-sample
          w[holdoutIndex,notOther] <- yPredict[holdoutIndex,notOther] <-
            .rMVN(holdoutN,mue[holdoutIndex,notOther],
                  sg[notOther,notOther]) #out-of-sample RE
        }
      }
      
      if(length(corCols) > 0){   # corr scale
        
        css  <- sg*0
        css[notOther,notOther]  <- .cov2Cor(sg[notOther,notOther])
        
        muo <- x%*%alpha
        if(RANDOM)muo <- muo + groupRandEff
        
        if(holdoutN < n){
          mur <- muo 
          if(length(rndEff) > 1)mur <- mur + .sqrtRootMatrix(rndEff,sg,DIVIDE=T)
          SC  <- length(corCols)
          
          # includes RE on correlation scale
          w[,corCols] <- matrix( .tnorm(n*SC, as.vector(t(plo[,corCols])), 
                                        as.vector(t(phi[,corCols])), 
                                        as.vector(t(mur[,corCols])),1),
                                 n,SC, byrow=T) 
          yPredict[,corCols] <- rnorm(n*SC,mur[,corCols],1)
        }                              
        if(holdoutN > 0){ # out-of-sample
          if(holdoutN < n){
            wHold[,corCols] <- matrix( .tnorm(holdoutN*SC, 
                                              as.vector(ploHold[,corCols]), 
                                              as.vector(phiHold[,corCols]),
                                              as.vector(t(mur[holdoutIndex,corCols])),
                                              1),holdoutN,SC) 
          }
          #      w[holdoutIndex,corCols] <- yPredict[holdoutIndex,corCols]  <- 
          #                                  .rMVN(holdoutN,muo,css[corCols,corCols])
          w[holdoutIndex,corCols] <- yPredict[holdoutIndex,corCols]  <- 
            rmvnormRcpp(holdoutN,rep(0,length(corCols)),
                        css[corCols,corCols]) + muo
        }
      }
      
      if(!is.null(sampleW))w[sampleW == 0] <- y[sampleW == 0]
      if(holdoutN > 0){   # in-sample to sample X out-out-sample
        wHold[sampleWhold == 0] <- y[holdoutIndex,][sampleWhold == 0]  
      }
      
      FCgroups  <- attr(typeNames,'FCgroups')
      CCgroups  <- attr(typeNames,'CCgroups')
      CATgroups <- attr(typeNames,'CATgroups')
      
      for(k in allTypes){
        
        wk <- which(typeCols == k)
        wo <- which(wk %in% notOther)
        nk <- length(wk)
        wu <- which(typeCols[notOther] == k)
        wp <- w[, wk, drop=F]
        yp <- yPredict[, wk, drop=F]
        
        groups <- NULL
        if(typeFull[wk[1]] == 'countComp')  groups <- CCgroups[wk]
        if(typeFull[wk[1]] == 'fracComp')   groups <- FCgroups[wk]
        
        if( typeFull[wk[1]] == 'categorical' ){
          groups <- CATgroups[wk]
          
          if(holdoutN < n){
            tmp <- .gjamWcatLoop2(y, ws = wp, mus = muf, sgs = sigvec, 
                                  notOther = notOther, plo, phi, 
                                  groups = CATgroups, REDUCT=T)
            wp[,wo] <- tmp$w[,wo]
            plo     <- tmp$plo
            phi     <- tmp$phi
          }
          
          if(holdoutN > 0){
            ws <- w[, wk, drop=F]
            ws[holdoutIndex,] <- wHold[, wk, drop=F]
            if(holdoutN < n)wHold[,wo] <- .gjamWcatLoop2(y, ws, mus = muf, 
                                                         sgs = sigvec, 
                                                         notOther = notOther, ploHold, phiHold, 
                                                         groups = CATgroups, REDUCT=T) 
          }
        }
        
        glist <- list(wo = wo, type = typeFull[wk[1]], yy = y[,wk,drop=F],
                      wq = wp, yq = yp, cutg = cutg, censor = censor,
                      censorCA = censorCA, censorDA = censorDA, censorCON = censorCON,
                      eff = effMat[rows,wk,drop=F],groups = groups, k = k, 
                      typeCols = typeCols, notOther = notOther, wk = wk, 
                      sampW = sampleW[,wk])
        
        if(holdoutN < n){
          
          tmp <- .gjamWLoopTypes( glist )  # if PA, yPredict on probit scale
          w[,wk]        <- tmp[[1]]
          yPredict[inSamp,wk] <- tmp[[2]][inSamp,]  # not holdouts
        }
        if(holdoutN > 0){
          
          glist$wq <- wHold[,wk,drop=F]
          glist$yq <- yPredict[holdoutIndex, wk, drop=F]
          glist$yy <- y[holdoutIndex,wk,drop=F]
          glist$eff <- effMat[holdoutIndex, wk, drop=F]
          glist$sampW <- sampleW[,wk]
          
          tmp <- .gjamWLoopTypes( glist )
          
          if(holdoutN < n)wHold[,wk] <- tmp[[1]] #in-sample for x prediction
          yPredict[holdoutIndex,wk] <- tmp[[2]]  #out-of-sample prediction
        }
        yPredict[,wk] <- .censorValues(censor,y,yPredict)[,wk]
      }
      
      if(!is.null(sampleW))w[sampleW[rows,] == 0] <- y[sampleW[rows,] == 0]
      if(holdoutN > 0){
        wHold[sampleWhold == 0] <- y[holdoutIndex,][sampleWhold == 0]
      }
      
      list(w = w, wHold = wHold, yp = yPredict, plo = plo, phi = phi )
    }
    
  } else {
    
    function(rows=1:nrow(x), x, w, y, bg, sg, alpha, cutg, plo, phi, 
             rndEff  = NULL, groupRandEff, sigmaerror = NULL, wHold){
      
      # for holdouts: wHold - w in-sample for sampling x out-out-sample
      #               w[holdoutIndex,] - predict out-of-sample
      
      n     <- nrow(y)
      sampW <- sampleW[rows,notOther]
      w[sampleW[rows,] == 0] <- y[sampleW[rows,] == 0]
      
      if(holdoutN > 0){
        wHold[sampleWhold == 0] <- y[holdoutIndex,][sampleWhold == 0]
      }
      
      yPredict <- w*0
      if(length(notCorCols) > 0){
        muw <- x%*%bg
        if(RANDOM)muw <- muw + groupRandEff
        #      yPredict[,notOther] <- .rMVN(n,muw[,notOther],sg[notOther,notOther])
        yPredict[,notOther] <- rmvnormRcpp(n,rep(0,length(notOther)),
                                           sg[notOther,notOther]) + 
          muw[,notOther]
      }
      
      if( length(corCols) > 0 ){    #expanded w on this scale
        wss  <- w*0
        css  <- .cov2Cor(sg[notOther,notOther])
        muss <- x%*%alpha
        if(RANDOM)muss <- muss + groupRandEff
        ypred <- yPredict
        ypred[,notOther]   <- rmvnormRcpp(n,rep(0,length(notOther)),css) + 
          muss[,notOther]
        yPredict[,corCols] <- ypred[,corCols]
      } 
      
      FCgroups  <- attr(typeNames,'FCgroups')
      CCgroups  <- attr(typeNames,'CCgroups')
      CATgroups <- attr(typeNames,'CATgroups')
      
      for(k in allTypes){
        
        wk <- which(typeCols == k)
        nk <- length(wk)
        wo <- which(wk %in% notOther)
        wu <- which(typeCols[notOther] == k)
        wp <- w[, wk, drop=F]
        yp <- yPredict[, wk, drop=F]
        
        if( typeFull[wk[1]] %in% c('presenceAbsence','ordinal') ) {
          
          wss[,notOther] <- .sqrtRootMatrix(w[,notOther],sg[notOther,notOther],
                                            DIVIDE=T)     
          llist <- list(ws = wss[,notOther], mus = muss[,notOther], 
                        sgs = css, wkk = wu, 
                        lo = plo[,notOther], hi = phi[,notOther],
                        sampW = sampW, indexW = indexW)
          wp[,wo] <- .gjamWLoop( llist )[,wu]
          
          if(holdoutN > 0){
            if(holdoutN < n){
              llist <- list(ws = wss[drop=F,holdoutIndex,notOther], 
                            mus = muss[drop=F,holdoutIndex,notOther], 
                            sgs = css, wkk = wu, 
                            lo = ploHold[drop=F,,notOther], 
                            hi = phiHold[drop=F,,notOther],
                            sampW = sampleWhold[,notOther], indexW=wo)
              wHold[,wo] <- .gjamWLoop( llist )[,wu] 
            }
            wp[holdoutIndex,wo] <- yp[holdoutIndex,wo]
          }
        }
        
        if( !typeFull[wk[1]] %in% c('presenceAbsence','ordinal','categorical') ){
          
          llist <- list(ws = w[,notOther], mus = muw[,notOther],
                        sgs = sg[notOther,notOther], wkk = wu,
                        lo = plo[,notOther], hi = phi[,notOther],sampW = sampW, 
                        indexW = indexW, byCol= byCol, byRow = byRow)
          wp[,wo] <- .gjamWLoop( llist )[,wu]
          
          if(holdoutN > 0){
            if(holdoutN < n){
              llist <- list(ws = w[drop=F,holdoutIndex,notOther], 
                            mus = muw[drop=F,holdoutIndex,notOther], 
                            sgs = sg[notOther,notOther], wkk = wu,
                            lo = ploHold[drop=F,,notOther],
                            hi = phiHold[drop=F,,notOther],
                            sampW = sampleWhold[,notOther], indexW = wo, 
                            byCol = byCol, byRow = byRow)
              wHold[,wo] <- .gjamWLoop( llist )[,wu] 
            }
            wp[holdoutIndex,wo] <- yp[holdoutIndex,wo]
          }
        }
        
        if( typeFull[wk[1]] == 'categorical' ){
          wss[,notOther] <- .sqrtRootMatrix(w[,notOther],sg[notOther,notOther],
                                            DIVIDE=T)
          yy  <- y
          if(holdoutN > 0)yy[holdoutIndex,] <- yp[holdoutIndex,]
          tmp <- .gjamWcatLoop2(yy, ws = wss, mus = muss, sgs = css, 
                                notOther, plo, phi, groups = CATgroups)
          wp      <- tmp$w[,wk]
          plo     <- tmp$plo
          phi     <- tmp$phi
          
          if(holdoutN > 0){
            if(holdoutN < n){
              wHold[,wk] <- .gjamWcatLoop2(yp[drop=F,holdoutIndex,],
                                           wss[drop=F,holdoutIndex,], 
                                           muss[drop=F,holdoutIndex,], sgs = css, 
                                           notOther, ploHold, phiHold, 
                                           groups = CATgroups)$w[,wk]
            }
            wp[holdoutIndex,wo] <- yp[holdoutIndex,wo]
          }
        }
        
        groups <- NULL
        if(typeFull[wk[1]] == 'countComp')  groups <- CCgroups[wk]
        if(typeFull[wk[1]] == 'fracComp')   groups <- FCgroups[wk]
        if(typeFull[wk[1]] == 'categorical')groups <- CATgroups[wk]
        
        glist <- list(wo = wo, type = typeFull[wk[1]], yy = y[,wk,drop=F],
                      wq = wp, yq = yp, cutg = cutg, censor = censor,
                      censorCA = censorCA, censorDA = censorDA, 
                      censorCON = censorCON,
                      eff = effMat[rows,wk,drop=F], groups = groups, k = k, 
                      typeCols = typeCols, notOther = notOther, wk = wk, 
                      sampW = sampleW[,wk])
        tmp <- .gjamWLoopTypes( glist )
        w[,wk]        <- tmp[[1]]
        yPredict[,wk] <- tmp[[2]]
        
        if(holdoutN > 0){
          
          # predict for actual sample size
          ys <- yp[holdoutIndex,,drop=F]
          ys[ys < 0] <- 0
          ys <- rowSums(y[holdoutIndex,wk,drop=F])*ys
          
          glist <- list(wo = wo, type = typeFull[wk[1]], yy = ys, 
                        wq = wp[drop=F,holdoutIndex,], yq = yp[drop=F,holdoutIndex,], 
                        cutg = cutg, censor = censor, censorCA = censorCA, 
                        censorDA = censorDA, censorCON = censorCON, 
                        eff = effMat[drop=F,holdoutIndex,wk], groups = groups, 
                        k = k, typeCols = typeCols, notOther = notOther, wk = wk, 
                        sampW = sampleW[drop=F,holdoutIndex,wk] )
          tmp <- .gjamWLoopTypes( glist )
          w[holdoutIndex,wk] <- tmp[[1]]
          yPredict[holdoutIndex,wk] <- tmp[[2]]    
        }
        
        yPredict[,wk] <- .censorValues(censor, y, yPredict)[,wk]
      }
      
      if(!is.null(sampleW))w[sampleW[rows,] == 0] <- y[sampleW[rows,] == 0]
      
      list(w = w, wHold = wHold, yp = yPredict, plo = plo, phi = phi )
    }
  }
}

.binaryScore <- function(p, x){
  
  #brier and logarithmic score, prediction prob p, event x = 0 or 1
  
  a <- mean((x - p)^2)
  b   <- -mean( x*log(p) + (1 - x)*log(1 - p))
  
  list(brierScore = a, logScore = b)
}


.betaWrapper <- function(REDUCT, TIME, BPRIOR, notOther, IXX, betaLim=50){
  
  # betaLim - outer prior limit for beta
  
  if(REDUCT){
    
    function(X, Y, sig, beta, lo, hi, rows=NULL, pattern=NULL, ixx=F,...){
      
      SS   <- ncol(Y)
      
      w0 <- which(colSums(X) == 0)
      if(length(w0) > 0){
        X <- X[,-w0]
        beta <- beta[-w0,]
        IXX  <- NULL
        rows[rows %in% w0] <- NA
      }
      
      if(is.null(IXX) | !ixx){
        tiny <- 1e-5
        XX   <- crossprod(X)
        diag(XX) <- tiny + diag(XX)          ## ridge here
        IXX   <- try( solve(XX), T )
        if( inherits(IXX,'try-error') ){
          diag(XX) <- diag(XX) + 1.01*diag(XX)
          IXX <- solve(XX)
        }
      }
      
      omega <- sig*IXX
      muB   <- t(omega%*%crossprod((1/sig)*X, Y))
      
      if(!BPRIOR){
        #       B   <- .rMVN( SS, 0, omega) + muB
        
        B   <- rmvnormRcpp( SS, rep(0,nrow(omega)), omega) + muB
        
        ws <- which(abs(B) > betaLim, arr.ind=T)
        
        if(length(ws) > 0){
          ws <- unique(ws[,1])
          bs <- B[drop=F,ws,]
          B[ws,] <- .tnormMVNmatrix(avec = bs, muvec = muB[drop=F,ws,], 
                                    smat = omega, lo = bs*0 - betaLim, 
                                    hi = bs*0 + betaLim)
        }
        return(t(B))
      }
      
      if(!TIME){
        
        tmp <- .tnormMVNmatrix(avec = t(beta), muvec = muB, 
                               smat = omega, lo = t(lo), 
                               hi = t(hi))
        return( t(tmp) )
      }
      
      B  <- t(beta)
      QX <- ncol(X)
      
      for(k in 1:nrow(rows)){
        krow <- rows[k,]
        krow <- krow[is.finite(krow)]
        notk <- c(1:QX)[-krow]
        if(length(notk) == 1){
          M1 <- omega[krow,notk, drop=F]/omega[notk,notk]
        }else{
          OI <- try( solveRcpp(omega[notk,notk]), T)
          if( inherits(OI,'try-error') ){
            OI <- diag(1/diag(omega[notk,notk]))
          }
          M1 <- omega[krow,notk, drop=F]%*%OI
        }
        pk  <- pattern[k,]
        pk  <- pk[is.finite(pk)]
        muk <- muB[pk, krow, drop=F] - muB[pk,notk]%*%t(M1)
        Mk  <- omega[krow,krow] - M1%*%omega[notk,krow]
        
        if(length(Mk) == 1){
          B[pk,krow] <- .tnorm(length(pk),lo[krow,pk],hi[krow,pk],muk,sqrt(Mk))
        } else {
          ll <- t(lo)[pk,krow,drop=F]
          hh <- t(hi)[pk,krow,drop=F]
          test <- try( .tnormMVNmatrix( avec=muk, muvec=muk, smat=Mk,
                                        lo=ll, hi=hh), T)
          if( inherits(test,'try-error') ){
            mm <- diag(Mk)
            mm[mm < tiny] <- tiny
            test <- .tnorm(length(ll),ll,hh,muk,sqrt(mm))
          }
          B[pk,krow] <- test
        }
      }
      return( t(B) )
    }
    
  }else{
    
    if(!BPRIOR){
      
      function(X, Y, sig,...){
        
        if(is.null(IXX)){
          XX    <- crossprod(X)
          IXX <- chol2inv(chol( XX ) )
        }
        WX  <- crossprod(X,Y)
        WIX <- IXX%*%WX
        bg  <- matrix( .rMVN(1,as.vector(WIX),
                             kronecker(sig,IXX)),nrow(IXX),ncol(WIX) )
        return(bg)
      }
      
    } else{
      
      function(X, Y, sig, beta, lo, hi, ...){
        
        if(is.null(IXX)){
          XX    <- crossprod(X)
          IXX <- chol2inv(chol( XX ) )
        }
        WX  <- crossprod(X,Y)
        WIX <- IXX%*%WX
        smat <- kronecker(sig,IXX)
        tmp <- .tnormMVNmatrix(avec = matrix(beta,1), muvec = matrix(WIX,1), 
                               smat = smat, lo = matrix(lo,1), 
                               hi = matrix(hi,1))
        tmp <- matrix(tmp,nrow(beta),ncol(beta))
        tmp[!is.finite(tmp)] <- beta[!is.finite(tmp)]
        return(tmp)
      }
    }
  }
}

.paramWrapper <- function(REDUCT, inSamples,SS){   
  
  if(REDUCT){    
    
    function(CLUST, x,beta,Y,otherpar){
      
      N  <- otherpar$N
      r  <- otherpar$r
      D  <- otherpar$D
      Z  <- otherpar$Z
      sigmaerror <- otherpar$sigmaerror
      K          <- otherpar$K
      pvec       <- otherpar$pvec
      alpha.DP   <- otherpar$alpha.DP
      tmp        <- .getPars(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta), 
                             D = D, Z = Z, sigmaerror = sigmaerror,
                             K = K, pvec = pvec, alpha.DP = alpha.DP,
                             inSamples = inSamples, SELECT = F)
      
      sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
                                   K = tmp$K, REDUCT=T))
      
      otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
                       sigmaerror = tmp$sigmaerror,
                       pvec = tmp$pvec, K = tmp$K, alpha.DP = alpha.DP)
      
      return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
    }
    
  } else {
    
    function(CLUST, x, beta,Y,otherpar){
      
      sigmaDf  <- otherpar$sigmaDf
      XX  <- crossprod(x[inSamples,])
      IXX <- solveRcpp(XX)
      WX  <- crossprod(x[inSamples,], Y[inSamples,])
      WIX <- IXX%*%WX
      
      sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
                                   beta = beta, IXX = IXX, WX = WX, WIX = WIX,
                                   TRYPRIOR = T)$sigma
      otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
      
      return(list(sg = sg, otherpar = otherpar))
    }
  }
}

.rwish <- function(df,SS){
  z  <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%chol(SS)
  crossprod(z)
}

.riwish <- function(df,S){
  solveRcpp(.rwish(df,solveRcpp(S)))
}

.expandSigmaChains <- function(snames, sgibbs, otherpar, 
                               simIndex = sample(nrow(sgibbs),50,replace=T), 
                               sigErrGibbs, kgibbs=NULL, 
                               REDUCT=F, CHAINSONLY=F){
  tiny <- 1e-8
  
  S <- otherpar$S
  K <- otherpar$K
  N <- otherpar$N
  r <- otherpar$r
  if(length(simIndex) > 1000)simIndex <- sample(simIndex,1000)
  ns     <- length(simIndex)
  xnames <- otherpar$xnames
  
  if(CHAINSONLY & !REDUCT){  #only return expanded sgibbs
    
    imat   <- matrix(1:(S*S),S,S)
    jmat   <- matrix(1:(S*S),S,S,byrow=T)
    tmp    <- matrix(NA,nrow(sgibbs),S*S)
    sindex <- imat[lower.tri(imat,diag=T)]
    tmp[,sindex] <- sgibbs
    sindex <- jmat[lower.tri(imat,diag=T)]
    tmp[,sindex] <- sgibbs
    
    sMu <- matrix( colMeans(tmp),S,S)
    sSe <- matrix( apply(tmp,2,sd),S,S)
    
    chainList <- list(cchain = NULL, schain = tmp, kchain = NULL)
    
    return( list(chainList = chainList, rMu = NULL, rSe = NULL, 
                 sMu = sMu, sSe = sSe) )
  }
  
  # summarize chains
  
  other    <- grep('other',snames)
  notOther <- c(1:S)
  if(length(other) > 0)notOther <- notOther[-other]
  
  Kindex <- which(lower.tri( diag(S),diag=T ) )
  kchain <- NULL
  
  schain <- cchain <- matrix(0,ns,length(Kindex))
  if(REDUCT)kchain <- matrix(0,ns,ncol(kgibbs))
  colnames(schain) <- colnames(cchain) <- .multivarChainNames(snames,snames)[Kindex]
  
  snames <- otherpar$snames
  s1 <- diag(S)*0
  s2 <- r1 <- r2 <- s1
  
  message('expanding covariance chains')
  
  pbar <- txtProgressBar(min=1,max=ns,style=1)
  
  sinvPlus <-  sinvMinus <- matrix(0,S,S)   # different from zero
  
  k <- 1
  
  for(j in simIndex){
    if(REDUCT){
      Z  <- matrix(sgibbs[j,],N,r)
      ss <- .expandSigma(sigErrGibbs[j], S, Z = Z, kgibbs[j,], REDUCT = REDUCT)
      si <- invWbyRcpp(sigErrGibbs[j], Z[kgibbs[j,],])
      cc <- .cov2Cor(ss)
      dc <- diag(sqrt(diag(ss)))
      ci <- dc%*%si%*%dc
    } else {
      ss <- .expandSigma(sgibbs[j,], S = S, REDUCT = REDUCT)
      si <- ci <- diag(1,S)
      si[notOther,notOther] <- solveRcpp(ss[notOther,notOther])
      cc <- .cov2Cor(ss)
      ci[notOther,notOther] <- solveRcpp(cc[notOther,notOther])
    }
    
    s1 <- s1 + ss
    s2 <- s2 + ss^2
    r1 <- r1 + cc
    r2 <- r2 + cc^2
    
    if(!CHAINSONLY){
      schain[k,]    <- ss[Kindex]
      cchain[k,]    <- cc[Kindex]
      if(REDUCT)kchain[k,] <- kgibbs[j,]
    }
    
    sinvPlus[si > 0]  <- sinvPlus[si > 0] + 1
    sinvMinus[si < 0] <- sinvMinus[si < 0] + 1
    
    setTxtProgressBar(pbar,k)
    k <- k + 1
  }
  diag(sinvPlus) <- diag(sinvMinus) <- 0
  sigInvPos <- which(sinvPlus > .95*length(simIndex),arr.ind=T)
  sigInvNeg <- which(sinvMinus > .95*length(simIndex),arr.ind=T)
  
  ssi <- sort( unique(c( sigInvPos[,1], sigInvNeg[,1]) ) )
  
  sMu  <- s1/ns
  vv   <- s2/ns - sMu^2
  vv[vv < tiny] <- tiny
  sSe  <- sqrt( vv )
  rMu  <- r1/ns
  vv   <- r2/ns - rMu^2
  vv[vv < tiny] <- tiny
  rSe  <- sqrt( vv )
  
  rownames(sMu)    <- colnames(sMu) <- snames
  rownames(sSe)    <- colnames(rSe) <- snames
  colnames(cchain) <- colnames(schain)
  
  chainList <- list(cchain = cchain, schain = schain, kchain = kchain)
  
  list(chainList = chainList, rMu = rMu, rSe = rSe, 
       sMu = sMu, sSe = sSe)
}

.expandSigma <- function(sigma, S, Z = NULL, K = NULL, REDUCT = F){
  
  if(REDUCT) return( sigma*diag(S) + tcrossprod(Z[K,]) )
  
  ss <- diag(S)
  ss[lower.tri(ss,diag=T)] <- sigma
  ss[upper.tri(ss)] <- t(ss)[upper.tri(ss)]
  ss
}

.ordTraitsFromWts <- function(yWt,ordTraits){
  
  # yWt - n by S species weights
  # ordTraits - S by p ordinal traits
  # returns n by p modal ordinal values
  
  if(!is.matrix(ordTraits))ordTraits <- matrix(ordTraits)
  
  n <- nrow(yWt)
  s <- ncol(yWt)
  
  ii <- rep(c(1:n),s)
  omat <- matrix(NA,n,ncol(ordTraits))
  
  for(j in 1:ncol(ordTraits)){
    
    PLUS <- F
    
    oj  <- ordTraits[,j]
    if(min(oj) < 0)stop('ordinal scores cannot be < 0')
    if(min(oj) == 0){
      PLUS <- T
      oj   <- oj + 1
    }
    
    rj  <- range(oj, na.rm=T)
    mm  <- matrix(0, n, rj[2] )
    jj  <- as.vector( matrix(oj, n, s, byrow=T) )
    tmp <- .byGJAM(as.vector(yWt),ii,jj,mm,mm,fun='sum')
    w0  <- which( apply(tmp,1,sum) == 0)
    
    m1  <- apply(tmp,1,which.max)
    m1  <- (rj[1]:rj[2])[m1]
    if(PLUS)m1 <- m1 - 1
    omat[,j] <- m1
    if(length(w0) > 0)omat[w0,j] <- 0
  }
  colnames(omat) <- colnames(ordTraits)
  omat
}

.incidence2Grid <- function(specs, lonLat, nx = NULL, ny = NULL, dx = NULL, 
                            dy = NULL, predGrid = NULL, effortOnly=TRUE){
  
  # must have either ngrid X 2 prediction grid, or 
  #   numbers of points nx, ny, or
  #   densities of points dx, dy
  
  ngrid <- length(predGrid)
  mapx  <- range(lonLat[,1])
  mapy  <- range(lonLat[,2])
  
  specs  <- as.character(specs)
  
  ynames <- sort(unique(specs))
  nspec  <- length(ynames)
  jj     <- match(specs,ynames)
  
  if(ngrid == 0){
    if(!is.null(dx)){
      xseq <- seq(mapx[1], mapx[2], by = dx)
      yseq <- seq(mapy[1], mapy[2], by = dy)
    } else {
      xseq <- seq(mapx[1], mapx[2], length = nx)
      yseq <- seq(mapy[1], mapy[2], length = ny)
    }
    predGrid <- as.matrix( expand.grid(lon = xseq, lat = yseq) )
    ngrid    <- nrow(predGrid)
  }
  
  ii <- RANN::nn2(predGrid, lonLat, k = 1  )$nn.idx
  mm <- matrix(0, ngrid, nspec )
  
  gridBySpec <- .byGJAM(ii*0 + 1, ii, jj, mm, mm, fun='sum')
  colnames(gridBySpec) <- ynames
  effort <- rowSums(gridBySpec)
  
  if(effortOnly){
    wk <- which(effort > 0)
    effort     <- effort[wk]
    gridBySpec <- gridBySpec[wk,]
    predGrid   <- predGrid[wk,]
  }
  list(gridBySpec = gridBySpec, predGrid = predGrid)
}

.spec2Trait <- function(pbys, sbyt, tTypes){
  
  # plotBySpec  - n by S numeric matrix
  # specByTrait - S by M data.frame
  # traitTypes  - data types for traits
  # FC can be factors that will be categorical
  
  n <- nrow(pbys)
  S <- ncol(pbys)
  M <- ncol(sbyt)
  
  ttt <- numeric(0)
  
  y2t  <- match(colnames(pbys),rownames(sbyt))  
  y2tf <- which(is.finite(y2t))
  t2y  <- match(rownames(sbyt),colnames(pbys))
  t2yf <- which(is.finite(t2y))
  
  if(is.data.frame(pbys))pbys <- as.matrix(pbys)
  
  ywt <- sweep(pbys,1,rowSums(pbys,na.rm=T),'/')
  ywt[is.na(ywt)] <- 0
  
  newTypes <- character(0)
  tmat     <- ttt <- numeric(0)
  
  ###################### neither ordinal nor factors (FC)
  
  wf   <- which(!tTypes %in% c('OC','CAT')) 
  
  if(length(wf) > 0){
    newTypes <- tTypes[wf]
    ttt <- sbyt[y2t,wf, drop=F]
    tmat <- ywt%*%as.matrix(sbyt[y2t,wf, drop=F])
  }
  
  ###################### ordinal classes
  
  ordNames <- which(tTypes == 'OC')
  
  if(length(ordNames) > 0){
    ordTraits <- as.matrix( round(sbyt[y2t[y2tf],ordNames],0) )
    ordCols   <- .ordTraitsFromWts(ywt,ordTraits)
    if(is.null(colnames(ordCols)))colnames(ordCols) <- colnames(ordTraits) <- 
      colnames(sbyt)[ordNames]
    ttt <- cbind(ttt, ordTraits )
    tmat <- cbind(tmat,ordCols)
    newTypes <- c(newTypes,tTypes[ordNames])
  }
  
  ##################### CAT to FC
  
  censor <- NULL
  mcol   <- ncol(tmat)
  if(is.null(mcol))mcol <- 0
  xx     <- numeric(0)
  FCgroups <- rep(0,mcol)   
  
  wf <- numeric(0)
  for(j in 1:ncol(sbyt))if(is.factor(sbyt[,j]))wf <- c(wf,j)
  
  wf <- union(wf,which(tTypes %in% 'CAT'))
  
  if(length(wf) > 0){
    
    xx <- sbyt[,wf,drop=F]
    xc <- numeric(0)
    kg <- 0
    
    for(kk in 1:length(wf)){
      
      xkk  <- xx[[kk]]            #rare type is reference
      xtab <- table(xkk)
      if(length(xtab) == 1){
        stop( paste('CAT trait _',names(xx)[kk],
                    '_ has only 1 level, need at least 2',sep='') )
      }
      
      xtab <- xtab[order(xtab)]
      xkk  <- relevel(xkk,ref=names(xtab)[1])
      cont <- contrasts(xkk,contrasts = F)
      xk   <- cont[xkk,,drop=F]
      tmp  <- ywt[,t2y]%*%xk[t2y,]
      
      if(ncol(tmp) == 2){
        mc    <- mcol + 1
        ktype <- 'CA'
        tmp <- tmp[,1,drop=F]
        gk  <- 0
        tc <- gjamCensorY( values = c(0,1), 
                           intervals = cbind( c(-Inf,0), c(1,Inf) ),
                           y = tmp)
        ttt <- cbind(ttt,xk[,1,drop=F])
        
        if(is.null(censor)){
          censor <- c(censor, tc$censor)
          censor$CA$columns <- mc
        } else {
          censor$CA$columns <- c(censor$CA$columns,mc)
        }
        
      } else {
        
        mc    <- ncol(tmp)
        cname <- colnames(tmp)
        cname[1] <- 'other'
        cname <- paste(colnames(xx)[kk],cname,sep='')
        colnames(tmp) <- colnames(xk) <- cname
        
        ttt   <- cbind(ttt,xk)
        ktype <- rep('FC',ncol(tmp))
        
        kg    <- kg + 1
        gk    <- rep(kg,mc)
      }
      
      mcol <- mcol + ncol(tmp)
      
      FCgroups <- c(FCgroups,gk)
      xc   <- cbind(xc,tmp)
      newTypes <- c(newTypes,ktype)
      
    }
    tmat <- cbind(tmat,xc)
  }
  
  colnames(tmat) <- colnames(ttt)
  
  attr(newTypes,'FCgroups') <- FCgroups
  
  list(plotByCWM = tmat, traitTypes = newTypes, censor = censor,
       specByTrait = ttt)
}

.boxplotQuant <- function( xx, ..., boxfill=NULL ){
  
  tmp <- boxplot( xx, ..., plot=F)
  ss  <- apply( xx, 2, quantile, pnorm(c(-1.96,-1,0,1,1.96)) ) 
  tmp$stats <- ss
  
  pars <- list(...)
  if( 'col' %in% names(pars) )boxfill <- pars$col
  
  bxp( tmp, ..., boxfill = boxfill )
  
  tmp
}

.gjamOrd <- function( output, SPECLABS, col, cex, PLOT, method ){
  
  # method can be 'PCA' or 'NMDS'
  
  ematrix   <- output$parameters$ematrix
  ematAlpha <- output$modelList$ematAlpha
  
  whConZero <- output$fit$whConZero
  whichZero <- output$fit$whichZero
  
  y <- output$inputs$y
  S <- SO <- ncol(y)
  snames  <- colnames(y)
  
  if(is.null(col))col <- rep('black',S)
  
  other <- output$inputs$other
  notOther <- output$inputs$notOther
  
  SO <- length(notOther)
  
  plab <- c('Axis I', 'Axis II', 'Axis III')
  
  if (method == 'NMDS') {
    tmp    <- isoMDS(.cov2Dist(ematrix[notOther,notOther]), k = 3)
    eVecs  <- tmp$points
    colnames(eVecs) <- paste('NMDS',c(1:3),sep = '_')
    eValues <- lambda <- cl <- NULL
  } else {
    tmp <- eigen(ematrix[notOther,notOther])    # PCA
    eVecs   <- tmp$vectors
    eValues <- tmp$values
    lambda  <- eValues/sum(eValues)
    cl      <- cumsum(lambda)
    clab    <- paste(' (',round(100*lambda,0),'%)',sep='')
    plab    <- paste(plab, clab, sep='')
  }
  rownames(eVecs) <- snames[notOther]
  
  if(!PLOT) return( list(eVecs = eVecs, eValues = eValues) )
  
  cbord <- .getColor(col[notOther],.6)
  
  par(mfcol=c(2,2), bty='n', cex = cex, mar=c(4,4,1,1))
  
  plot(eVecs[,1],eVecs[,2],cex=1,col=cbord, bg = cbord, pch=16,
       xlab=plab[1], ylab = plab[2]) 
  abline(h=0,col=.getColor('black',.3),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.3),lwd=2,lty=2)
  
  if(length(SPECLABS) > 0){
    mmm <- match(SPECLABS,rownames(eVecs))
    text(eVecs[mmm,2],eVecs[mmm,3],SPECLABS,col=cbord[notOther][mmm])
  }
  
  plot(eVecs[,1],eVecs[,3],cex=1,col=cbord, bg = cbord, pch=16,
       xlab=plab[1], ylab = plab[3]) 
  abline(h=0,col=.getColor('black',.3),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.3),lwd=2,lty=2)
  
  if(length(SPECLABS) > 0){
    mmm <- match(SPECLABS,rownames(eVecs))
    text(eVecs[mmm,2],eVecs[mmm,3],SPECLABS,col=cbord[notOther][mmm])
  }
  
  plot(eVecs[,2],eVecs[,3],cex=1,col=cbord, bg = cbord, pch=16,
       xlab=plab[2], ylab = plab[3])
  abline(h=0,col=.getColor('black',.3),lwd=2,lty=2)
  abline(v=0,col=.getColor('black',.3),lwd=2,lty=2)
  
  if(length(SPECLABS) > 0){
    mmm <- match(SPECLABS,rownames(eVecs))
    text(eVecs[mmm,2],eVecs[mmm,3],SPECLABS,col=cbord[notOther][mmm])
  }
  
  if(method == 'PCA'){
    plot(cl,type='s',xlab='Rank',ylab='Proportion of variance',xlim=c(.9,S),
         ylim=c(0,1),log='x',lwd=2)
    lines(c(.9,1),c(0,cl[1]),lwd=2,type='s')
    for(j in 1:length(lambda))lines(c(j,j),c(0,cl[j]),col='grey')
    lines(cl,lwd=2,type='s')
    abline(h=1,lwd=2,col=.getColor('grey',.5),lty=2)
  }
  
  list(eVecs = eVecs, eValues = eValues)
} 


columnSplit <- function(vec, sep='_', ASFACTOR = F, ASNUMERIC=F,
                        LASTONLY=F){
  
  vec <- as.character(vec)
  nc  <- length( strsplit(vec[1], sep)[[1]] )
  mat <- matrix( unlist( strsplit(vec, sep) ), ncol=nc, byrow=T )
  if(LASTONLY & ncol(mat) > 2){
    rnn <- mat[,1]
    for(k in 2:(ncol(mat)-1)){
      rnn <- columnPaste(rnn,mat[,k])
    }
    mat <- cbind(rnn,mat[,ncol(mat)])
  }
  if(ASNUMERIC){
    mat <- matrix( as.numeric(mat), ncol=nc )
  }
  if(ASFACTOR){
    mat <- data.frame(mat)
  }
  mat
}

columnPaste <- function(c1, c2, sep='-'){
  
  FACT <- T
  if(!is.factor(c1))FACT <- F
  c1    <- as.character(c1)
  c2    <- as.character(c2)
  #  c1    <- .fixNames(c1)
  #  c2    <- .fixNames(c2)
  c12   <- apply( cbind(c1, c2) , 1, paste0, collapse=sep)
  c12   <- .replaceString(c12, ' ', '')
  if(FACT) c12 <- as.factor(c12)
  c12
}


# 
# 
# 
# .getPars_1 <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
#                        alpha.DP, inSamples,shape,rate,...){      
#   
#   # Y includes all terms but x%*%beta
#   
#   nn   <- length(inSamples)
#   p    <- ncol(x)
#   S    <- ncol(Y)
#   ntot <- nrow(Y)
#   nn   <- length(inSamples)
#   
#   covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
#   z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
#   RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
#   if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
#   rndEff <- RR%*%t(Z[K,])
#   
#   res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
#   sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
#   
#   if(CLUST){   #only until convergence
#     avec <- 1/rgamma(r, shape = (2 + r )/2, 
#                      rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
#     
#     D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
#     Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
#                     Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
#     
#     pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
#                          Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
#                          sigmasqk = sigmaerror)
#     K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
#     
#     #pvec <- .sampleP(N = N, avec = rep(alpha.DP/N,(N-1)),
#     #                 bvec = ((N-1):1)*alpha.DP/N, K = K)  
#     pvec <- .sampleP(N=N, avec=rep(1,(N-1)),
#                      bvec=rep(alpha.DP,(N-1)), K=K)
#     
#     alpha.DP<-rgamma(1, shape=N+shape-1, rate = rate-log(pvec[N]))
#     print(c(shape.rate))
#   }
#   
#   list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
#        sigmaerror = sigmaerror, rndEff = rndEff,alpha.DP=alpha.DP,shape=shape,rate=rate)
# } 

# .paramWrapper_1 <- function(REDUCT, inSamples,SS){   
#   
#   if(REDUCT){    
#     
#     function(CLUST, x,beta,Y,otherpar){
#       
#       N  <- otherpar$N
#       r  <- otherpar$r
#       D  <- otherpar$D
#       Z  <- otherpar$Z
#       sigmaerror <- otherpar$sigmaerror
#       K          <- otherpar$K
#       pvec       <- otherpar$pvec
#       alpha.DP   <- otherpar$alpha.DP
#       rate       <- otherpar$rate
#       shape       <- otherpar$shape
#       tmp        <- .getPars_1(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta), 
#                                D = D, Z = Z, sigmaerror = sigmaerror,
#                                K = K, pvec = pvec, alpha.DP = alpha.DP,
#                                inSamples = inSamples, SELECT = F,shape=shape,rate=rate)
#       
#       sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
#                                    K = tmp$K, REDUCT=T))
#       
#       otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
#                        sigmaerror = tmp$sigmaerror,
#                        pvec = tmp$pvec, K = tmp$K, alpha.DP = tmp$alpha.DP,shape=tmp$shape,rate=tmp$rate)
#       
#       return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
#     }
#     
#   } else {
#     
#     function(CLUST, x, beta,Y,otherpar){
#       
#       sigmaDf  <- otherpar$sigmaDf
#       XX  <- crossprod(x[inSamples,])
#       IXX <- solveRcpp(XX)
#       WX  <- crossprod(x[inSamples,], Y[inSamples,])
#       WIX <- IXX%*%WX
#       
#       sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
#                                    beta = beta, IXX = IXX, WX = WX, WIX = WIX,
#                                    TRYPRIOR = T)$sigma
#       otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
#       
#       return(list(sg = sg, otherpar = otherpar))
#     }
#   }
# }
# 
# 

.getPars_1 <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                       alpha.DP, inSamples,rate,shape,alpha.DP_vec,...){      
  
  # Y includes all terms but x%*%beta
  
  nn   <- length(inSamples)
  p    <- ncol(x)
  S    <- ncol(Y)
  ntot <- nrow(Y)
  nn   <- length(inSamples)
  
  covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
  rndEff <- RR%*%t(Z[K,])
  
  res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
  
  if(CLUST){   #only until convergence
    avec <- 1/rgamma(r, shape = (2 + r )/2, 
                     rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
    
    D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
    Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
                    Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
    
    pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                         Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                         sigmasqk = sigmaerror)
    K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
    
    pvec <- .sampleP(N = N, avec = rep(alpha.DP/N,(N-1)),
                     bvec = ((N-1):1)*alpha.DP/N, K = K)  
    #pvec <- .sampleP(N=N, avec=rep(1,(N-1)),
    #                bvec=rep(alpha.DP,(N-1)), K=K)
    
    alpha.DP<-metrop_DP(theta=alpha.DP,pvec=pvec,lik.fun=lik.alpha.DP.fun,N=N,rate=rate,shape=shape,alpha.DP_vec=alpha.DP_vec)
    
  }
  
  list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
       sigmaerror = sigmaerror, rndEff = rndEff,alpha.DP=alpha.DP,rate,shape,alpha.DP_vec=c(alpha.DP_vec,alpha.DP))
} 


.paramWrapper_1 <- function(REDUCT, inSamples,SS){   
  
  if(REDUCT){    
    
    function(CLUST, x,beta,Y,otherpar){
      
      N  <- otherpar$N
      r  <- otherpar$r
      D  <- otherpar$D
      Z  <- otherpar$Z
      sigmaerror <- otherpar$sigmaerror
      K          <- otherpar$K
      pvec       <- otherpar$pvec
      alpha.DP   <- otherpar$alpha.DP
      rate       <- otherpar$rate
      shape       <- otherpar$shape
      alpha.DP_vec <-otherpar$alpha.DP_vec
      
      tmp        <- .getPars_1(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta), 
                               D = D, Z = Z, sigmaerror = sigmaerror,
                               K = K, pvec = pvec, alpha.DP = alpha.DP, shape=shape,rate=rate, alpha.DP_vec=alpha.DP_vec,
                               inSamples = inSamples, SELECT = F)
      tmp
      sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
                                   K = tmp$K, REDUCT=T))
      
      otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
                       sigmaerror = tmp$sigmaerror,
                       pvec = tmp$pvec, K = tmp$K, alpha.DP = tmp$alpha.DP, shape= shape,rate= rate,alpha.DP_vec=tmp$alpha.DP_vec)
      
      return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
    }
    
  } else {
    
    function(CLUST, x, beta,Y,otherpar){
      
      sigmaDf  <- otherpar$sigmaDf
      XX  <- crossprod(x[inSamples,])
      IXX <- solveRcpp(XX)
      WX  <- crossprod(x[inSamples,], Y[inSamples,])
      WIX <- IXX%*%WX
      
      sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
                                   beta = beta, IXX = IXX, WX = WX, WIX = WIX,
                                   TRYPRIOR = T)$sigma
      otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
      
      return(list(sg = sg, otherpar = otherpar))
    }
  }
}


.paramWrapper_2 <- function(REDUCT, inSamples,SS){   
  
  if(REDUCT){    
    
    function(CLUST, x,beta,Y,otherpar){
      
      N  <- otherpar$N
      r  <- otherpar$r
      D  <- otherpar$D
      Z  <- otherpar$Z
      sigmaerror <- otherpar$sigmaerror
      K          <- otherpar$K
      pvec       <- otherpar$pvec
      alpha_py   <- otherpar$alpha.PY
      sigma_py   <- otherpar$discount.PY
      matrix_cnk <-  otherpar$matrixCnk
      ptr_logv_comp_mat <- otherpar$fun_pointer
      tmp        <- .getPars_2(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta), 
                               D = D, Z = Z, sigmaerror = sigmaerror,
                               K = K, pvec = pvec, alpha.py = alpha_py,sigma.py=sigma_py,
                               inSamples = inSamples, ,matrixCnk = matrix_cnk, fun_pointer =ptr_logv_comp_mat,SELECT = F)
      
      sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
                                   K = tmp$K, REDUCT=T))
      
      otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
                       sigmaerror = tmp$sigmaerror,
                       pvec = tmp$pvec, K = tmp$K, alpha.PY = alpha_py,discount.PY=sigma_py,matrixCnk = matrix_cnk,fun_pointer = ptr_logv_comp_mat  )
      
      return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
    }
    
  } else {
    
    function(CLUST, x, beta,Y,otherpar){
      
      sigmaDf  <- otherpar$sigmaDf
      XX  <- crossprod(x[inSamples,])
      IXX <- solveRcpp(XX)
      WX  <- crossprod(x[inSamples,], Y[inSamples,])
      WIX <- IXX%*%WX
      
      sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
                                   beta = beta, IXX = IXX, WX = WX, WIX = WIX,
                                   TRYPRIOR = T)$sigma
      otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
      
      return(list(sg = sg, otherpar = otherpar))
    }
  }
}


.getPars_2 <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                       alpha.py, sigma.py, inSamples,matrixCnk,fun_pointer,...){      
  
  # Y includes all terms but x%*%beta
  
  nn   <- length(inSamples)
  p    <- ncol(x)
  S    <- ncol(Y)
  ntot <- nrow(Y)
  nn   <- length(inSamples)
  
  covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
  rndEff <- RR%*%t(Z[K,])
  
  res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
  
  if(CLUST){   #only until convergence
    avec <- 1/rgamma(r, shape = (2 + r )/2, 
                     rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
    
    D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
    Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
                    Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
    
    pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                         Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                         sigmasqk = sigmaerror)
    K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
   # pvec <- .sampleP(N = N, avec = rep(1-sigma.py,(N-1)),
  #                   bvec = ((1:(N-1))*sigma.py + alpha.DP), K = K)  
   
    pvec <- .sampleP_PYM(N = N, alpha_val = alpha.py, sigma_val = sigma.py, K = K, Mat = matrixCnk, func =fun_pointer )  
    
     #alphaDP_g<- rgamma(1+N , 1/2 - log(pvec[N]))
    
  }
  
  list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
       sigmaerror = sigmaerror, rndEff = rndEff)
} 




.paramWrapper_2 <- function(REDUCT, inSamples,SS){   
  
  if(REDUCT){    
    
    function(CLUST, x,beta,Y,otherpar){
      
      N  <- otherpar$N
      r  <- otherpar$r
      D  <- otherpar$D
      Z  <- otherpar$Z
      sigmaerror <- otherpar$sigmaerror
      K          <- otherpar$K
      pvec       <- otherpar$pvec
      alpha_py   <- otherpar$alpha.PY
      sigma_py   <- otherpar$discount.PY
      matrix_cnk <-  otherpar$matrixCnk
      ptr_logv_comp_mat <- otherpar$fun_pointer
      tmp        <- .getPars_2(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta), 
                               D = D, Z = Z, sigmaerror = sigmaerror,
                               K = K, pvec = pvec, alpha.py = alpha_py,sigma.py=sigma_py,
                               inSamples = inSamples, ,matrixCnk = matrix_cnk, fun_pointer =ptr_logv_comp_mat,SELECT = F)
      
      sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
                                   K = tmp$K, REDUCT=T))
      
      otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
                       sigmaerror = tmp$sigmaerror,
                       pvec = tmp$pvec, K = tmp$K, alpha.PY = alpha_py,discount.PY=sigma_py,matrixCnk = matrix_cnk,fun_pointer = ptr_logv_comp_mat  )
      
      return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
    }
    
  } else {
    
    function(CLUST, x, beta,Y,otherpar){
      
      sigmaDf  <- otherpar$sigmaDf
      XX  <- crossprod(x[inSamples,])
      IXX <- solveRcpp(XX)
      WX  <- crossprod(x[inSamples,], Y[inSamples,])
      WIX <- IXX%*%WX
      
      sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
                                   beta = beta, IXX = IXX, WX = WX, WIX = WIX,
                                   TRYPRIOR = T)$sigma
      otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
      
      return(list(sg = sg, otherpar = otherpar))
    }
  }
}


.getPars_2 <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                       alpha.py, sigma.py, inSamples,matrixCnk,fun_pointer,...){      
  
  # Y includes all terms but x%*%beta
  
  nn   <- length(inSamples)
  p    <- ncol(x)
  S    <- ncol(Y)
  ntot <- nrow(Y)
  nn   <- length(inSamples)
  
  covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
  rndEff <- RR%*%t(Z[K,])
  
  res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
  
  if(CLUST){   #only until convergence
    avec <- 1/rgamma(r, shape = (2 + r )/2, 
                     rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
    
    D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
    Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
                    Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
    
    pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                         Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                         sigmasqk = sigmaerror)
    K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
    # pvec <- .sampleP(N = N, avec = rep(1-sigma.py,(N-1)),
    #                   bvec = ((1:(N-1))*sigma.py + alpha.DP), K = K)  
    
    pvec <- .sampleP_PYM(N = N, alpha_val = alpha.py, sigma_val = sigma.py, K = K, Mat = matrixCnk, func =fun_pointer )  
    
    #alphaDP_g<- rgamma(1+N , 1/2 - log(pvec[N]))
    
  }
  
  list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
       sigmaerror = sigmaerror, rndEff = rndEff)
} 



.paramWrapper_3 <- function(REDUCT, inSamples,SS){   
  
  if(REDUCT){    
    
    function(CLUST, x,beta,Y,otherpar){
      
      N  <- otherpar$N
      r  <- otherpar$r
      D  <- otherpar$D
      Z  <- otherpar$Z
      sigmaerror <- otherpar$sigmaerror
      K          <- otherpar$K
      pvec       <- otherpar$pvec
      alpha_py   <- otherpar$alpha.PY
      sigma_py   <- otherpar$discount.PY
      tmp        <- .getPars_3(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta), 
                               D = D, Z = Z, sigmaerror = sigmaerror,
                               K = K, pvec = pvec, alpha.py = alpha_py,sigma.py=sigma_py,
                               inSamples = inSamples, SELECT = F)
      
      sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
                                   K = tmp$K, REDUCT=T))
      
      otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
                       sigmaerror = tmp$sigmaerror,
                       pvec = tmp$pvec, K = tmp$K, alpha.PY = alpha_py,discount.PY=sigma_py )
      
      return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
    }
    
  } else {
    
    function(CLUST, x, beta,Y,otherpar){
      
      sigmaDf  <- otherpar$sigmaDf
      XX  <- crossprod(x[inSamples,])
      IXX <- solveRcpp(XX)
      WX  <- crossprod(x[inSamples,], Y[inSamples,])
      WIX <- IXX%*%WX
      
      sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
                                   beta = beta, IXX = IXX, WX = WX, WIX = WIX,
                                   TRYPRIOR = T)$sigma
      otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
      
      return(list(sg = sg, otherpar = otherpar))
    }
  }
}


.getPars_3 <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                       alpha.py, sigma.py, inSamples,...){      
  
  # Y includes all terms but x%*%beta
  
  nn   <- length(inSamples)
  p    <- ncol(x)
  S    <- ncol(Y)
  ntot <- nrow(Y)
  nn   <- length(inSamples)
  
  covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
  rndEff <- RR%*%t(Z[K,])
  
  res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
  
  if(CLUST){   #only until convergence
    avec <- 1/rgamma(r, shape = (2 + r )/2, 
                     rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
    
    D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
    Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
                    Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
    
    pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                         Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                         sigmasqk = sigmaerror)
    K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
    pvec <- .sampleP(N = N, avec = rep(1-sigma.py,(N-1)),
                     bvec = ((1:(N-1))*sigma.py + alpha.py), K = K)  
     
  }
  
  list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
       sigmaerror = sigmaerror, rndEff = rndEff)
} 



# 
# 
# .getPars_4 <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
#                        alpha.PY,discount.PY, inSamples,rate,shape,ro.disc,alpha.PY_vec,...){      
#   
#   # Y includes all terms but x%*%beta
#   
#   nn   <- length(inSamples)
#   p    <- ncol(x)
#   S    <- ncol(Y)
#   ntot <- nrow(Y)
#   nn   <- length(inSamples)
#   
#   covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
#   z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
#   RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
#   if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
#   rndEff <- RR%*%t(Z[K,])
#   
#   res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
#   sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
#   
#   if(CLUST){   #only until convergence
#     avec <- 1/rgamma(r, shape = (2 + r )/2, 
#                      rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
#     
#     D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
#     Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
#                     Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
#     
#     pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
#                          Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
#                          sigmasqk = sigmaerror)
#     K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
#     
#     pvec <- .sampleP(N = N, avec = rep(1-discount.PY,(N-1)),
#                      bvec = ((1:(N-1))*discount.PY+alpha.PY), K = K)  
#     #pvec <- .sampleP(N=N, avec=rep(1,(N-1)),
#     #                bvec=rep(alpha.PY,(N-1)), K=K)
#     
#     alpha.PY<-metrop_PY_alpha(theta=alpha.PY,pvec=pvec,lik.fun=lik.alpha.fun,N=N,rate=rate,shape=shape,discount=discount.PY,alpha.PY_vec=alpha.PY_vec)
#     discount.PY<-metrop_PY_discount(theta=discount.PY,pvec=pvec,lik.fun=lik.disc.fun,ro.disc=ro.disc,N=N,alpha.PY=alpha.PY)
#     
#   }
#   
#   list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
#        sigmaerror = sigmaerror, rndEff = rndEff,alpha.PY=alpha.PY,discount.PY=discount.PY,rate,shape,ro.disc=ro.disc,alpha.PY_vec=c(alpha.PY_vec,alpha.PY))
# } 
# 
# 
# .paramWrapper_4 <- function(REDUCT, inSamples,SS){   
#   
#   if(REDUCT){    
#     
#     function(CLUST, x,beta,Y,otherpar){
#       
#       N  <- otherpar$N
#       r  <- otherpar$r
#       D  <- otherpar$D
#       Z  <- otherpar$Z
#       sigmaerror <- otherpar$sigmaerror
#       K          <- otherpar$K
#       pvec       <- otherpar$pvec
#       alpha.PY   <- otherpar$alpha.PY
#       discount.PY   <- otherpar$discount.PY
#       rate       <- otherpar$rate
#       shape      <- otherpar$shape
#       ro.disc    <- otherpar$ro.disc
#       alpha.PY_vec<-otherpar$alpha.PY_vec
#  
# 
#       
#       
#       tmp        <- .getPars_4(CLUST, x = x, N = N, r = r, Y = Y, B = t(beta), 
#                                D = D, Z = Z, sigmaerror = sigmaerror,
#                                K = K, pvec = pvec, alpha.PY = alpha.PY, discount.PY=discount.PY, shape=shape,rate=rate,
#                                ro.disc=ro.disc,alpha.PY_vec=alpha.PY_vec,
#                                inSamples = inSamples, SELECT = F)
#       
#       sg <- with(tmp, .expandSigma(sigma = tmp$sigmaerror, SS, Z = tmp$Z, 
#                                    K = tmp$K, REDUCT=T))
#       
#       otherpar <- list(A = tmp$A, N = N, r = r, D = tmp$D, Z = tmp$Z, 
#                        sigmaerror = tmp$sigmaerror,
#                        pvec = tmp$pvec, K = tmp$K, alpha.PY = tmp$alpha.PY,discount.PY=tmp$discount.PY,shape= shape,rate= rate,ro.disc=ro.disc,alpha.PY_vec=tmp$alpha.PY_vec)
#       
#       return(list(sg = sg, rndEff = tmp$rndEff, otherpar = otherpar))
#     }
#     
#   } else {
#     
#     function(CLUST, x, beta,Y,otherpar){
#       
#       sigmaDf  <- otherpar$sigmaDf
#       XX  <- crossprod(x[inSamples,])
#       IXX <- solveRcpp(XX)
#       WX  <- crossprod(x[inSamples,], Y[inSamples,])
#       WIX <- IXX%*%WX
#       
#       sg <- .updateWishartNoPrior( x[inSamples,], Y[inSamples,], sigmaDf,
#                                    beta = beta, IXX = IXX, WX = WX, WIX = WIX,
#                                    TRYPRIOR = T)$sigma
#       otherpar=list(Z = NA, K = NA, sigmaDf = sigmaDf)
#       
#       return(list(sg = sg, otherpar = otherpar))
#     }
#   }
# }

metrop_DP <- function(theta, #previous iteration alpha.DP
                      pvec,
                      lik.fun,
                      #prior.fun,
                      V=diag(theta), 
                      N, shape,rate,alpha.DP_vec
) {
  accept<-FALSE
  
  last.lik<-lik.fun(alpha=theta,pvec=pvec,N=N,shape=shape,rate=rate)
  last=theta
  if(length(alpha.DP_vec)<50) {ad_var=1}else{ad_var=(2.38^2)*var(alpha.DP_vec)}
  proposal <-rnorm(1,mean=theta,sd=sqrt(ad_var))
  if(proposal<0){return(last)
  }else{
  proposal.prior <-  dnorm(last,mean=proposal,sd=sqrt(ad_var),log=TRUE) #q(x,y)
  last.prior <-  dnorm(proposal,mean=last,sd=sqrt(ad_var),log=TRUE) #q(y,x)
  proposal.lik <- lik.fun(proposal,pvec,N,shape,rate)
  alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
  if (alpha > runif(1) & !is.nan(alpha) ) accept <- TRUE
  if (accept) {
    last <- proposal
  }
  return(last)
  }
}






.bisec<-function (f, a, b, num = 10, eps = 1e-05) 
{
  h = abs(b - a)/num
  i = 0
  j = 0
  a1 = b1 = 0
  while (i <= num) {
    a1 = a + i * h
    b1 = a1 + h
    if (f(a1) == 0) {
      print(a1)
      print(f(a1))
    }
    else if (f(b1) == 0) {
      print(b1)
      print(f(b1))
    }
    else if (f(a1) * f(b1) < 0) {
      repeat {
        if (abs(b1 - a1) < eps) 
          break
        x <- (a1 + b1)/2
        if (f(a1) * f(x) < 0) 
          b1 <- x
        else a1 <- x
      }
      print(j + 1)
      j = j + 1
      print((a1 + b1)/2)
      print(f((a1 + b1)/2))
    }
    i = i + 1
  }
  if (j == 0) {
    print("finding root is fail")
    return(0)
  }
  else return(x)
}


#likelihood function for MH steps

g_func<- function(alpha, sigma, N){
  alpha_vec1<- (0:(N-2))*sigma+ alpha +1
  alpha_vec2<- (1:(N-1))*sigma+ alpha
  gamma_vec<- gamma(alpha_vec1)/gamma(alpha_vec2)
  return(prod(gamma_vec))
}

lik.alpha.fun<-function(alpha,pvec,N,shape,rate,discount){
  if(alpha<0){ stop("alpha is negative!")
  }else{
    tmp<-sum(lgamma(alpha+1+discount*(c(1:(N-1))-1))-lgamma((alpha+discount*c(1:(N-1)))))+alpha*log(pvec[N])+(shape-1)*log(alpha)-rate*alpha
    #tmp<-g_func(alpha,discount,N)*pvec[length(pvec)]^(alpha)*alpha^(shape-1)*exp(-rate*alpha)
    return(tmp)
  }
}


lik.disc.fun<-function(discount,pvec,N,ro.disc,alpha){
  tmp<- (-N*lgamma(1-discount))+sum(lgamma(alpha+1+discount*(c(1:(N-1))-1))-lgamma((alpha+discount*c(1:(N-1)))))-discount*sum(log(pvec[1:(N-1)]))+discount*(N-1)*log(pvec[N])+log(ro.disc*ifelse(discount==0,1,0)+2*(1-ro.disc)*ifelse((discount<=0.5 & discount>0),1,0))
  #tmp<-(1/(gamma(1-discount)^N))*g_func(alpha,discount,N)*prod(pvec[1:(N-1)]^(-discount))*(pvec[length(pvec)]^(discount*(N-1)))*(ro.disc*ifelse(discount==0,1,0)+2*(1-ro.disc)*ifelse((discount<=0.5 & discount>0),1,0))
  return(tmp)
}

lik.alpha.DP.fun<-function(alpha,pvec,N,shape,rate){
  if(alpha<0){ stop("alpha is negative!")
  }else{
    tmp<-log(gamma(alpha)) - N*log(gamma(alpha/N)) + sum(((alpha/N)-1)*log(pvec)) + (shape-1)*log(alpha) - rate*alpha
    
    
    return(tmp)
  }
}


metrop_PY_alpha <- function(theta, #previous iteration alpha.DP
                            pvec,
                            lik.fun,
                            V=diag(theta), 
                            N, shape,rate,discount,alpha.PY_vec) {
  accept<-FALSE

  last.lik<-lik.fun(alpha=theta,pvec=pvec,N=N,shape=shape,rate=rate,discount=discount) 
  last=theta
  if(length(alpha.PY_vec)<50) {ad_var=2.38^2}else{ad_var=(2.38^2)*var(alpha.PY_vec)}
  proposal <-rnorm(1,mean=last,sd=sqrt(ad_var))
  if(proposal<0){return(last)
    }else{
  proposal.prior <-  dnorm(last,mean=proposal,sd=sqrt(ad_var),log=TRUE) #q(x,y)
  last.prior <-  dnorm(proposal,mean=last,sd=sqrt(ad_var),log=TRUE) #q(y,x)
  proposal.lik <- lik.fun(proposal,pvec,N,shape,rate,discount)
  alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
  if (alpha > runif(1) & !is.nan(alpha)) accept <- TRUE
  if (accept) {
    last <- proposal
  }
  return(last)
    }
}


metrop_PY_discount <- function(theta, #previous iteration alpha.DP
                               pvec,
                               lik.fun,
                               ro.disc,
                               V=diag(theta), 
                               N, alpha.PY) { 
  accept<-FALSE

  last.lik<-lik.fun(theta,pvec=pvec,N=N,ro.disc,alpha.PY) 
  last=theta
  #sample from proposal p(a)=1/2dirac(0)+1/2U[0,0.5]
  c<-rbinom(n=1, size=1,prob=0.5)
  if(c==1){proposal<-0}else{proposal<-runif(1,min=0,max=0.5)}
  dproposal<-function(x){if(x==0) {return(0.5)}else{return(1)}}
  proposal.prior <-  log(dproposal(last)) #q(x)
  last.prior <-  log(dproposal(proposal)) #q(y)
  
  proposal.lik <- lik.fun(proposal,pvec,N,ro.disc,alpha.PY)
  alpha <- exp(proposal.lik+proposal.prior-last.lik-last.prior)
  if (alpha > runif(1) & !is.nan(alpha)) accept <- TRUE
  if (accept) {
    last <- proposal
  }
  return(last)
}


.compute_tau_mean<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  gamma<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ alpha/(1-alpha) + theta))
  N<- N_eps*gamma
  return(N)
}



.compute_tau_mean_large_dim<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  log_val<- lgamma(1+theta) + lgamma(1+ theta/alpha + 1/(1-alpha) ) - lgamma(1+theta/alpha) - lgamma(1+ alpha/(1-alpha) + theta)
  #gamma<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ alpha/(1-alpha) + theta))
  gamma_val<- exp(log_val)
  N<- N_eps*gamma_val
  return(N)
}


.compute_tau_var<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  gamma2<- (gamma(1+theta)*gamma(1+ theta/alpha + 2/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (2*alpha)/(1-alpha) + theta))
  gamma1<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (1*alpha)/(1-alpha) + theta))
  gamma<-(gamma2-gamma1*gamma1)
  N<- (N_eps^2)*gamma
  return(sqrt(N))
}



.compute_tau_var_large_dim<- function(alpha,theta, eps=0.1){
  N_eps<- (eps/alpha)^(-alpha/(1-alpha))
  log_gamma_2<- lgamma(1+theta) + lgamma(1+ theta/alpha + 2/(1-alpha)) - lgamma(1+theta/alpha) - lgamma(1+ (2*alpha)/(1-alpha) + theta)
  gammaval_2<- exp(log_gamma_2)
  # gamma2<- (gamma(1+theta)*gamma(1+ theta/alpha + 2/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (2*alpha)/(1-alpha) + theta))
  log_gamma_1<- lgamma(1+theta) + lgamma(1+ theta/alpha + 1/(1-alpha)) - lgamma(1+theta/alpha) - lgamma(1+ (1*alpha)/(1-alpha) + theta)
  # gamma1<- (gamma(1+theta)*gamma(1+ theta/alpha + 1/(1-alpha)))/(gamma(1+theta/alpha)*gamma(1+ (1*alpha)/(1-alpha) + theta))
  gammaval_1<- exp(log_gamma_1)
  gamma_val<-(gammaval_2 -gammaval_1*gammaval_1)
  N<- (N_eps^2)*gamma_val
  return(sqrt(N))
}


##### Added functions to compute the hyperparameters for Ga(nu1,nu2) for alpha in DP and PY cases

#Function to define prior expected mean in number of groups for DP
funcDP<-function(x, S, K) {sum(x/(x+(1:S)-1))- K}
#Function to define prior expected mean in number of groups
funcPY<-function(x, S, K,sigma_py=0.25) {(x/sigma_py)*(prod((x+sigma_py+c(1:S) -1)/(x+c(1:S) -1))-1) - K}

#Function for sampling alpha
simulatuion_function_GD<- function(nu_ratio,ft,ns,Sn, K){
  nu2<-nu_ratio/20
  nu1<- nu2*nu_ratio
  alpha_s<- rgamma(ns, nu1,nu2)
  alpha_s_mod<- replace(alpha_s, alpha_s< 10^(-8), 10^(-8)) #to avoid small values for alpha, which could lead ti inf values in funcDP/PY
  sum_list<- sapply(alpha_s_mod,ft,S=Sn, K=K)
  return(mean(sum_list))
}

## function to obtain gamma parameters
gamma_parameters_for_K<- function(fn,Ktr,S_p,n_s){
  ratio<-.bisec(f=function(x) simulatuion_function_GD(x,ft=fn,ns=n_s,Sn=S_p,K=Ktr),0.01,1000, num=10)
  nu2<-ratio/20
  nu1<- ratio*nu2
  alpha_s<- rgamma(n_s, nu1,nu2)
  alpha_s_mod<- replace(alpha_s, alpha_s< 10^(-8), 10^(-8))
  sum_list<- sapply(alpha_s_mod,fn,S=S_p,K=Ktr)
  #plot(density(sum_list))
  alpha_fixed<-.bisec(f=function(x) fn(x,S=S_p,K=Ktr),0.01,1000)
  dif<- alpha_fixed- mean(alpha_s_mod)
  ## Version of function that contains warnings
  cat(dif," difference between fixed and simulated \n")
  if(abs(mean(sum_list))>0.05) cat(mean(sum_list),"big deviation! \n")
  else cat(mean(sum_list),"deviation \n")
  return(list(nu1=nu1, nu2=nu2))
}


##################### BNP functions

##Expectation for DP
functionDP<-function(x, n) {sum(x/(x+(1:n)-1))}
##Derivative of expectation for DP
functionDP_deriv<-function(x, n,K ) {sum(((1:n)-1)/(x+(1:n)-1)^2) -K}
### Expectation for DP multinomial
functionDPM<-function(x, n,N) {
  vec<- 0:(n-2)
  E<- N - (N-1)*(prod(x + 1 - x/N + vec)/(prod(x + 1 +vec)))
  return(E)
}

### Expectation for PY
functionPY<-function(x, n,sigma_py=0.25) {(x/sigma_py)*(prod((x+sigma_py+c(1:(n))-1)/(x+c(1:(n))-1))-1)}

### Compute variance for Pitman--Yor process

Var_PY <-  function(alpha, sigma, n) {
  if (n==1) {
    return(0)
  } else {
    El_prev=functionPY(alpha,n-1,sigma)
    exp_term<- (El_prev*((n -1)*sigma - alpha*sigma) + (n-1)*alpha -  sigma*sigma*((El_prev)^2)) /(n-1+ alpha)^2
    return (Var_PY(alpha, sigma,n-1)*(n-1+ alpha + 2*sigma)/(n-1+alpha) + exp_term)
  }
}

##### Simulation functions
#Function for sampling alpha
simulatuion_function_PY<- function(nu_ratio,variance=20,funct,ns,Sn){
  nu2<-nu_ratio/variance
  nu1<- nu2*nu_ratio
  alpha_s<- rgamma(ns, nu1,nu2)
  alpha_s_mod<- replace(alpha_s, alpha_s< 10^(-10), 10^(-10)) #to avoid small values for alpha, which could lead to inf values in funcDP/PY
  sum_list<- sapply(alpha_s_mod,funct,n=Sn)
  return(mean(sum_list, na.rm = TRUE))
}
simulatuion_function_DPM<- function(nu_ratio,variance=20,funct,ns,Sn,N_tr){
  nu2<-nu_ratio/variance
  nu1<- nu2*nu_ratio
  alpha_s<- rgamma(ns, nu1,nu2)
  alpha_s_mod<- replace(alpha_s, alpha_s< 10^(-10), 10^(-10)) #to avoid small values for alpha, which could lead to inf values in funcDP/PY
  sum_list<- sapply(alpha_s_mod,funct,n=Sn,N=N_tr)
  return(mean(sum_list, na.rm = TRUE))
}
#####

newton2 <- function(f,f_der, tol=1E-12,x0=1,N=50) {
  i <- 1; x1 <- x0
  p <- numeric(N)
  while (i<=N) {
    x1 <- (x0 - (f(x0)/f_der(x0)))
    p[i] <- x1
    i <- i + 1
    if (abs(x1-x0) < tol) break
    x0 <- x1
  }
  return(p[1:(i-1)])
}


compute_gamma_parameters<- function(fun,K,var_gamma=20){
  x<- seq(0.01,300,1)
  y=sapply(x, function(x) fun(x)) - K
  f_spline_smooth=smooth.spline(x, y) 
  roots <- newton2(f = function(x) predict(f_spline_smooth, x,deriv = 0)$y ,f_der=  function(x) predict(f_spline_smooth, x,deriv = 1)$y,x0=1,N=50)
  root<-  uniroot(function(x) predict(f_spline_smooth, x, deriv = 0)$y - 0, interval = c(0, 100))$root
  #print(root)
  nu2<- roots[length(roots)]/ var_gamma
  nu1<- roots[length(roots)]*nu2
  return(list(ratio=roots[length(roots)],nu1=nu1,nu2=nu2))
}


compute_fixed_parameters_1d<- function(fun,K){
  x<- seq(0.000001,300,0.1)
  y=sapply(x, function(x) fun(x)) - K
  f_spline_smooth=smooth.spline(x, y) 
  roots <- newton2(f = function(x) predict(f_spline_smooth, x,deriv = 0)$y ,f_der=  function(x) predict(f_spline_smooth, x,deriv = 1)$y,x0=1,N=50)
  #root<-  uniroot(function(x) predict(f_spline_smooth, x, deriv = 0)$y - 0, interval = c(0, 100))$root
  #print(roots)
  return(roots[length(roots)])
}

## using rootSolve and multiroot package!
compute_fixed_parameters_PY_2d<- function(K,V,n){
  model<- function(x, K,V,n){
    F1<- functionPY(x[1], n,x[2]) - K
    F2<- Var_PY(x[1], x[2],n) -V
    c(F1 = F1, F2 = F2)
  }
  roots_values <- multiroot(f = function(x) model(x,K,V,n), start=c(1,0.2), positive=TRUE)
  return(list(alpha = roots_values$root[1],sigma=roots_values$root[2]))
}


