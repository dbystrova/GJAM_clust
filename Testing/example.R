

library(inline)
library(rbenchmark)
library(frailtySurv)
library(rust)

library("gtools")

Rcpp::sourceCpp('src/user_fns.cpp')




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




Prior_on_K_PYM<- function(alpha, sigma, H,ns, runs=10^4 ){
  array_nc<-c()
  i=1
  while (i<=runs){
    weights_PYM<- mult_PY(alpha,sigma,H)
    c <- sample(1:H,size=ns, replace=TRUE, prob=weights_PYM )
    n_c<- length(unique(c))
    array_nc[i]<- n_c
    i=i+1
  }
  p_k = tibble(k=as.numeric(names(table(array_nc))), 
               p_k=as.vector(table(array_nc))/sum(table(array_nc)))
  p_zeros= tibble(k=(1:H)[!1:H%in%p_k$k], 
                  p_k=rep(0,length((1:H)[!1:H%in%p_k$k])))
  return(rbind(p_k, p_zeros))
}


#A = Prior_on_K_PYM(0.18, 0.75, 50, 100,100000)
# 
# p = tibble(k=A$k,pPY=  A$p_k)%>%
#   gather(Process_type, distribution, pPY)%>%
#   ggplot(aes( x = k, y = distribution, color=Process_type )) +
#   theme_bw() +
#   geom_point() +
#   geom_line()
# p


incltxt <- '
double reccursive_Cnk(const int n, const int k, const double sigma) {
if (k == 0) {
if ( n == 0) {
return(1);
}
else{ 
return (0);
} 
} 
else{
if (k>n){
return(0);
}
else{
return((n - 1 - sigma*k)*reccursive_Cnk(n-1, k, sigma) + sigma*reccursive_Cnk(n-1, k-1, sigma));
}
}
}'



Cnk_Rcpp <- cxxfunction(signature(ns="integer", ks="integer", sigmas ="numeric"),
                        plugin="Rcpp",
                        incl=incltxt,
                        body='
                        int n = Rcpp::as<int>(ns);
                        int k = Rcpp::as<int>(ks);
                        float sigma= Rcpp::as<float>(sigmas);
                        return Rcpp::wrap( reccursive_Cnk(n,k,sigma));
                        ')



pdf_lk<- function(l,v, n_k, sigma,H){
  return( ((v/H)^l)*Cnk_Rcpp(n_k,l,sigma))
}

sample_lk<- function(nk_vec,v,sigma,H){
  l_post<-c()
  k<- length(nk_vec)
  for (i in 1:k){
    l_vec<- 1:nk_vec[i]
    if (length(l_vec)==1){
      l_post[i]=l_vec
    }
    else{
      p_v<- sapply(l_vec, function(x) pdf_lk(x,v,nk_vec[i],sigma,H))
      pv_norm<- p_v/sum(p_v)
      l_post[i]<- sample(1:(nk_vec[i]),size=1, replace=TRUE, prob=pv_norm)
    }
  }
  return(l_post)
}

######## Example from Baysian nonparametric data analysis 


read.dta <- function()
{
  X <- read.table("EIG.txt",header=1)
  y <- X$recorded
  y <- log( y[!is.na(y)] )  # log EIG121 expression
  n <- length(y)
  return(dta=list(y=y, n=n))
}

## hyperparameters
a <- 1;   b <- 1     # 1/sig ~ Ga(a,b)
m0 <- -3;  B0 <- 4    # G0 = N(m0,B0)
M <- 1
alpha=1
sigma=0.25
H <- 10

# ##################################################################
# initialize clustering..
# ##################################################################

init.DPk <- function()
{ ## inital EDA estimate of G = sum_{h=1..10} w_h delta(m_h)
  ## returns:
  ##   list(mh,wh)
  ## use (mh,wh) to initialize the blocked Gibbs
  
  ## cluster data, and cut at height H=10, to get 10 clusters
  hc <- hclust(dist(y)^2, "cen")
  r  <- cutree(hc, k = 10)
  ## record cluster specific means, order them 
  mh1 <- sapply(split(y,r),mean)    # cluster specific means == m_h
  wh1 <- table(r)/n
  idx <- order(wh1,decreasing=T)    # re-arrange in deceasing order
  mh <- mh1[idx]
  wh <- wh1[idx]
  return(list(mh=mh,wh=wh,r=r))
}   


# ##################################################################
# 2. Blocked GS
# ##################################################################


gibbs.H <- function(n.iter=500)
{
  
  G <- init.DPk()
  sig <- 0.11
  
  
  ## data structures to save imputed F ~ p(F | ...)
  #xgrid <- seq(from= -10, to=2,length=50)
  xgrid <- seq(from= -4, to=4,length=50)
  fgrid <- NULL
  plot(density(y),xlab="X",ylab="Y",bty="l",type="l",
       xlim=c(-4, 4),ylim=c(0,0.4), main="")
  #-8, 2
  ## Gibbs
  for(iter in 1:n.iter){
    G$r <-  sample.r(G$wh,G$mh,sig)   # 1. r_i ~ p(r_i | ...), i=1..n
    G$mh <- sample.mh(G$wh,G$r,sig)   # 2. m_h ~ p(m_h | ...), h=1..H
    #G$vh <- sample.vh(G$r)            # 3. v_h ~ p(v_h | ...), h=1..H
    G$vh <- sample.vh2(G$r)       
    th <- G$mh[G$r]                   # record implied th[i] = mh[r[i]]
    sig <- sample.sig(th)       # 4. sig ~ p(sig | ...)
    
    ## record draw F ~ p(F | th,sig,y) (approx)
    f   <- fbar.H(xgrid,G$wh,G$mh,sig)
    lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
  }
  ## add overall average (= posterior mean) to the plot
  fbar <- apply(fgrid,2,mean)
  lines(xgrid,fbar,lwd=3,col=2)
  return(fgrid)
}

sample.r <- function(wh,mh,sig)
{ ## samle allocation indicators
  
  r <- rep(0,n)
  for(i in 1:n){
    ph <-   dnorm(y[i],m=mh,sd=sig)*wh # likelihood   * prior
    ## p(yi | ri=h) * w_h
    r[i] <- sample(1:H,1,prob=ph)
  }
  return(r)
}


sample.mh <- function(wh,r,sig)
{ ## sample mh ~ p(mh | ...)
  ##
  
  mh <- rep(0,H)     # initialize
  for(h in 1:H){
    if(any(r==h)){      # some data assigned to h-th pointmass
      Sh <- which(r==h) # Sh = {i: r[i]=h
      nh <- length(Sh)
      ybarh <- mean(y[Sh])
      varh   <- 1.0/(1/B0 + nh/sig^2)
      meanh  <- varh*(1/B0*m0 + nh/sig^2*ybarh)
    } else {            # no data assinged to h-th pointmass
      varh  <- B0       # sample from base measure
      meanh <- m0
    }
    mh[h] <- rnorm(1,m=meanh,sd=sqrt(varh))
  }
  return(mh)
}



### This is DP multinomial variant


sample.vh <- function(r)
{## sample vh ~ p(vh | ...)
  ## returns: wh
  
  vh <- rep(0,H)  # initialize
  wh <- rep(0,H)
  V <-  1         # record prod_{g<h} (1-vh_h)
  for(h in 1:(H-1)){
    Ah <- which(r==h)
    Bh <- which(r>h)
    vh[h] <-  rbeta(1, 1+length(Ah), M+length(Bh))
    wh[h] <- vh[h]*V
    V <- V*(1-vh[h])
  }
  vh[H] <- 1.0
  wh[H] <- V
  return(wh)
}

## This is PY multinomial version


sample.vh2 <- function(r){
  ## sample vh ~ p(vh | ...)
  ## returns: wh
  #r=c
  n_k<- table(r)
  lh<-  rep(0,H)
  #sample v
  v_s = ru_rcpp(logf = ptr_N01,alpha=alpha, sigma=sigma,H=H,k = length(n_k), nk_vec=n_k, n=1,  d=1, init=1)
  #sample lk
  lk <- sample_lk(n_k,v_s$sim_vals[1],sigma,H)
  lh[c(as.numeric(c(names(n_k))))]= lk
  vh <- rep(0,H)  # initialize
  W_h <- rep(0,H)
  P_h<-  rep(0,H)
  p_vec<- n_k - lk*sigma
  W_h<- rdirichlet(1,c(p_vec, sum(lk)*sigma + alpha))
  ### R
  alpha_post<- alpha + sum(lk)*sigma 
  Uv<- rgamma(1,alpha_post/sigma,alpha_post/sigma)
  U<- (Uv)^(1/sigma)
  x.rlap <- rlaptrans(H, lt.temp_st_pdf, c=alpha_post/(sigma*H), sigma, k=U)
  R_h<- x.rlap /sum(x.rlap)
  P_h[c(as.numeric(c(names(n_k))))]<- W_h[1:length(n_k)] + W_h[length(n_k)]* R_h[1:length(n_k)]
  P_h[-c(as.numeric(c(names(n_k))))] <-  W_h[length(n_k)]* R_h[(length(n_k)+1):H]
  return(P_h)
}






fbar.H <- function(xgrid,wh,mh,sig)
{ ## return a draw F ~ p(F | ...) (approx)
  
  fx <- rep(0,length(xgrid))
  for(h in 1:H)
    fx <- fx + wh[h]*dnorm(xgrid,m=mh[h],sd=sig)
  return(fx)
}


sample.sig <- function(th)
{ ## sample
  ##   sig ~ p(sig | ...)
  ## returns: sig
  
  s2 <- sum( (y-th)^2 )    # sum of squared residuals
  a1 <- a+0.5*n
  b1 <- b+0.5*s2
  s2.inv <- rgamma(1,shape=a1,rate=b1)
  return(1/sqrt(s2.inv))
}

plt.all <- function(fgrid,sim=T,dens=T)
{
  #xgrid <- seq(from= -10, to=2,length=50)
  xgrid <- seq(from= -4, to=4,length=50)
  M <- nrow(fgrid)
  idx0 <- 21:M
  fgrid0 <- fgrid[idx0,]            # drop initial transient
  idx1 <- which(idx0 %% 5 == 0 )    # thin out for plotting
  fgrid1 <- fgrid[idx1,]
  fbar <- apply(fgrid0,2,mean)
  plot(xgrid,fbar,xlab="log(EIG121)", ylab="G",
       type="l",lwd=3,bty="l",ylim=c(0,0.35))
  if(sim){
    matlines(xgrid,t(fgrid1),col=1)
    ## matlines(xgrid,t(fgrid0),col=2)
    lines(xgrid,fbar,type="l",lwd=4,col="grey")
  }
  if (dens){
    lines(density(y),col="yellow",lty=2,lwd=4)
  }
}

plt.dta <- function()
{
  hist(y,main="",xlab="log(EIG121)",ylab="FREQ",prob=T)
}


ex <- function()
{
  dta <- read.dta()
  attach(dta)
  
  ## run MCMC
  fgrid <- gibbs.H()
  plt.dta()
  plt.all(fgrid)
}


#ex()




ex2 <- function()
{
  c<- c(1/4,1/8,1/4,1/8,1/4)
  mu<-c(-2,-1, 0, 1, 2)
  sd<- c(0.2,0.2, 0.2,0.2, 0.2)
  plot(1:5,c, ylim=c(0,1))
  N<-50
  components<- sample(1:5, prob=c, size=N, replace=TRUE)
  mu<-c(-2,-1, 0, 1, 2)
  sigma_vec<- sd
  
  samples<- rnorm(n=N, mean=mu[components],sd=sigma_vec[components])
  plot(density(samples))
  dta$y <- samples
  dta$n<- N
  attach(dta)
  y <- samples
  n <- N
  ## run MCMC
  fgrid <- gibbs.H()
  plt.dta()
  plt.all(fgrid)
}

ex2()
