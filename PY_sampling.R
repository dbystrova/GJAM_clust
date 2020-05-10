#library(tweedie)


library(frailtySurv)

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


A = Prior_on_K_PYM(0.18, 0.75, 50, 100,100000)

p = tibble(k=A$k,
           pPY=  A$p_k)%>%
  gather(Process_type, distribution, pPY)%>%
  ggplot(aes( x = k, y = distribution, color=Process_type )) +
  theme_bw() +
  geom_point() +
  geom_line()
p



##############################################################################################################

reccursive_Cnk<- function(n, k, sigma){
  if (k==0){
    if(n==0){
      return(1)
    }
    else{
      return(0)
    }
  }
  else{
    if (k>n){
      return(0)
    }
    else{
      return((n - 1 - sigma*k)*reccursive_Cnk(n-1, k, sigma) + sigma*reccursive_Cnk(n-1, k-1, sigma))
    }
  }
}





library(inline)
library(rbenchmark)


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




Cnk_Rcpp(39,13,0.2)

library(JuliaCall)



res <- benchmark(fibR(N),
                 fibRcpp(N),
                 columns=c("test", "replications", "elapsed",
                           "relative", "user.self", "sys.self"),
                 order="relative",
                 replications=1)
print(res)  ## show result



#a = reccursive_Cnk(50, 10, 0.5)

log_v_pdf<- function(v, alpha, sigma, nk_vec=c_tab, H){
  k<- length(nk_vec)
  sum<- 0
  for (j in (1:k)){
    val<- 0
    for (l in (1:nk_vec[j])){
      val0 = (v/H)^l * Cnk_Rcpp(nk_vec[j], l, sigma)
      val<- val + val0
    }
    sum<- sum + log(val)
  }
  result <- -v + (alpha/sigma -1)*log(v) + sum 
  
  return(result)
}

sigma=0.2
alpha=1 
H=111
ns=111
weights_PYM<- mult_PY(alpha,sigma,H)
c <- sample(1:H,size=ns, replace=TRUE, prob=weights_PYM )
c_tab =table(c)

x = log_v_pdf(50, 1.0, 0.2,c_tab, H=112)
v_vec<- seq(0.5, 112, 1)

x<- sapply(v_vec, log_v_pdf,alpha=1, sigma=0.2, nk_vec=c_tab, H=112)

plot(v_vec, exp(x)/sum(exp(x)))

bound1_vec<- sqrt(exp(x)) #10
bound1<- max(sqrt(x)) #10

plot(v_vec, bound1_vec)
bound2_vec<- x*sqrt(x) #650
bound2<- max(x*sqrt(x)) #650
plot(v_vec, bound2_vec)


i=1
x_vec<- c()
while (i<=100){
  u<- runif(1, 0, bound1)
  v <- runif(1, 0, bound2 )
 if (u <= sqrt(exp(log_v_pdf(v/u, alpha=1, sigma=0.2, nk_vec=c_tab, H=50)))){
   x_vec[i]<-v/u 
   i=i+1
 }
}

plot(density(x_vec))
plot(v_vec, bound1)


library(rust)
x <- ru(logf =log_v_pdf,alpha=1, sigma=0.2, nk_vec=c_tab, H=20, n=90,  d=1, init=1)


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

m <- sample_lk(c_tab,x1$sim_vals[1],sigma,H)
print(m)


sigma=0.2
alpha=1 
H=111
ns=111
weights_PYM<- mult_PY(alpha,sigma,H)
c <- sample(1:H,size=ns, replace=TRUE, prob=weights_PYM )
c_tab =table(c)

x <- ru(logf =log_v_pdf,alpha=alpha, sigma=sigma, nk_vec=c_tab, H=H, n=1,  d=1, init=1)
sample_lk(c_tab,x1$sim_vals[1],sigma, H)


alpha <- 0.1
max_phi <- qgamma(0.999, shape = alpha)
ptr_gam <- create_xptr("logdgamma")
lambda <- find_lambda_one_d_rcpp(logf = ptr_gam, alpha = alpha,
                                 max_phi = max_phi)
# Box-Cox transformation parameter
lambda$lambda
#> [1] 0.06758891
gam <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = 1000, trans = "BC",
               lambda = lambda)


plot(gam, xlab = "x")
plot(gam, ru_scale = TRUE, xlab = "y")

# Jumps from TS(1/H, sigma, U)
# U ~ Ga(alpha/sigma, 1)

## rescaled 
# J_h| U ~ TS(alpha/(sigma*H),sigma, U ) , U^(sigma)~ Ga(alpha/sigma, alpha/sigma) 









ratioU <- function(nvals)
{
  h_x = function(x) exp(-x)
  # u- is b-, u+ is b+ and v+ is a in the example:
  uminus = 0
  uplus = 2/exp(1)
  vplus = 1
  X.vals <- NULL
  i <- 0
  repeat {
    i <- i+1
    u <- runif(1,0,vplus)
    v <- runif(1,uminus,uplus)
    X <- u/v
    if(v^2 <= h_x(X)) {
      tmp <- X
    }
    else {
      next
    }
    X.vals <- c(X.vals,tmp)
    if(length(X.vals) >= nvals) break
  }
  answer <- X.vals  
  answer
}

sol = ratioU(1000) 
par(mfrow=c(1,2))
hist(sol,breaks=50, main= "using ratioU",freq=F)
hist(rexp(1000),breaks = 50, main="using rexp from R",freq=F)
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot(density(sol))
plot(density(rexp(1000)))
par(mfrow=c(1,1))













library(JuliaCall)
julia <- julia_setup()
#> Julia version 1.0.3 at location /Applications/Julia-1.0.app/Contents/Resources/julia/bin will be used.
#> Loading setup script for JuliaCall...
#> Finish loading setup script for JuliaCall.

## If you want to use `Julia` at a specific location, you could do the following:
## julia_setup(JULIA_HOME = "the folder that contains Julia binary").
## You can also set JULIA_HOME in command line environment or use `options(...)`.

## Different ways of using Julia to calculate sqrt(2)

# julia$command("a = sqrt(2);"); julia$eval("a")
julia_command("a = sqrt(2);"); julia_eval("a")
#> [1] 1.414214
julia_eval("sqrt(2)")
#> [1] 1.414214
julia_call("sqrt", 2)
#> [1] 1.414214
julia_eval("sqrt")(2)
#> [1] 1.414214
julia_assign("x", sqrt(2)); julia_eval("x")
#> [1] 1.414214
julia_assign("rsqrt", sqrt); julia_call("rsqrt", 2)
#> [1] 1.414214
2 %>J% sqrt
#> [1] 1.414214





