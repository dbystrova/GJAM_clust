#library(tweedie)

library(Rcpp)
library(RcppArmadillo)
library(arm)
library(inline)
library(rbenchmark)
library(frailtySurv)
library(rust)
library(gtools)
source("rlaptrans.r")
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
# p = tibble(k=A$k,
#            pPY=  A$p_k)%>%
#   gather(Process_type, distribution, pPY)%>%
#   ggplot(aes( x = k, y = distribution, color=Process_type )) +
#   theme_bw() +
#   geom_point() +
#   geom_line()
# p



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




Cnk_Rcpp <- cxxfunction(signature(ns="integer", ks="integer", sigmas ="numeric"),
                        plugin="Rcpp",
                        incl=incltxt,
                        body='
                        int n = Rcpp::as<int>(ns);
                        int k = Rcpp::as<int>(ks);
                        float sigma= Rcpp::as<float>(sigmas);
                        return Rcpp::wrap( reccursive_Cnk(n,k,sigma));
                        ')




#Cnk_Rcpp(39,10,0.2)




#a = reccursive_Cnk(50, 10, 0.5)

log_v_pdf<- function(v, alpha, sigma, nk_vec=c_tab, H){
  k<- length(nk_vec)
  sum<- 0
  for (j in (1:k)){
    val<- 0
    for (l in (1:nk_vec[j])){
      val0 = (v/H)^l * Cnk_Rcpp(nk_vec[j], l, sigma)
      #val0 = exp(l*log(v) + log(Cnk_Rcpp(nk_vec[j], l, sigma))  - l*log(H))
      val<- val + val0
    }
    sum<- sum + log(val)
  }
  result <- -v + (alpha/sigma -1)*log(v) + sum 
  
  return(result)
}



log_v_pdf_matrix<- function(v, alpha, sigma, nk_vec=c_tab, H, Mat){
  k<- length(nk_vec)
  sum<- 0
  for (j in (1:k)){
    val<- 0
    for (l in (1:nk_vec[j])){
      val0 = exp(l*log(v) + Mat[nk_vec[j],l])
      val<- val + val0
      print(Mat[nk_vec[j],l])
    }
    sum<- sum + log(val)
  }
  result <- -v + (alpha/sigma -1)*log(v) + sum 
  
  return(result)
}


v_pdf_matrix<- function(v, alpha, sigma, nk_vec=c_tab, H, Mat){
  k<- length(nk_vec)
  sum<- 0
  for (j in (1:k)){
    val<- 0
    for (l in (1:nk_vec[j])){
      val0 = exp(l*log(v) + Mat[nk_vec[j],l])
      val<- val + val0
    }
    sum<- sum + log(val)
  }
  result <- -v + (alpha/sigma -1)*log(v) + sum 
  
  return(exp(result))
}

log_v_pdf_grad<- function(v, alpha, sigma, nk_vec=c_tab, H){
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

Rcpp::cppFunction('arma::mat timesThree(arma::mat x) {return x * 3.0;}', depends = c("RcppArmadillo","BH"))
timesThree(matrix(c(1,7,14),3))


load("Cnk_mat_112.Rdata")


cppFunction('double log_v_pdf_mat(const Rcpp:: NumericVector& x, const Rcpp::List& pars){
  double alpha = pars["alpha"];
            double sigma = pars["sigma"];
            int k = pars["k"];
            Rcpp::NumericVector nk_vec = pars["nk_vec"];
            arma::mat Cnk_mat =  as<arma::mat>(pars["Cnk_mat"]);
            double sum = 0;
            double val0;
            for (int j=0; j<k;++j){
            Rcpp::Rcout << "nk matrix is" << std::endl << nk_vec[j] << std::endl;
            double val = 0;
            for (int l=1; l<(nk_vec[j]+1); ++l){
             val0 = exp(l*log(x[0]) + Cnk_mat(nk_vec[j]-1,l-1));
            //val0 = pow(x[0], l) * exp(Cnk_mat(nk_vec[j]-1,l-1));
            val =val+ val0;
            Rcpp::Rcout << "Cnk value" << std::endl <<  Cnk_mat(nk_vec[j]-1,l-1) << std::endl;
            }
            sum  = sum +log(val);
            }
            return ( -x[0] + (alpha/sigma -1)*log(x[0]) + sum);
            }' ,depends = c("RcppArmadillo","BH"))

pars_v = list(alpha =alpha, sigma=sigma,k =length(c_tab),nk_vec=c_tab, Cnk_mat =Cnk_112_112)
log_v_pdf_mat(1, pars_v)
log_v_pdf_matrix(1,alpha=alpha, sigma=sigma, nk_vec=c_tab, H=H, Mat=Cnk_112_112)

sigma=0.25
alpha=4
H=112
ns=112
weights_PYM<- mult_PY(alpha,sigma,H)
c <- sample(1:H,size=ns, replace=TRUE, prob=weights_PYM )
c_tab =table(c)

v_vec<- seq(0.5, 2*H, 1)


x2<- sapply(v_vec, log_v_pdf_matrix,alpha=alpha, sigma=sigma, nk_vec=c_tab, H=H, Mat=Cnk_112_112)

x2_s<- sapply(v_vec, log_v_pdf_mat,pars = list(alpha =alpha, sigma=sigma,k =length(c_tab),nk_vec=c_tab, Cnk_mat =Cnk_112_112))

#x3<- sapply(v_vec, v_pdf_matrix,alpha=alpha, sigma=sigma, nk_vec=c_tab, H=H, Mat=Cnk_112_112)

#log_v_pdf_matrix(1, alpha, sigma, nk_vec=c_tab, H, Mat=M)

plot(v_vec, exp(x2)/sum(exp(x2)))
lines(v_vec, exp(x2_s)/sum(exp(x2_s)))

#lines(v_vec, exp(x), col="2")


library(rust)
ptr_N01 <- create_xptr("log_v_pdf_C")
ptr_logv_mat <- create_xptr("log_v_pdf_mat")
ptr_logv_comp_mat <- create_xptr("log_v_pdf_comp_mat")

x_new <- ru(logf =log_v_pdf_matrix,alpha=4, sigma=0.25, nk_vec=c_tab, H=H,Mat=Cnk_112_112, n=100,  d=1, init=1)


#new2 = ru_rcpp(logf = ptr_N01,alpha=alpha, sigma=sigma,H=H,k = length(c_tab), nk_vec=c_tab, n=100,  d=1, init=1)
#new = ru_rcpp(logf = ptr_logv_mat,alpha=alpha, sigma=sigma,H=H,k = length(c_tab), nk_vec=c_tab,Cnk_mat=Cnk_112_112, n=1000,  d=1, init=1)

#new3 = ru_rcpp(logf = ptr_logv_mat,alpha=alpha, sigma=sigma,H=H,k = length(c_tab), nk_vec=c_tab,Cnk_mat=Cnk_112_112, n=1000,  d=1, init=1)
new3 = ru_rcpp(logf = ptr_logv_comp_mat,alpha=alpha, sigma=sigma,H=H,k = length(c_tab), nk_vec=c_tab,Cnk_mat=Cnk_112_112, n=1000,  d=1, init=1)

#plot(density(new$sim_vals))
plot(v_vec, exp(x2)/sum(exp(x2)), col="red")
#lines(density(new$sim_vals))
lines(density(x_new$sim_vals))
lines(density(new3$sim_vals))


plot(v_vec, exp(x2)/sum(exp(x2)))
lines(x_new$sim_vals, col="red")
lines(new3$sim_vals, col="b")



##########

ptr_N01 <- create_xptr("logdN01")

# Use ru and ru_rcpp starting from the same random number seed and check
# that the simulated values are the same.
n=100
set.seed(47)
x_old <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = n, init = 0.1)
head(x_old$sim_vals)
#>            [,1]
#> [1,]  0.7764728
#> [2,]  0.5310434
#> [3,] -0.1046049
#> [4,]  1.2111509
#> [5,]  1.1391379
#> [6,]  0.5180914

logf = function(x) -x ^ 2 / 2
y_vec<- seq(-20, 20, 1)

y_v<- sapply(y_vec,logf)
x_new <- ru_rcpp(logf = ptr_N01, d = 1, n = n, init = 0.1)

#plot(y_vec, exp(y_v)/sum(exp(y_v)), col="red")
lines(density(x_old$sim_vals))
#lines(density(x_new$sim_vals))




set.seed(47)
x_new <- ru_rcpp(logf = ptr_N01, d = 1, n = n, init = 0.1)
head(x_new$sim_vals)


########













sigma=0.25
alpha=4
H=112
ns=50
weights_PYM<- mult_PY(alpha,sigma,H)
c <- sample(1:H,size=ns, replace=TRUE, prob=weights_PYM )
c_tab =table(c)

v_vec<- seq(0.5, 2*H, 1)

x<- sapply(v_vec, log_v_pdf,alpha=alpha, sigma=sigma, nk_vec=c_tab, H=H)

x2<- sapply(v_vec, log_v_pdf_matrix,alpha=alpha, sigma=sigma, nk_vec=c_tab, H=H, Mat=M)

x3<- sapply(v_vec, v_pdf_matrix,alpha=alpha, sigma=sigma, nk_vec=c_tab, H=H, Mat=Cnk_112_112)

#log_v_pdf_matrix(1, alpha, sigma, nk_vec=c_tab, H, Mat=M)

plot(v_vec, exp(x2)/sum(exp(x2)))
#lines(v_vec, exp(x), col="2")

plot(v_vec, x3)

library(rust)

ptr_logv_mat <- create_xptr("log_v_pdf_mat")

#x <- ru(logf =log_v_pdf,alpha=4, sigma=0.25, nk_vec=c_tab, H=50, n=10,  d=1, init=1)
ptr_N01 <- create_xptr("log_v_pdf_C")


#new = ru_rcpp(logf = ptr_N01,alpha=alpha, sigma=sigma,H=H,k = length(c_tab), nk_vec=c_tab, n=10,  d=1, init=1)

new = ru_rcpp(logf = ptr_logv_mat,alpha=alpha, sigma=sigma,H=H,k = length(c_tab), nk_vec=c_tab,Cnk_mat=Cnk_112_112, n=1000,  d=1, init=1)

plot(density(new$sim_vals))
lines(v_vec, exp(x2)/sum(exp(x2)), col="red")




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










bound1_vec<- sqrt(exp(x)) #10
bound1<- max(sqrt(exp(x))) #10

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

library(multicool)
library(CryptRndTest)

load("IJulia_part/Cnk_mat_112_H025.Rdata")
PY_prior<- function(k,H, n, alpha, sigma, Cnk_mat){
   n_vec<- 0:(n-2)
  fal_fact <- prod(alpha +1 +n_vec)
  coef = exp(log(factorial(H)) -  log(factorial(H - k)) - log(fal_fact))
  sum<- 0
  for (l in (k:n)){
    if (k==l){
      val0 = exp( lgamma(alpha/sigma +l ) - lgamma(alpha/sigma + 1) - l*log(H)  + Cnk_mat[n,l])
    }
    else{
    val0 = exp( lgamma(alpha/sigma +l ) - lgamma(alpha/sigma + 1) - l*log(H) + Strlng2(l, k, log = TRUE) + Cnk_mat[n,l])
    }
     sum<- sum + val0
  }
 return( (coef*sum)/sigma)
}

exp_var_PYM<- function(alpha,sigma,H=H, n=n, Mat_prior){
  x_vec<- 1:H
  pks<- sapply(x_vec, PY_prior,H=H, n=n, alpha=alpha, sigma=sigma,Cnk_mat=Mat_prior)
  plot(x_vec, pks)
  Exp <- sum(pks*x_vec)
  Var<- sum(((x_vec- Exp)^2)*pks)
 return(list(E= Exp, V=Var ))
}




compute_alpha_PYM<- function(H,n,sigma, Mat_prior, K){
  x<- seq(0.001,300,0.3)
  y=sapply(x, function(x) exp_var_PYM(x, sigma=sigma,H=H, n=n,Mat_prior=Mat_prior)$E) - K
  f_spline_smooth=smooth.spline(x, y) 
  roots <- newton2(f = function(x) predict(f_spline_smooth, x,deriv = 0)$y ,f_der=  function(x) predict(f_spline_smooth, x,deriv = 1)$y,x0=1,N=50)
  #root<-  uniroot(function(x) predict(f_spline_smooth, x, deriv = 0)$y - 0, interval = c(0, 100))$root
  #print(roots)
  return(roots[length(roots)])
}


par = compute_alpha_PYM(H=112,n=112,sigma=0.25,Mat_prior= Cnk_112_112_H025, K=16)

x_vec<- 1:112
pks<- sapply(x_vec, PY_prior,H=112, n=112, alpha=2.577792, sigma=0.25,Cnk_mat=Cnk_112_112_H025)
plot(x_vec, pks)
Exp <- sum(pks*x_vec)
Exp
Var<- sum(((x_vec- Exp)^2)*pks)
#PY_prior(2,H=112, n=112, alpha=0.88, sigma=0.25,Cnk_mat=Cnk_112_112_025)


x_vec<- 1:112

pks<- sapply(x_vec, PY_prior,H=112, n=112, alpha=0.4827001, sigma=0.5,Cnk_mat=Cnk_112_112_H05)
plot(x_vec, pks)

exp<-sum(x_vec*pks)
exp



tibble(it= 1: length(apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x)))),
       DP= apply(fit_gjam$chains$kgibbs,1,function(x) length(unique(x))),
       DP2 =apply(fit_gjamDP1$chains$kgibbs,1,function(x) length(unique(x))),
       PY1=apply(fit_gjamPY1$chains$kgibbs,1,function(x) length(unique(x))),
       PY2=apply(fit_gjamPY2$chains$kgibbs,1,function(x) length(unique(x)))
) %>%
  gather(Model, trace, DP:PY2)%>%
  ggplot(aes(x=it,y=trace,col=Model))+geom_line(alpha=0.8)+ scale_color_viridis(discrete=TRUE)+
  labs(title="Traceplots of the posterior of the number of clusters")+xlab("iterations")+ylab("Number of clusters") +theme_bw()+geom_hline(yintercept = 16,color = "red")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 15),legend.position = "top", plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 20)) +theme(legend.text=element_text(size=15))



