########## BNP functions
########## Compute Expected Number of clusters

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


### With hyperparameters

#DP_par<- compute_gamma_parameters(fun=function(x) simulatuion_function_PY(x,funct=functionDP,ns=30000,Sn=111), K=16)
#DP_mult_par<- compute_gamma_parameters(fun=function(x) simulatuion_function_DPM(x,funct=functionDPM,ns=30000,Sn=111,N_tr = 111), K=16)
#PY_par<- compute_gamma_parameters(function(x) simulatuion_function_PY(x,funct=functionPY,ns=30000,Sn=111), K=16)
#PY_fixed<- compute_fixed_parameters_1d(fun= function(x) functionPY(x, 111,sigma_py=0.25),K=8)


### With hyperparameters
#############################################################################


StickBreakingPY <- function(alpha,sigma, N) {
  a<- rep(1-sigma,N-1)
  b<- alpha+ sigma*(1:N-1)
  V <- rbeta(N-1 , a, b)
  p    <- vector("numeric",length=N)
  p[1] <- V[1]
  for(l in 2:(N - 1))p[l] <- prod(1 - V[1:(l - 1)])*V[l]
  p[N] <- prod(1 - V)   
  p
  return(p)
}

Prior_on_K_SB<- function(alpha, sigma, N_s,N_tr, runs=10^4 ){
  array_nc<-c()
  i=1
  while (i<=runs){
    weights_PY<- StickBreakingPY(alpha,sigma,N_tr)
    c <- sample(1:N_tr,size=N_s, replace=TRUE, prob=weights_PY )
    n_c<- length(unique(c))
    array_nc[i]<- n_c
    i=i+1
  }
  p_k = tibble(k=as.numeric(names(table(array_nc))), 
               p_k=as.vector(table(array_nc))/sum(table(array_nc)))
  p_zeros= tibble(k=(1:N_tr)[!1:N_tr%in%p_k$k], 
                  p_k=rep(0,length((1:N_tr)[!1:N_tr%in%p_k$k])))
  pk_padded= rbind(p_k, p_zeros)
  E_k=sum((1:N_tr) *(pk_padded$p_k))
  V_k = sum((((1:N_tr) - E_k)^2) *(pk_padded$p_k))
  return(list(pk=pk_padded$p_k, Ek= E_k,Vk=V_k ))
}




compute_parameters_SB_1d<- function(K,n,N_tr,ns=10^5){
  x<- seq(1,10,0.1)
  y=sapply(x, function(x) Prior_on_K_SB(x, 0.25, n,N_tr, runs=ns)$Ek) - K
  f_spline_smooth=smooth.spline(x, y) 
  roots <- newton2(f = function(x) predict(f_spline_smooth, x,deriv = 0)$y ,f_der=  function(x) predict(f_spline_smooth, x,deriv = 1)$y,x0=1,N=50)
  root<-  uniroot(function(x) predict(f_spline_smooth, x, deriv = 0)$y - 0, interval = c(0, 5))$root
  print(root)
  return(roots[length(roots)])
}

#par <- compute_parameters_SB_1d(16,112,112,10^4)
#PkSB<- Prior_on_K_SB(par,0.25,112,112, runs=10^5 )
#PkSB


#p = tibble(k=1:112,
#           pPY=  PkSB$pk)%>%
#  gather(Process_type, distribution, pPY)%>%
#  ggplot(aes( x = k, y = distribution, color=Process_type )) +
#  theme_bw() +
#  geom_point() +
#  geom_line()
#p


#library(akima)
#x= seq(1,10,0.1)
#y = seq(0.1,0.5,0.05)
#grid= expand.grid(x, y)
#f = function(x,y) Prior_on_K_SB(x, y, 112,112, runs=100)$Ek 
#z= mapply( function(x,y) Prior_on_K_SB(x, y, 112,112, runs=100)$Ek , grid$Var1, grid$Var2)
#spline <- interp(x,y,z,linear=FALSE)
#library(rgl)
#open3d(scale=c(1/diff(range(x)),1/diff(range(y)),1/diff(range(z))))
#with(spline,surface3d(x,y,z,alpha=.2))
#points3d(x,y,z)
#title3d(xlab="rating",ylab="complaints",zlab="privileges")
#axes3d()


#library(MBA)
#spline <- mba.surf(data.frame(x,y,z),100,100)


## using rootSolve and multiroot package!
#compute_fixed_parameters_SBPY_2d<- function(K,V,n, Ntrunc, n_sim=1000){
#  model<- function(x, K,V,n, Ntrunc){
#    F1<- Prior_on_K_SB(x[1], x[2], n,Ntrunc, runs=n_sim)$Ek - K
#    F2<- Prior_on_K_SB(x[1], x[2], n,Ntrunc, runs=n_sim)$Vk -V
#    c(F1 = F1, F2 = F2)
#  }
#  roots_values <- multiroot(f = function(x) model(x,K,V,n, Ntrunc), rtol = 1e-6, start=c(1,0.2))
#  return(list(alpha = roots_values$root[1],sigma=roots_values$root[2]))
#}
#compute_fixed_parameters_SBPY_2d(16, 20, 100,50,  100)










################################################################################  Test
## Test DP
#ratio<- DP_par$ratio
#nu2<-DP_par$nu2
#nu1<- DP_par$nu1

#nu2<-DP_par$nu2
#nu1<- DP_par$nu1
#nu2<-DP_par$nu2
#nu1<- DP_par$nu1
#aa<- rgamma(100000, nu1,nu2)

#x<- sapply(aa, functionDPM,n=112, N=112)
#mean(x)
#plot(density(x))

## Test DPM

## Test PY fixed
#functionPY(PY_fixed, 111,sigma_py=0.25)
#a = .bisec(f=function(x) functionDPM(x,111,111)- 16, a=0.01, b=10)


#### Test for Expectation and Variance
#values<- compute_fixed_parameters_PY_2d(16,20,111)


#####################################################################################################
