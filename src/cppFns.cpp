#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List byRcpp(const int nr, 
                     const arma::mat frommat,
                     arma::mat totmat, 
                     arma::mat summat, 
                     arma::mat minmat, 
                     arma::mat maxmat){
  int i, j;
  double s;
  
  for(int k = 0; k < nr; k++){
    
    i = frommat(k,0) - 1;
    j = frommat(k,1) - 1;
    s = frommat(k,2);
    totmat(i,j) = totmat(i,j) + 1;
    summat(i,j) = summat(i,j) + s;
    
    if(s > maxmat(i,j))
      maxmat(i,j) = s;   
    
    if(s < minmat(i,j))
      minmat(i,j) = s;
  }
  
  return Rcpp::List::create(Rcpp::Named("total")=totmat,
                            Rcpp::Named("sum")=summat,
                            Rcpp::Named("min")=minmat,
                            Rcpp::Named("max")=maxmat);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List condMVNRcpp(const arma::uvec cdex, 
                                 const arma::uvec gdex, 
                                 const arma::mat xx, arma::mat mu, 
                                 const arma::mat sigma) {
  
  arma::mat sinv = arma::inv_sympd(sigma.submat(gdex,gdex));
  arma::mat p1 = sigma.submat(cdex, gdex) * sinv;
  arma::mat mu1 = mu.cols(cdex) + trans(p1 * trans(xx.cols(gdex) - mu.cols(gdex)));
  arma::mat vr1 = sigma.submat(cdex, cdex) - p1 * sigma.submat(gdex,cdex);
  
  return Rcpp::List::create(Rcpp::Named("mu")=mu1,
                            Rcpp::Named("vr")=vr1);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double tnormRcpp(double lo, double hi, double mu, double sig){
  
  double q1, q2, z;
  
  q1 = Rf_pnorm5(lo,mu,sig,1,0);
  q2 = Rf_pnorm5(hi,mu,sig,1,0);
  z = q1 + unif_rand()*(q2-q1);
  z = Rf_qnorm5(z, mu, sig, 1, 0);
  
  if(z > hi){
    z = lo;
  }
  
  if(z < lo){
    z = hi;
  }
  return(z);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat trMVNmatrixRcpp(arma::mat avec, arma::mat muvec, 
                                arma::mat smat, arma::mat lo,
                                arma::mat hi, arma::uvec whichSample, 
                                arma::uvec idxALL){
  int cindex;
  arma::rowvec av;
  arma::rowvec mv;
  arma::vec mAs(2);
  int nm = smat.n_rows;
  int nr = muvec.n_rows;
  arma::rowvec p1(nm-1);
  arma::mat sin(nm-1, nm-1);
  arma::uvec cid(1);
  arma::uvec idx;
  arma::mat m1(1,1);
  arma::mat s1(1,1);
  double tiny = min(smat.diag())*.0001;
  int nk = whichSample.n_elem;
  
  arma::mat A(nr, nm); A.fill(NA_REAL);
  arma::umat idxALLm(nm-1, nm);
  
  for(int j=0; j < nm; j++)
    
    idxALLm.col(j) = idxALL.elem( find(idxALL != j) );
  
  for(int i = 0; i < nr ; i++){
    
    for(int k = 0; k < nk; k++){
      
      cindex = whichSample[k] - 1;
      
      av = avec.row(i);
      mv = muvec.row(i);
      
      cid(0) = cindex;
      idx = idxALLm.col(cindex);
      sin = arma::inv_sympd(smat.submat(idx, idx));
      p1 = trans(smat.submat(idx, cid)) * sin;
      
      m1 = mv[cindex] + dot(p1, (av.elem(idx) - mv.elem(idx)));
      s1 = smat(cindex,cindex) - dot(p1, smat.submat(cid, idx)) ;
      
      mAs[0] = m1(0,0);
      mAs[1] = s1(0,0);
      if(mAs[1] < 0) mAs[1] = tiny;  
      
      double sss = pow(mAs[1],.5);
      
      avec(i,cindex) = tnormRcpp(lo(i,cindex), hi(i,cindex), mAs[0], sss);
      A(i,cindex) = avec(i,cindex);
    
    }
  }
  return A;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rmvnormRcpp(int n, arma::vec mu, arma::mat sigma) {
  
  int ncols = sigma.n_cols;
  bool success = false;
  arma::mat S = sigma;
  int i = 1;
  
  arma::mat Y = randn(n, ncols);
  
  while(success == false && i < 5){
    
    success = chol(S, sigma);
    
    if(success == false){
      sigma += eye(ncols,ncols) * 1e-5;
    }
    
    i = i + 1;
  }
  
  if(success == false){
    //    throw std::range_error("sigma not positive definite");
    return arma::repmat(mu*0, 1, n).t();
  }
  
  return arma::repmat(mu, 1, n).t() + Y * chol(sigma);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat solveRcpp(arma::mat A) {
  arma::mat AA(A);
  arma::mat Ainv = arma::inv_sympd(AA);
  return Ainv;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat getPmatKRcpp(arma::vec pveck, arma::mat Yk, arma::mat Zk, arma::mat Xk,
                         arma::mat Bk, arma::mat Wk, double sigmasqk) {
  vec pvec = pveck;
  mat Y(Yk);
  mat Z(Zk);
  mat X(Xk);
  mat B(Bk);
  mat W(Wk);
  double sigmasq = sigmasqk;
  
  int N = Z.n_rows;
  int nn = Y.n_rows;
  int q = Y.n_cols;
  
  vec lpvec = log(pveck);
  double mxpval;
  mat pmat(q,N);
  mat emat(q,N);
  mat epmat(q,N);
  vec prdvec(nn);
  
  for(int i = 0; i < q; ++i){
    
    for(int j = 0; j < N; ++j){
      prdvec = Y.col(i) - X*B.row(i).t() - W*trans(Z.row(j));
      pmat(i,j) = lpvec(j) - (0.5/sigmasq)*dot(prdvec,prdvec);
    }
    mxpval = max(pmat.row(i));
    pmat.row(i) = pmat.row(i) - mxpval;
    
  }
  emat = exp(pmat);
  epmat = normalise( emat, 1, 1 );
  
  return epmat;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fnZRcpp(arma::vec kk, arma::mat Yk, arma::mat Xk, arma::mat Dk,
                      arma::mat Bk, arma::mat Wk, double sigmasqk, int Nz) {
  vec k = kk;
  mat Y(Yk);
  mat X(Xk);
  mat B(Bk);
  mat W(Wk);
  mat D(Dk);
  double sigmasq = sigmasqk;
  
  int N = Nz;
  int r = W.n_cols;
  int nn = Y.n_rows;
  //int q = Y.n_cols;
  int s = 0;
  
  vec kstar = unique(kk) - 1;
  vec knotstar(N - kstar.size());
  for(int j = 0; j < N; ++j){
    if(all(kstar != j)){
      knotstar[s] = j;
      s = s+1;
    }
  }
  
  mat Z(N,r);
  mat WtW = W.t()*W;
  mat Dinv = arma::inv_sympd(D);
  int nkk = kstar.size();
  
  uvec J;
  int js;
  mat CovZj(r,r);
  vec ssY(nn);
  vec meanZj(r);
  mat tempmat;
  vec tempvec(r);
  
  for(int i = 0; i < nkk; ++i){
    J = find(k == (kstar(i) + 1));
    js = J.size();
    CovZj = arma::inv_sympd(as_scalar(js/sigmasq)*WtW + Dinv);
    tempmat = Y.cols(J) - X*B.rows(J).t();
    if(J.size() > 1){
      ssY = arma::sum(tempmat,1);
    } else {
      ssY = tempmat;
    }
    
    tempvec = vectorise(as_scalar(1/sigmasq)*CovZj*W.t()*ssY);
    meanZj = tempvec;
    Z.row(kstar(i)) = rmvnormRcpp(1, meanZj, CovZj);
    //arma::repmat(meanZj, 1, 1).t() + randn(1, r) * chol(CovZj);
    J.reset();
  }
  
  int nss = knotstar.size();
  
  for(int i = 0; i < nss; ++i) {
    Z.row(knotstar(i)) = rmvnormRcpp(1, zeros<arma::vec>(r), D);
    //repmat(zeros<arma::vec>(r), 1, 1).t() + randn(1, r) * chol(D);
  }
  
  //rmvnormArma(knotstar.size(),zeros<vec>(r),D);
  
  return Z;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat invWbyRcpp(double sigsq, arma::mat A) {
  int s = A.n_rows; 
  int r = A.n_cols; 
  arma::mat Ds = eye<arma::mat>(s,s);
  arma::mat Dr = eye<arma::mat>(r,r);
  arma::mat Q = arma::inv_sympd(Dr + as_scalar(1/sigsq)*trans(A)*A);
  return arma::mat (as_scalar(1/sigsq)*(Ds - as_scalar(1/sigsq)*A*Q*trans(A)));
}

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnormRcpp(arma::mat x,  
                      arma::mat mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean.row(i));    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat randEffRcpp(arma::mat v, arma::vec ngvec, 
                      arma::mat sinv, arma::mat CImat){
  
  int nv = ngvec.n_elem;
  int s = v.n_rows;
  arma::mat Z(s,nv); Z.fill(0);
  
  for(int k = 0; k < nv; k++){
    arma::mat VI = as_scalar(ngvec(k))*sinv + CImat;
    arma::mat V = arma::inv_sympd(VI);
    arma::vec mu =  V * v.col(k);
    Z.col(k) = trans( rmvnormRcpp(1, mu, V) );
  }
  return(Z);
}
