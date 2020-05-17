// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
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
}




double log_v_pdf_C(const Rcpp:: NumericVector& x, const Rcpp::List& pars){
  double alpha = pars["alpha"];
  double sigma = pars["sigma"];
  int H = pars["H"];
  int k = pars["k"];
  arma::vec nk_vec = as<arma::vec>(pars["nk_vec"]);
  double sum = 0;
  double val0;
  for (int j=0; j<k;++j){
    double val = 0;
    for (int l=1; l<(nk_vec[j]+1); ++l){
      val0 = pow(x[0]/H, l) * reccursive_Cnk(nk_vec[j], l, sigma);
      val += val0;
    }
    sum  += log(val);
  }
  return ( -x[0] + (alpha/sigma -1)*log(x[0]) + sum);
}

double log_v_pdf_comp_mat(const Rcpp:: NumericVector& x, const Rcpp::List& pars){
  double alpha = pars["alpha"];
  double sigma = pars["sigma"];
  int k = pars["k"];
  Rcpp::NumericVector nk_vec = pars["nk_vec"];
  arma::mat Cnk_mat =  as<arma::mat>(pars["Cnk_mat"]);
  double sum = 0;
  double val0;
  for (int j=0; j<k;++j){
   // Rcpp::Rcout << "nk matrix is" << std::endl << nk_vec[j] << std::endl;
    double val = 0;
    for (int l=1; l<(nk_vec[j]+1); ++l){
      val0 = exp(l*log(x[0]) + Cnk_mat(nk_vec[j]-1,l-1));
      //val0 = pow(x[0], l) * exp(Cnk_mat(nk_vec[j]-1,l-1));
      val =val+ val0;
    }
    sum  = sum +log(val);
  }
  return ( -x[0] + (alpha/sigma -1)*log(x[0]) + sum);
}

// A function to create external pointers for any of the functions above.  
// See http://gallery.rcpp.org/articles/passing-cpp-function-pointers/  
// If you write a new function above called new_name then add the following
//
// else if (fstr == "new_name")  
//   return(Rcpp::XPtr<funcPtr>(new funcPtr(&new_name))) ;  


// [[Rcpp::export]]  
SEXP create_xptr(std::string fstr) {  
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,  
                  const Rcpp::List& pars) ;  
  if (fstr == "log_v_pdf_C")  
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&log_v_pdf_C))) ;  
  else if (fstr == "log_v_pdf_comp_mat")  
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&log_v_pdf_comp_mat))) ; 
  else  
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ;  
}  

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
