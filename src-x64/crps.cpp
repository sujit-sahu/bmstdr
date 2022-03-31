#include <Rcpp.h>
using namespace Rcpp;
double crpscpp(NumericMatrix tmp);
double crps_one(NumericVector p);

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
NumericMatrix gibbs_cpp(int N, int thin) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rgamma(1, 3, 1 / (y * y + 4))[0];
      y = rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  
  return(mat);
}
// [[Rcpp::export]]
double crps_one(NumericVector p){
  // p is just one row of tmp 
  int N = p.size() - 1; 
  double result, dsum=0, dout=0; 
 
  for(int i=0; i<N; i++){
    dout = dout + fabs(p(i+1) -p(1));
    for(int j=(i+1); j<N; j++){
      dsum =  dsum + fabs(p(i+1) - p(j+1));
    }
  }
  result = dout/N - dsum/(N*N); 
  return(result);
}
// for one observation crps

// [[Rcpp::export]]
double crpscpp(NumericMatrix tmp){
    // tmp is r by (N+1)
    // first column is observation 
    // Next N columns are MCMC samples 
    //  tmp <- cbind(xval,xits)
    //  tmp <- na.omit(tmp)
    //  its <- length(tmp[1,])-1
    
    int N = tmp.ncol(), r = tmp.nrow(); 
    NumericVector onerow(N+1), resvec(r);
   
    
    for (int j =0; j <r; j++) {
      for(int i=0; i<N; i++){
        onerow(i) = tmp(j, i); 
      }
      resvec(j) = crps_one(onerow);
    }
    return (mean(resvec));
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//