// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat powCPP(arma::mat x, double alpha) {
  uword ncols = x.n_cols;
  uword nrows = x.n_rows;
  uword i; uword j;
  
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      x(i, j) = exp(alpha * log(x(i, j)));
    }
  }
  
  return(x);
}