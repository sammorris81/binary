// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat ifelsematCPP(arma::mat x, double tol) {
  int ncols = x.n_cols;
  int nrows = x.n_rows;
  
  for (uword i = 0; i < nrows; i++) {
    for (uword j = 0; j < ncols; j++) {
      x(i, j) = fabs(x(i, j)) < tol ? 0 : x(i, j);
    }
  }
  
  return x;
}

// [[Rcpp::export]]
arma::vec ifelsevecCPP(arma::vec x, double tol) {
  int n = x.n_elem;
  
  for (uword i = 0; i < n; i++) {
    x[i] = fabs(x[i]) < tol ? 0 : x[i];
  }
  
  return x;
}

