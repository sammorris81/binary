// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

double getJointPtr(arma::mat *kernel, int i, int j, double alpha) {
  double joint = 0;
  int nknots = (*kernel).n_cols;

  for (uword k = 0; k < nknots; k++) {
    joint += pow((*kernel)(i, k) + (*kernel)(j, k), alpha);
  }

  return joint;
}

// [[Rcpp::export]]
double pairwiseCPP(arma::mat kernel, double alpha, arma::vec z, arma::vec y,
                  int threads) {

#ifdef _OPENMP
  if(threads > 0) {
    omp_set_num_threads(threads);
  }
#endif

  double ll = 0; double joint; uword j; int ns = y.n_elem;

#pragma omp parallel for reduction(+ : ll)  \
  private(j, joint) shared(kernel, z, y, alpha) schedule(dynamic, 10)
    // did try static scheduling, and dynamic appears to be a little faster
    for (uword i = 0; i < (ns - 1); i++) {
      for (j = (i + 1); j < ns; j++) {
        joint = -getJointPtr(&kernel, i, j, alpha);
        if (y[i] == 0 && y[j] == 0) {
          ll += joint;
        } else if (y[i] == 1 && y[j] == 0) {
          ll += log(exp(-1 / z[j]) - exp(joint));
        } else if (y[i] == 0 && y[j] == 1) {
          ll += log(exp(-1 / z[i]) - exp(joint));
        } else if (y[i] == 1 && y[j] == 1) {
          ll += log(1 - exp(-1 / z[i]) - exp(-1 / z[j]) + exp(joint));
        }
      }
    }

    return ll;
}
