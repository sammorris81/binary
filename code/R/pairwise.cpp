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
                   double rho, arma::mat d, double max_dist, int threads) {

#ifdef _OPENMP
  if(threads > 0) {
    omp_set_num_threads(threads);
  }
#endif

  double ll = 0; double joint; uword j; int ns = y.n_elem;

#pragma omp parallel for reduction(+ : ll)  \
  private(j, joint) shared(kernel, z, y, alpha, d) schedule(dynamic, 10)
    // did try static scheduling, and dynamic appears to be a little faster
    for (uword i = 0; i < (ns - 1); i++) {
      for (j = (i + 1); j < ns; j++) {
        // if (d(i, j) < max_dist) {
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
        // }
      }
    }

    return ll;
}

// [[Rcpp::export]]
double pairwiseCPPij(arma::mat kernel, double alpha, double z1, double z2,
                     int y1, int y2) {
  /* This is the likelihood evaluated for a single pair. Will be using this
     to evaluate the hessian and gradient in R */

  double ll; double joint;

  joint = -getJointPtr(&kernel, 0, 1, alpha);
  if (y1 == 0 && y2 == 0) {
    ll = joint;
  } else if (y1 == 1 && y2 == 0) {
    ll = log(exp(-1 / z2) - exp(joint));
  } else if (y1 == 0 && y2 == 1) {
    ll = log(exp(-1 / z1) - exp(joint));
  } else if (y1 == 1 && y2 == 1) {
    ll = log(1 - exp(-1 / z1) - exp(-1 / z2) + exp(joint));
  }

  return ll;
}