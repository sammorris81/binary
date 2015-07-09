// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat getKernelCPPwithID(SEXP wz_star, arma::mat a, Rcpp::List IDs,
                       uword nt, uword nknots, uword ns) {
  /* Function is way slower than just doing all the multiplication. May need *
   * to explore a more efficient way to access the elements in the IDs list. */

  // Hack for Rcpp to get arma::cube
  Rcpp::NumericVector wz_star_v(wz_star);
  Rcpp::IntegerVector arrayDims = wz_star_v.attr("dim");
  arma::cube wz_star_c(wz_star_v.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

  // temporary and return storage
  arma::mat kernel(ns, nt);
  arma::colvec a_t(nknots);
  arma::rowvec wz_star_it(nknots);
  arma::uvec these;
  uword this_a;

  for (uword t = 0; t < nt; t++) {
    a_t = a.col(t);
    for (uword i = 0; i < ns; i++) {
      // Hack to extract these from list
      SEXP temp = IDs[i];
      Rcpp::IntegerVector these_r(temp);
      these = Rcpp::as<arma::uvec>(these_r) - 1;  // vectors in R start with 1

      wz_star_it = wz_star_c.slice(t).row(i);  // select slice t and row i
      kernel(i, t) = 0;
      for (uword j = 0; j < these.n_elem; j++) {
        this_a = these[j];
        kernel(i, t) += wz_star_it(this_a) * a_t(this_a);
      }
    }
  }

  return kernel;
}

// [[Rcpp::export]]
arma::mat getKernelCPP(SEXP wz_star, arma::mat a) {

  // Hack for Rcpp to get arma::cube
  Rcpp::NumericVector wz_star_v(wz_star);
  Rcpp::IntegerVector arrayDims = wz_star_v.attr("dim");
  arma::cube wz_star_c(wz_star_v.begin(), arrayDims[0], arrayDims[1],
                       arrayDims[2], false);

  uword ns = arrayDims[0];
  // uword nknots = arrayDims[1];
  uword nt = arrayDims[2];

  // return storage
  arma::mat kernel(ns, nt);

  for (uword t = 0; t < nt; t++) {
    kernel.col(t) = wz_star_c.slice(t) * a.col(t);
  }

  return kernel;
}

// [[Rcpp::export]]
arma::cube getwzStarCPP(arma::mat z, arma::mat w, double alpha) {
  uword nknots = w.n_cols;
  uword ns = w.n_rows;
  uword nt = z.n_cols;

  // return storage
  arma::cube wz_star(ns, nknots, nt);

  // temp
  double wz_star_ilt;

  for (uword t = 0; t < nt; t++) {
    for (uword l = 0; l < nknots; l++) {
      for (uword i = 0; i < ns; i++) {
        if (w(i, l) > 0) {  // only include sites that are close to the knot
          wz_star_ilt = exp((log(w(i, l)) - log(z(i, t))) / alpha);
          wz_star(i, l, t) = wz_star_ilt;
          // wz_star(i, l, t) = wz_star_ilt < 1e-6 ? 0 : wz_star_ilt;
          // Rcout << wz_star(i, l, t) << std::endl;
          // Rcout << wz_star_ilt << std::endl;
        } else {
          wz_star(i, l, t) = 0;
        }
      }
    }
  }

  return wz_star;
}