#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
SEXP dPS_wrapper(SEXP a, SEXP alpha, SEXP mid_points, SEXP bin_width) {
  arma::mat a_c = Rcpp::as<arma::mat>(a);
  int ns = a_c.n_rows; int nt = a_c.n_cols;
  double alpha_c = Rcpp::as<double>(alpha);
  arma::vec mid_points_c = Rcpp::as<arma::vec>(mid_points);
  arma::vec bin_width_c = Rcpp::as<arma::vec>(bin_width);
  int nbins = mid_points_c.n_elem;
  double integral; double llst; double psi; double c; double logint;
  double a_cst;
  arma::mat ll(ns, nt);

  for (int s = 0; s < ns; s++) {
    for (int t = 0; t < nt; t++) {
      a_cst = a_c(s, t);
      llst = log(alpha_c) - log(1 - alpha_c) - (1 / (1 - alpha_c)) *
             log(a_cst);
      integral = 0;
      for (int i = 0; i < nbins; i++) {
        psi = PI * mid_points_c[i];
        c = pow((sin(alpha_c * psi) / sin(psi)), (1 / (1 - alpha_c)));
        c *= sin((1 - alpha_c) * psi) / sin(alpha_c * psi);
        logint = log(c) - c * (1 / pow(a_cst, (alpha_c / (1 - alpha_c))));
        integral += exp(logint) * bin_width_c[i];
      }
      ll(s, t) = llst + log(integral);
    }
  }

  return ll;
}