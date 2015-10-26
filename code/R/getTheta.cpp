// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat getThetaCPP(SEXP wz, arma::mat a_star, double alpha) {
  /* a_star = a^alpha
   * from some quick benchmarks, it would appear that R is quicker with the
   * element-wise power operation. So, we're just passing in a_star.
   */
  
  // Hack for Rcpp to get arma::cube
  Rcpp::NumericVector wz_v(wz);
  Rcpp::IntegerVector arrayDims = wz_v.attr("dim");
  arma::cube wz_c(wz_v.begin(), arrayDims[0], arrayDims[1],
                  arrayDims[2], false);
  
  uword ns = arrayDims[0];
  uword nknots = arrayDims[1];
  uword nt = arrayDims[2];
  
  double alpha_inv = 1 / alpha;
  
  // temp storage
  arma::mat awz_lt(ns, nknots);
  
  // return storage
  arma::mat theta(ns, nt);
  
  for (uword t = 0; t < nt; t++) {
    // take the slice for the day and elementwise multiply a_star for each knot
    awz_lt = wz_c.slice(t);
    awz_lt.each_row() %= a_star.col(t).t();  // transpose to a row vector  
    awz_lt = pow(awz_lt, alpha_inv);
    theta.col(t) = sum(awz_lt, 1); // add up across the row
  }
  
  return theta;
}



// // [[Rcpp::export]]
// arma::cube getwzStarCPP(arma::mat z, arma::mat w, double alpha) {
//   uword nknots = w.n_cols;
//   uword ns = w.n_rows;
//   uword nt = z.n_cols;
// 
//   // return storage
//   arma::cube wz_star(ns, nknots, nt);
// 
//   // temp
//   double wz_star_ilt;
// 
//   for (uword t = 0; t < nt; t++) {
//     for (uword l = 0; l < nknots; l++) {
//       for (uword i = 0; i < ns; i++) {
//         // Rcout << w(i, l) << std::endl;
//         if (w(i, l) > 0) {  // only include sites that are close to the knot
//           wz_star_ilt = exp((log(w(i, l)) - log(z(i, t))) / alpha);
//           wz_star(i, l, t) = wz_star_ilt;
//           // wz_star(i, l, t) = wz_star_ilt < 1e-6 ? 0 : wz_star_ilt;
//           // Rcout << wz_star(i, l, t) << std::endl;
//           // Rcout << wz_star_ilt << std::endl;
//         } else {
//           wz_star(i, l, t) = 0;
//         }
//       }
//     }
//   }
// 
//   return wz_star;
// }

// [[Rcpp::export]]
arma::cube getwzCPP(arma::mat z, arma::mat w) {
  uword nknots = w.n_cols;
  uword ns = w.n_rows;
  uword nt = z.n_cols;
  
  // return storage
  arma::cube wz(ns, nknots, nt);
  
  // temp
  double wz_ilt;
  
  for (uword t = 0; t < nt; t++) {
    for (uword l = 0; l < nknots; l++) {
      for (uword i = 0; i < ns; i++) {
        // Rcout << w(i, l) << std::endl;
        if (w(i, l) > 0) {  // only include sites that are close to the knot
          wz_ilt = exp(log(w(i, l)) - log(z(i, t)));
          wz(i, l, t) = wz_ilt;
          // wz_star(i, l, t) = wz_star_ilt < 1e-6 ? 0 : wz_star_ilt;
          // Rcout << wz_star(i, l, t) << std::endl;
          // Rcout << wz_star_ilt << std::endl;
        } else {
          wz(i, l, t) = 0;
        }
      }
    }
  }
  
  return wz;
}

// [[Rcpp::export]]
arma::mat getawCPP(arma::mat a_star, arma::mat w, double alpha) {
  /* when beta is being updated, it is quicker to work with 
   * aw_it = sum(a_lt w_li^1 / alpha)
   * 
   * this avoids unnecessary calculation and summing since only z changes.
   * using a_star to help with numerical stability in multiplication.
   */ 
  uword ns = w.n_rows;
  uword nt = a_star.n_cols;
  
  // return storage
  arma::mat aw = zeros<mat>(ns, nt);
  
  // temp
  
  for (uword t = 0; t < nt; t++) {
    // Take the tth column of A and multiply each row of w, then add up the row
    aw.col(t) = sum(pow(w.each_row() % a_star.col(t).t(), 1 / alpha), 1);
  }
  
  return aw;
}

// [[Rcpp::export]]
arma::mat getawCPP2(arma::mat a, arma::mat w, double alpha) {
  /* when beta is being updated, it is quicker to work with 
   * aw_it = sum(a_lt w_li^1 / alpha)
   * 
   * this avoids unnecessary calculation and summing since only z changes.
   * using a_star to help with numerical stability in multiplication.
   */ 
//   uword ns = w.n_rows;
//   uword nt = a.n_cols;
  
  // return storage
  arma::mat aw = pow(w, 1 / alpha) * a;
  
  return aw;
}