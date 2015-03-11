library(inline)
library(microbenchmark)

# used when evaluating the postive stable density
ld <- function(u, a, alpha) {
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * log(a) +
          log(c) - c * (1 / a^(alpha / (1 - alpha)))

  return(exp(logd))
}

dPS.R <- function(a, alpha, mid.points, bin.width) {
  if (is.null(dim(a))) {
    ns <- length(a)
    nt <- 1
    a <- matrix(a, ns, nt)  # turn it into a matrix
  } else {
    ns <- nrow(a)
    nt <- ncol(a)
  }

  ll <- matrix(-Inf, ns, nt)
  for (t in 1:nt) {
    for (s in 1:ns) {
      ll[s, t] <- log(sum(bin.width * ld(mid.points, a[s, t], alpha)))
    }
  }
  return(ll)
}

code <- '
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

  return Rcpp::List::create(
    Rcpp::Named("ll") = ll
  );
'

dPS.wrapper <- cxxfunction(signature(a="numeric", alpha="numeric",
                                     mid_points="numeric", bin_width="numeric"),
                                     code, plugin="RcppArmadillo")


dPS.Rcpp <- function(a, alpha, mid.points, bin.width) {
  if (is.null(dim(a))) {
    ns <- length(a)
    nt <- 1
    a <- matrix(a, ns, nt)  # turn it into a matrix
  } else {
    ns <- nrow(a)
    nt <- ncol(a)
  }

  results <- dPS.wrapper(a=a, alpha=alpha, mid_points=mid.points,
                         bin_width=bin.width)
  return(results$ll)
}

npts <- 70
u.beta     <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
mid.points <- (u.beta[-1] + u.beta[-(npts + 1)]) / 2
bin.width  <- u.beta[-1] - u.beta[-(npts + 1)]
alpha      <- 0.5

a <- matrix(1:81, 81, 20)

res.R <- dPS.R(a, alpha, mid.points, bin.width)
res.Rcpp <- dPS.Rcpp(a, alpha, mid.points, bin.width)
sum( (res.R / res.Rcpp) == 1)

microbenchmark(dPS.R(a=a, alpha=alpha, mid.points=mid.points,
                     bin.width=bin.width),
               dPS.Rcpp(a=a, alpha=alpha, mid.points=mid.points,
                        bin.width=bin.width),
               times = 100)