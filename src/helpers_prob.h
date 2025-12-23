// helpers_prob.h
#ifndef HELPERS_PROB_H
#define HELPERS_PROB_H

#include <RcppArmadillo.h>

inline arma::vec rmvnorm_vec(const arma::vec& mean, const arma::mat& covar) {
  arma::vec out = mean;
  const arma::vec rand = Rcpp::rnorm(covar.n_cols, 0, 1);
  out += covar * rand;
  return out;
}

inline arma::mat irgamma_mat(const double shape, const arma::mat& scale) {
  const arma::uword nr = scale.n_rows;
  const arma::uword nc = scale.n_cols;
  arma::mat x(nr, nc, arma::fill::none);
  for (arma::uword r = 0; r < nr; r++) {
    for (arma::uword c = 0; c < nc; c++) {
      x(r, c) = 1.0 / R::rgamma(shape, scale(r, c));
    }
  }
  return x;
}

inline arma::mat irtgamma_mat(const double shape, const arma::mat& scale, const arma::mat& thres) {
  const arma::uword nr = scale.n_rows;
  const arma::uword nc = scale.n_cols;
  arma::mat x(nr, nc, arma::fill::none);
  for (arma::uword r = 0; r < nr; r++) {
    for (arma::uword c = 0; c < nc; c++) {
      const double max = R::pgamma(thres(r, c), shape, scale(r, c), true, false);
      double u = R::runif(0, max);
      x(r, c) = 1.0 / R::qgamma(u, shape, scale(r, c), true, false);
    }
  }
  return x;
}

inline arma::mat rtnorm_mat(const arma::mat& mean, const arma::mat& sd, const arma::mat& thres) {
  const arma::uword nr = mean.n_rows;
  const arma::uword nc = mean.n_cols;
  arma::mat x(nr, nc, arma::fill::none);
  for (arma::uword r = 0; r < nr; ++r) {
    for (arma::uword c = 0; c < nc; ++c) {
      double max = R::pnorm(thres(r, c), mean(r, c), sd(r, c), true, false);
      if (max > 0) {
        double u = R::runif(0, max);
        x(r, c) = R::qnorm(u, mean(r, c), sd(r, c), true, false);
      }
    }
  }
  return x;
}

inline arma::mat rnorm_mat(const arma::mat& mean, const arma::mat& sd) {
  const arma::uword nr = mean.n_rows;
  const arma::uword nc = mean.n_cols;
  arma::mat x(nr, nc, arma::fill::none);
  for (arma::uword r = 0; r < nr; ++r) {
    for (arma::uword c = 0; c < nc; ++c) {
      x(r, c) = R::rnorm(mean(r, c), sd(r, c));
    }
  }
  return x;
}


#endif // HELPERS_PROB_H
