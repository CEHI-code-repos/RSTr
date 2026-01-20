#ifndef RCPPEXTRACT_H
#define RCPPEXTRACT_H

#include <RcppArmadillo.h>

inline arma::cube get_cube(const Rcpp::List& x, const char* name) {
  Rcpp::NumericVector nv(x[name]);
  return Rcpp::as<arma::cube>(nv);
}

inline arma::mat get_mat(const Rcpp::List& x, const char* name) {
  Rcpp::NumericMatrix nm(x[name]);
  return Rcpp::as<arma::mat>(nm);
}

inline arma::uvec get_uvec(const Rcpp::List& x, const char* name) {
  Rcpp::IntegerVector iv(x[name]);
  return Rcpp::as<arma::uvec>(iv);
}

inline std::string get_string(const Rcpp::List& x, const char* name) {
  return Rcpp::as<std::string>(x[name]);
}

#endif
