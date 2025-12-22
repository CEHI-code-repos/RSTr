#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <cmath>
#include "cpp_helpers.h"
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;
using std::string;

//[[Rcpp::export]]
void update_sig2_default(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat sig2 = sample["sig2"];
  const cube& Z = sample["Z"];
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec>& adjacency = sp_data["adjacency"];
  const vec& n_adj = sp_data["n_adj"];
  const field<uvec>& isl_region = sp_data["isl_region"];
  const List& priors = RSTr_obj["priors"];
  const double sig_a = priors["sig_a"];
  const double sig_b = priors["sig_b"];
  const uword n_region = Z.n_rows;
  const uword n_group = Z.n_cols;
  const uword n_time = Z.n_slices;
  const uword n_island = isl_region.n_elem;
  for (uword grp = 0; grp < n_group; grp++) {
    for (uword time = 0; time < n_time; time++) {
      double sum_adj = 0;
      for (uword reg = 0; reg < n_region; reg++) {
        sum_adj += Z(reg, grp, time) * arma::sum(get_subregs(Z, adjacency(reg), grp, time));
      }
      double sig_shape = (n_region - n_island) / 2 + sig_a;
      double sig_scale = 1.0 / ((arma::sum(arma::pow(get_row(Z, grp, time), 2) % n_adj) - sum_adj) / 2 + sig_b);
      sig2(grp, time) = 1.0 / R::rgamma(sig_shape, sig_scale);
    }
  }
  sample["sig2"] = sig2;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_sig2_rcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat sig2 = sample["sig2"];
  cube Z = sample["Z"];
  cube beta = sample["beta"];
  mat tau2 = sample["tau2"];
  List sp_data = RSTr_obj["sp_data"];
  field<uvec> adjacency = sp_data["adjacency"];
  vec n_adj = sp_data["n_adj"];
  field<uvec> isl_region = sp_data["isl_region"];
  uvec n_isl_region = sp_data["n_isl_region"];
  List params = RSTr_obj["params"];
  mat A = params["A"];
  double m0 = params["m0"];
  const string method = Rcpp::as<string>(params["method"]);
  List priors = RSTr_obj["priors"];
  double sig_a = priors["sig_a"];
  double sig_b = priors["sig_b"];
  uword n_region = Z.n_rows;
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  uword n_island = isl_region.n_elem;
  for (uword grp = 0; grp < n_group; grp++) {
    for (uword time = 0; time < n_time; time++) {
      double sum_adj = 0;
      for (uword reg = 0; reg < n_region; reg++) {
        sum_adj += Z(reg, grp, time) * arma::sum(get_subregs(Z, adjacency(reg), grp, time));
      }
      double sig_shape = (n_region - n_island) / 2 + sig_a;
      double sig_scale = 1.0 / ((arma::sum(arma::pow(get_row(Z, grp, time), 2) % n_adj) - sum_adj) / 2 + sig_b);
      double max = 0;
      if (method == "binomial") {
        double pi = arma::sum(get_row(beta, grp, time) % n_isl_region / n_region);
        pi = std::exp(pi) / (1 + std::exp(pi));
        max = (1.0 / ((A(grp, time) + pi) * (1 - pi)) - tau2(grp, time) * (1 + 1.0 / m0)) * m0;
      } else if (method == "poisson") {
        max = (log(1.0 / A(grp, time) + 1) - tau2(grp, time) * (1 + 1.0 / m0)) * m0;
      }
      max = (max < 0) ? 0 : max;
      double u = R::runif(0, R::pgamma(1.0 / max, sig_shape, sig_scale, true, false));
      sig2(grp, time) = 1.0 / R::qgamma(u, sig_shape, sig_scale, true, false);
    }
  }
  sample["sig2"] = sig2;
  RSTr_obj["sample"] = sample;
}
