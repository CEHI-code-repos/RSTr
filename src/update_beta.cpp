#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "helpers_indexing.h"
#include "helpers_prob.h"
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;
using std::string;

mat get_thres_beta(const mat& tau2, const mat& sig2, const double m0,
                   const mat& A) {
  const mat var_latent = tau2 + (tau2 + sig2) / m0;
  const mat pi_beta = arma::square(A - 1) + 4 * (A - 1.0 / var_latent);
  mat thres_beta = ((1 - A) + arma::sqrt(pi_beta)) / 2;
  thres_beta = log(thres_beta / (1 - thres_beta));
  thres_beta = arma::clamp(thres_beta, 0.0, arma::datum::inf);
  return thres_beta;
}

mat get_mean_beta(const cube& lambda, const cube& Z, const uvec& isl_idx) {
  const cube sub_diff = get_regs(lambda, isl_idx) - get_regs(Z, isl_idx);
  mat mean_beta = arma::mean(sub_diff, 0);
  if (lambda.n_slices == 1) {
    return mean_beta.t();
  }
  return mean_beta;
}

mat get_sd_beta(const mat& tau2, const double n_isl_region, const uword n_time) {
  mat sd_beta = arma::sqrt(tau2 / n_isl_region);
  if (tau2.n_cols != n_time) {
    return repmat(sd_beta, 1, n_time);
  }
  return sd_beta;
}

//[[Rcpp::export]]
void update_beta_default(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube beta = Rcpp::as<cube>(sample["beta"]);
  const cube lambda = Rcpp::as<cube>(sample["lambda"]);
  const cube Z = Rcpp::as<cube>(sample["Z"]);
  const mat tau2 = Rcpp::as<mat>(sample["tau2"]);
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec> isl_region = Rcpp::as<field<uvec>>(sp_data["isl_region"]);
  const uword n_island = isl_region.n_elem;
  const uword n_time = Z.n_slices;

  for (uword isl = 0; isl < n_island; ++isl) {
    const uvec& isl_idx = isl_region[isl];
    const uword n_isl_region = isl_idx.n_elem;
    const mat mean_beta = get_mean_beta(lambda, Z, isl_idx);
    const mat sd_beta = get_sd_beta(tau2, n_isl_region, n_time);
    beta.row(isl) = rnorm_mat(mean_beta, sd_beta);
  }

  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_beta_rcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  auto beta = Rcpp::as<cube>(sample["beta"]);
  const auto lambda = Rcpp::as<cube>(sample["lambda"]);
  const auto Z = Rcpp::as<cube>(sample["Z"]);
  const auto tau2 = Rcpp::as<mat>(sample["tau2"]);
  const auto sig2 = Rcpp::as<mat>(sample["sig2"]);
  const List& sp_data = RSTr_obj["sp_data"];
  const auto isl_region = Rcpp::as<field<uvec>>(sp_data["isl_region"]);
  const List& params = RSTr_obj["params"];
  const uword n_island = isl_region.n_elem;
  const uword n_time = Z.n_slices;
  const auto A = Rcpp::as<mat>(params["A"]);
  const double m0 = params["m0"];
  const auto method = Rcpp::as<string>(params["method"]);

  const mat thres_beta = get_thres_beta(tau2, sig2, m0, A);
  for (uword isl = 0; isl < n_island; ++isl) {
    const uvec& isl_idx = isl_region[isl];
    const uword n_isl_region = isl_idx.n_elem;
    const mat sd_beta = get_sd_beta(tau2, n_isl_region, n_time);
    const mat mean_beta = get_mean_beta(lambda, Z, isl_idx);
    if (method == "binomial") {
      beta.row(isl) = rtnorm_mat(mean_beta, sd_beta, thres_beta);
    } else if (method == "poisson") {
      beta.row(isl) = rnorm_mat(mean_beta, sd_beta);
    }
  }

  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}
