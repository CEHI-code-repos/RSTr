#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;
using Rcpp::String;

mat get_mean_beta(const cube& lambda, const cube& Z, const uvec& isl_idx) {
  cube sub_diff = get_regs(lambda, isl_idx) - get_regs(Z, isl_idx);
  mat mean_mat = arma::mean(sub_diff, 0);
  if (sub_diff.n_slices == 1) mean_mat = mean_mat.t();
  return mean_mat;
}

mat rtnorm_mat(const mat& mean_mat, const mat& sd_mat, const mat& beta_thres) {
  const uword n_group = mean_mat.n_rows;
  const uword n_time = mean_mat.n_cols;
  mat x(n_group, n_time, arma::fill::zeros);
  for (uword grp = 0; grp < n_group; grp++) {
    for (uword time = 0; time < n_time; time++) {
      double beta_max = R::pnorm(
        beta_thres(grp, time), 
        mean_mat(grp, time), 
        sd_mat(grp, time), 
        true, 
        false
      );
      if (beta_max > 0) {
        double u = R::runif(0, beta_max);
        x(grp, time) = R::qnorm(
          u, 
          mean_mat(grp, time), 
          sd_mat(grp, time), 
          true, 
          false
        );
      }
    }
  }
  return x;
}

mat rnorm_mat(const mat& mean_mat, const mat& sd_mat) {
  const uword n_group = mean_mat.n_rows;
  const uword n_time = mean_mat.n_cols;
  return (mat(n_group, n_time, arma::fill::randn) % sd_mat + mean_mat);
}

//[[Rcpp::export]]
void update_beta_default(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube beta = sample["beta"];
  const cube& lambda = sample["lambda"];
  const cube& Z = sample["Z"];
  const mat& tau2 = sample["tau2"];
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec>& isl_region = sp_data["isl_region"];
  const uword n_island = isl_region.n_elem;
  for (uword isl = 0; isl < n_island; isl++) {
    const uvec& isl_idx = isl_region[isl];
    const uword n_isl_region = isl_idx.n_elem;
    const mat mean_beta = get_mean_beta(lambda, Z, isl_idx);
    const mat sd_beta = arma::sqrt(tau2 / n_isl_region);
    beta.row(isl) = rnorm_mat(mean_beta, sd_beta);
  }
  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_beta_eucar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube beta = sample["beta"];
  const cube& lambda = sample["lambda"];
  const cube& Z = sample["Z"];
  const mat& tau2 = sample["tau2"];
  const mat& sig2 = sample["sig2"];
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec>& isl_region = sp_data["isl_region"];
  const List& params = RSTr_obj["params"];
  const mat& A = params["A"];
  const double& m0 = params["m0"];
  const String method = params["method"];
  const uword n_island = isl_region.n_elem;
  const mat var_latent = tau2 + (tau2 + sig2) / m0;
  for (uword isl = 0; isl < n_island; isl++) {
    const uvec& isl_idx = isl_region[isl];
    const uword n_isl_region = isl_idx.n_elem;
    const mat sd_beta = arma::sqrt(tau2 / n_isl_region);
    const mat mean_beta = get_mean_beta(lambda, Z, isl_idx);
    if (method == "binomial") {
      const mat pi_beta = arma::square(A - 1) + 4 * (A - 1.0 / var_latent);
      mat beta_thres = ((1 - A) + arma::sqrt(pi_beta)) / 2;
      beta_thres = log(beta_thres / (1 - beta_thres));
      beta_thres = arma::clamp(beta_thres, 0.0, arma::datum::inf);
      beta.row(isl) = rtnorm_mat(mean_beta, sd_beta, beta_thres);
    } else if (method == "poisson") {
      beta.row(isl) = rnorm_mat(mean_beta, sd_beta);
    }
  }
  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_beta_mstcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube beta = sample["beta"];
  const cube& lambda = sample["lambda"];
  const cube& Z = sample["Z"];
  const vec& tau2 = sample["tau2"];
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec>& isl_region = sp_data["isl_region"];
  const uword n_time = Z.n_slices;
  const uword n_island = isl_region.n_elem;
  for (uword isl = 0; isl < n_island; isl++) {
    const uvec& isl_idx = isl_region[isl];
    const uword n_isl_region = isl_idx.n_elem;
    const mat sd_beta = repmat(arma::sqrt(tau2 / n_isl_region), 1, n_time);
    const mat mean_beta = get_mean_beta(lambda, Z, isl_idx);
    beta.row(isl) = rnorm_mat(mean_beta, sd_beta);
  }
  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}
