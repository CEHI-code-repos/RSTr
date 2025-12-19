#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using arma::mat;
using arma::cube;
using arma::uword;
using arma::uvec;
using Rcpp::List;
using Rcpp::String;

mat get_tau_scale(const cube& lambda, const cube& beta_0, const cube& Z,
                  const double tau_b) {
  const cube square_resid = arma::square(lambda - beta_0 - Z);
  mat sum_sq_gt = arma::sum(square_resid, 0);
  mat tau_scale = 1.0 / (0.5 * sum_sq_gt + tau_b);
  if (lambda.n_slices == 1) tau_scale = tau_scale.t();
  return tau_scale;
}

mat get_tau_scale_mst(const cube& lambda, const cube& beta_0, const cube& Z,
  const double tau_b) {
  const cube square_resid = arma::square(lambda - beta_0 - Z);
  const mat sum_sq_grp = arma::sum(arma::sum(square_resid, 0), 2);
  const mat tau_scale = 1.0 / (0.5 * sum_sq_grp.t() + tau_b);
  return tau_scale;
}

mat irgamma_mat(const double shape, const mat& scale) {
  const uword nr = scale.n_rows;
  const uword nc = scale.n_cols;
  mat tau2(nr, nc, arma::fill::zeros);
  for (uword r = 0; r < nr; r++) {
    for (uword c = 0; c < nc; c++) {
      tau2(r, c) = 1.0 / R::rgamma(shape, scale(r, c));
    }
  }
  return tau2;
}

//[[Rcpp::export]]
void update_tau2_default(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat tau2 = sample["tau2"];
  const cube& lambda = sample["lambda"];
  const cube& beta = sample["beta"];
  const cube& Z = sample["Z"];
  const List& sp_data = RSTr_obj["sp_data"];
  const uvec& isl_id = sp_data["isl_id"];
  const List& priors = RSTr_obj["priors"];
  const double tau_a = priors["tau_a"];
  const double tau_b = priors["tau_b"];
  const uword n_region = Z.n_rows;
  const double tau_shape = n_region / 2.0 + tau_a;
  const mat tau_scale = get_tau_scale(lambda, get_regs(beta, isl_id), Z, tau_b);
  tau2 = irgamma_mat(tau_shape, tau_scale);
  sample["tau2"] = tau2;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_tau2_rcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat tau2 = sample["tau2"];
  const cube& lambda = sample["lambda"];
  const cube& beta = sample["beta"];
  const cube& Z = sample["Z"];
  const mat& sig2 = sample["sig2"];
  const List& sp_data = RSTr_obj["sp_data"];
  const uvec& n_isl_region = sp_data["n_isl_region"];
  const uvec& isl_id = sp_data["isl_id"];
  const List& params = RSTr_obj["params"];
  const mat& A = params["A"];
  const double m0 = params["m0"];
  const String& method = params["method"];
  const List& priors = RSTr_obj["priors"];
  const double tau_a = priors["tau_a"];
  const double tau_b = priors["tau_b"];
  const uword n_region = Z.n_rows;
  const uword n_group = Z.n_cols;
  const uword n_time = Z.n_slices;
  const double tau_shape = n_region / 2.0 + tau_a;
  const cube& square_resid = pow(lambda - get_regs(beta, isl_id) - Z, 2);
  for (uword grp = 0; grp < n_group; grp++) {
    for (uword time = 0; time < n_time; time++) {
      double tau_scale = 1.0 / (sum(get_row(square_resid, grp, time)) / 2.0 + tau_b);
      double max = 0;
      if (method == "binomial") {
        double pi = sum(get_row(beta, grp, time) % n_isl_region / n_region);
        pi = exp(pi) / (1 + exp(pi));
        max = (1.0 / ((A(grp, time) + pi) * (1 - pi)) - sig2(grp, time) / m0) / (1 + 1.0 / m0);
      } else if (method == "poisson") {
        max = (log(1.0 / A(grp, time) + 1) - sig2(grp, time) / m0) / (1 + 1.0 / m0);
      }
      max = (max < 0) ? 0 : max;
      double u = R::runif(0, R::pgamma(1.0 / max, tau_shape, tau_scale, true, false));
      tau2(grp, time) = 1.0 / R::qgamma(u, tau_shape, tau_scale, true, false);
    }
  }
  sample["tau2"] = tau2;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_tau2_mstcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat tau2 = sample["tau2"];
  const cube& lambda = sample["lambda"];
  const cube& beta = sample["beta"];
  const cube& Z = sample["Z"];
  const List& priors = RSTr_obj["priors"];
  const double tau_a = priors["tau_a"];
  const double tau_b = priors["tau_b"];
  const List& sp_data = RSTr_obj["sp_data"];
  const uvec& isl_id = sp_data["isl_id"];
  const uword n_region = Z.n_rows;
  const uword n_time = Z.n_slices;
  const double tau_shape = n_region * n_time / 2.0 + tau_a;
  const mat& tau_scale = get_tau_scale_mst(lambda, get_regs(beta, isl_id), Z, tau_b);
  tau2 = irgamma_mat(tau_shape, tau_scale);
  sample["tau2"] = tau2;
  RSTr_obj["sample"] = sample;
}
