#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

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
  const double& tau_a = priors["tau_a"];
  const double& tau_b = priors["tau_b"];
  const uword n_region = Z.n_rows;
  const double tau_shape = n_region / 2.0 + tau_a;
  const cube square_resid = arma::square(lambda - get_regs(beta, isl_id) - Z);
  mat sum_sq_gt = arma::sum(square_resid, 0);
  mat tau_scale = 1.0 / (0.5 * sum_sq_gt + tau_b);
  tau2 = irgamma_mat(tau_shape, tau_scale);
  sample["tau2"] = tau2;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_tau2_eucar(List& RSTr_obj) {
  Rcpp::List sample = RSTr_obj["sample"];
  mat tau2 = sample["tau2"];
  cube lambda = sample["lambda"];
  cube beta = sample["beta"];
  cube Z = sample["Z"];
  mat sig2 = sample["sig2"];
  Rcpp::List sp_data = RSTr_obj["sp_data"];
  uvec n_isl_region = sp_data["n_isl_region"];
  uvec isl_id = sp_data["isl_id"];
  Rcpp::List params = RSTr_obj["params"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  Rcpp::List priors = RSTr_obj["priors"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword n_region = Z.n_rows;
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  double tau_shape = n_region / 2.0 + tau_a;
  cube square_resid = pow(lambda - get_regs(beta, isl_id) - Z, 2);
  for (uword grp = 0; grp < n_group; grp++) {
    for (uword time = 0; time < n_time; time++) {
      double tau_scale = 1.0 / (sum(get_row(square_resid, grp, time)) / 2.0 + tau_b);
      double tau_thres = 0;
      if (method == "binomial") {
        double pi = sum(get_row(beta, grp, time) % n_isl_region / n_region);
        pi = exp(pi) / (1 + exp(pi));
        tau_thres = (1.0 / ((A(grp, time) + pi) * (1 - pi)) - sig2(grp, time) / m0) / (1 + 1.0 / m0);
      } else if (method == "poisson") {
        tau_thres = (log(1.0 / A(grp, time) + 1) - sig2(grp, time) / m0) / (1 + 1.0 / m0);
      }
      tau_thres = (tau_thres < 0) ? 0 : tau_thres;
      double u = R::runif(0, R::pgamma(1.0 / tau_thres, tau_shape, tau_scale, true, false));
      tau2(grp, time) = 1.0 / R::qgamma(u, tau_shape, tau_scale, true, false);
    }
  }
  sample["tau2"] = tau2;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_tau2_mstcar(List& RSTr_obj) {
  Rcpp::List sample = RSTr_obj["sample"];
  mat tau2 = sample["tau2"];
  cube lambda = sample["lambda"];
  cube beta = sample["beta"];
  cube Z = sample["Z"];
  Rcpp::List priors = RSTr_obj["priors"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  Rcpp::List sp_data = RSTr_obj["sp_data"];
  uvec isl_id = sp_data["isl_id"];
  uword n_region = Z.n_rows;
  uword n_group  = Z.n_cols;
  uword n_time   = Z.n_slices;
  double tau_shape = n_region * n_time / 2.0 + tau_a;
  cube square_resid = pow(lambda - get_regs(beta, isl_id) - Z, 2) / 2.0;
  for (uword grp = 0; grp < n_group; grp++) {
    double tau_scale = 1.0 / (accu(square_resid.col(grp)) + tau_b);
    tau2[grp] = 1.0 / R::rgamma(tau_shape, tau_scale);
  }
  sample["tau2"] = tau2;
  RSTr_obj["sample"] = sample;
}
