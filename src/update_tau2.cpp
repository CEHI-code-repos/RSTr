#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using arma::mat;
using arma::cube;
using arma::uword;
using arma::uvec;
using Rcpp::List;
using std::string;

mat get_tau_scale(const cube& lambda, const cube& beta_0, const cube& Z,
                  const double tau_b) {
  const cube square_resid = arma::square(lambda - beta_0 - Z);
  mat sum_sq_gt = arma::sum(square_resid, 0);
  mat scale = 1.0 / (0.5 * sum_sq_gt + tau_b);
  if (lambda.n_slices == 1) {
    return scale.t();
  }
  return scale;
}

mat get_tau_scale_mst(const cube& lambda, const cube& beta_0, const cube& Z,
                      const double tau_b) {
  const cube square_resid = arma::square(lambda - beta_0 - Z);
  const mat sum_sq_grp = arma::sum(arma::sum(square_resid, 0), 2);
  const mat scale = 1.0 / (0.5 * sum_sq_grp.t() + tau_b);
  return scale;
}

mat get_pi_tau(const cube& beta, const uvec& n_isl_reg, const uword n_reg) {
  const uword n_group = beta.n_cols;
  const uword n_time = beta.n_slices;
  mat pi(n_group, n_time);
  for (uword grp = 0; grp < n_group; ++grp) {
    for (uword time = 0; time < n_time; ++time) {
      pi(grp, time) = arma::sum(get_row(beta, grp, time) % n_isl_reg / n_reg);
    }
  }
  return (arma::exp(pi) / (1 + arma::exp(pi)));
}

mat get_tau_thres(const cube& beta, const uvec& n_isl_region, 
                  const uword n_region, const string method,
                  const mat& A, const mat& sig2, const double m0) {
  const uword n_group = beta.n_cols;
  const uword n_time = beta.n_slices;
  mat thres(n_group, n_time);
  if (method == "binomial") {
    const mat pi = get_pi_tau(beta, n_isl_region, n_region);
    thres = (1.0 / ((A + pi) % (1 - pi)) - sig2 / m0) / (1 + 1.0 / m0);
  } else if (method == "poisson") {
    thres = (log(1.0 / A + 1) - sig2 / m0) / (1 + 1.0 / m0);
  }
  thres = 1.0 / arma::clamp(thres, 0.0, arma::datum::inf);
  return thres;
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

  const double shape = 0.5 * n_region + tau_a;
  const mat scale = get_tau_scale(lambda, get_regs(beta, isl_id), Z, tau_b);
  tau2 = irgamma_mat(shape, scale);

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
  const string method = Rcpp::as<string>(params["method"]);
  const List& priors = RSTr_obj["priors"];
  const double tau_a = priors["tau_a"];
  const double tau_b = priors["tau_b"];
  const uword n_region = Z.n_rows;

  const double shape = 0.5 * n_region + tau_a;
  const mat scale = get_tau_scale(lambda, get_regs(beta, isl_id), Z, tau_b);
  const mat thres = get_tau_thres(beta, n_isl_region, n_region, method, A, sig2, m0);
  tau2 = irtgamma_mat(shape, scale, thres);

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

  const double shape = 0.5 * n_region * n_time + tau_a;
  const mat scale = get_tau_scale_mst(lambda, get_regs(beta, isl_id), Z, tau_b);
  tau2 = irgamma_mat(shape, scale);

  sample["tau2"] = tau2;
  RSTr_obj["sample"] = sample;
}
