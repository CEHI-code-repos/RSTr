#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_beta_default(List& RSTr_obj) {
  Rcpp::List sample = RSTr_obj["sample"];
  cube beta = sample["beta"];
  cube lambda = sample["lambda"];
  cube Z = sample["Z"];
  mat tau2 = sample["tau2"];
  Rcpp::List sp_data = RSTr_obj["sp_data"];
  field<uvec> isl_region = sp_data["isl_region"];
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  uword n_island = isl_region.n_elem;
  for (uword isl = 0; isl < n_island; isl++) {
    uword n_isl_region = isl_region[isl].n_elem;
    mat var_beta = sqrt(tau2 / n_isl_region);
    mat mean_beta = mean(get_regs(lambda, isl_region[isl]) - get_regs(Z, isl_region[isl]), 0);
    if (n_time == 1) mean_beta = mean_beta.t();
    beta.row(isl) = mat(n_group, n_time, arma::fill::randn) % var_beta + mean_beta;
  }
  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_beta_eucar(List& RSTr_obj) {
  Rcpp::List sample = RSTr_obj["sample"];
  cube beta = sample["beta"];
  cube lambda = sample["lambda"];
  cube Z = sample["Z"];
  mat tau2 = sample["tau2"];
  mat sig2 = sample["sig2"];
  Rcpp::List sp_data = RSTr_obj["sp_data"];
  field<uvec> isl_region = sp_data["isl_region"];
  Rcpp::List params = RSTr_obj["params"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  uword n_island = isl_region.n_elem;
  cube tmZ = lambda - Z;
  for (uword isl = 0; isl < n_island; isl++) {
    uword n_isl_region = isl_region(isl).n_elem;
    mat var_t = tau2 + (tau2 + sig2) / m0;
    for (uword grp = 0; grp < n_group; grp++) {
      for (uword time = 0; time < n_time; time++) {
        double sd_beta = sqrt(tau2(grp, time) / n_isl_region);
        double mean_beta = mean(get_subregs(tmZ, isl_region(isl), grp, time));
        if (method == "binomial") {
          double pi_beta = pow(A(grp, time) - 1, 2) + 4 * (A(grp, time) - 1 / var_t(grp, time));
          double beta_thres = ((1 - A(grp, time)) + sqrt(pi_beta)) / 2;
          beta_thres = log(beta_thres / (1 - beta_thres));
          beta_thres = (beta_thres < 0) ? 0 : beta_thres;
          double beta_max = R::pnorm(beta_thres, mean_beta, sd_beta, true, false);
          if (beta_max > 0) {
            double u = R::runif(0, beta_max);
            beta(isl, grp, time) = R::qnorm(u, mean_beta, sd_beta, true, false);
          }
        } else if (method == "poisson") {
          beta(isl, grp, time) = R::rnorm(mean_beta, sd_beta);
        }
      }
    }
    
  }
  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_beta_mstcar(List& RSTr_obj) {
  Rcpp::List sample = RSTr_obj["sample"];
  cube beta = sample["beta"];
  cube lambda = sample["lambda"];
  cube Z = sample["Z"];
  vec tau2 = sample["tau2"];
  Rcpp::List sp_data = RSTr_obj["sp_data"];
  field<uvec> isl_region = sp_data["isl_region"];
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  uword n_island = isl_region.n_elem;
  for (uword isl = 0; isl < n_island; isl++) {
    uword n_isl_region = isl_region[isl].n_elem;
    mat var_beta = repmat(sqrt(tau2 / n_isl_region), 1, n_time);
    mat mean_beta = mean(get_regs(lambda, isl_region[isl]) - get_regs(Z, isl_region[isl]), 0);
    if (n_time == 1) mean_beta = mean_beta.t();
    beta.row(isl) = mat(n_group, n_time, arma::fill::randn) % var_beta + mean_beta;
  }
  sample["beta"] = beta;
  RSTr_obj["sample"] = sample;
}
