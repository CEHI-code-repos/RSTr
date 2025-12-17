#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_lambda(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube lambda = sample["lambda"];
  cube Z = sample["Z"];
  cube beta = sample["beta"];
  mat tau2 = sample["tau2"];
  List data = RSTr_obj["data"];
  cube Y = data["Y"];
  cube n = data["n"];
  List priors = RSTr_obj["priors"];
  cube lambda_sd = priors["lambda_sd"];
  cube lambda_acpt = priors["lambda_acpt"];
  List sp_data = RSTr_obj["sp_data"];
  uvec isl_id = sp_data["isl_id"];
  List params = RSTr_obj["params"];
  String method = params["method"];
  uword n_region = Z.n_rows;
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  for (uword time = 0; time < n_time; time++) {
    for (uword reg = 0; reg < n_region; reg++) {
      for (uword grp = 0; grp < n_group; grp++) {
        double lambda_star = R::rnorm(lambda(reg, grp, time), lambda_sd(reg, grp, time));
        double rk1 = Y(reg, grp, time) * (lambda_star - lambda(reg, grp, time));
        double rk2 = 0;
        if (method == "binomial") {
          rk2 = n(reg, grp, time) * (log(1 + exp(lambda_star)) - log(1 + exp(lambda(reg, grp, time))));
        } 
        if (method == "poisson") {
          rk2 = n(reg, grp, time) * (exp(lambda_star) - exp(lambda(reg, grp, time)));
        }
        double rk3a = pow(lambda_star            - beta(isl_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk3b = pow(lambda(reg, grp, time) - beta(isl_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk = exp(rk1 - rk2 - 1 / (2 * tau2[grp]) * (rk3a - rk3b));
        if (rk >= R::runif(0, 1)) {
          lambda(reg, grp, time) = lambda_star;
          lambda_acpt(reg, grp, time) ++;
        }
      }
    }
  }
  priors["lambda_acpt"] = lambda_acpt;
  RSTr_obj["priors"] = priors;
  sample["lambda"] = lambda;
  RSTr_obj["sample"] = sample;
}
