#include <RcppArmadillo.h>
#include <RcppDist.h>

using namespace Rcpp;
using namespace arma;

inline double get_r2(const double l_0, const double l_star, const double n_0,
                     const std::string& method) {
  if (method == "binomial") {
    return n_0 * (std::log1p(std::exp(l_star)) - std::log1p(std::exp(l_0)));
  } else {
    return n_0 * (std::exp(l_star) - std::exp(l_0));
  }
}

inline double get_r(const double l_0, const double l_star, const double b_0,
                    const double Z_0, const double t_0, const double Y_0, 
                    const double n_0, const std::string& method) {
  const double r1 = Y_0 * (l_star - l_0);
  const double r2 = get_r2(l_0, l_star, n_0, method);
  const double delta = l_star - l_0;
  const double r3 = delta * (l_star + l_0 - 2.0 * (b_0 + Z_0));
  return std::exp(r1 - r2 - 1.0 / (2.0 * t_0) * r3);
}

//[[Rcpp::export]]
void update_lambda(List& RSTr_obj) {
  // 1. Access R objects directly (no copying yet)
  List sample = RSTr_obj["sample"];
  List data = RSTr_obj["data"];
  List priors = RSTr_obj["priors"];
  List sp_data = RSTr_obj["sp_data"];
  List params = RSTr_obj["params"];

  NumericVector r_lambda = sample["lambda"];
  IntegerVector dims = r_lambda.attr("dim");
  cube lambda(r_lambda.begin(), dims[0], dims[1], dims[2], false, true);

  NumericVector r_Z = sample["Z"];
  const cube Z(r_Z.begin(), dims[0], dims[1], dims[2], false, true);

  NumericVector r_beta = sample["beta"];
  IntegerVector b_dims = r_beta.attr("dim");
  const cube beta(r_beta.begin(), b_dims[0], b_dims[1], b_dims[2], false, true);

  const mat tau2 = as<mat>(sample["tau2"]);

  NumericVector r_Y = data["Y"];
  IntegerVector y_dims = r_Y.attr("dim");
  const cube Y(r_Y.begin(), y_dims[0], y_dims[1], y_dims[2], false, true);
  NumericVector r_n = data["n"];
  IntegerVector n_dims = r_n.attr("dim");
  const cube n(r_n.begin(), n_dims[0], n_dims[1], n_dims[2], false, true);
  
  NumericVector r_acpt = priors["lambda_acpt"];
  cube lambda_acpt(r_acpt.begin(), dims[0], dims[1], dims[2], false, true);
  
  const cube lambda_sd = as<cube>(priors["lambda_sd"]);
  const uvec isl_id = as<uvec>(sp_data["isl_id"]);
  const std::string method = as<std::string>(params["method"]);

  const uword n_region = dims[0];
  const uword n_group  = dims[1];
  const uword n_time   = dims[2];

  bool tau_is_vec = (tau2.n_cols == 1);

  for (uword time = 0; time < n_time; ++time) {
    for (uword grp = 0; grp < n_group; ++grp) {
      const double t_0 = tau_is_vec ? tau2(grp, 0) : tau2(grp, time);
      
      for (uword reg = 0; reg < n_region; ++reg) {
        const double l_0 = lambda(reg, grp, time);
        const double l_star = std::clamp(R::rnorm(l_0, lambda_sd(reg, grp, time)), -100.0, 15.0);
        
        const double r = get_r(l_0, l_star, beta(isl_id[reg], grp, time), 
                               Z(reg, grp, time), t_0, Y(reg, grp, time), 
                               n(reg, grp, time), method);
        
        if (r >= R::runif(0, 1)) {
          lambda(reg, grp, time) = l_star;
          lambda_acpt(reg, grp, time) += 1.0;
        }
      }
    }
  }
}
