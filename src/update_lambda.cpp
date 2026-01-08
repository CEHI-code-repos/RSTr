#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <cmath>
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;
using std::string;

inline double get_r2(const double l_0, const double l_star, const double n_0,
                     const string method) {
  double r2 = 0;
  if (method == "binomial") {
    r2 = n_0 * (std::log1p(std::exp(l_star)) - std::log1p(std::exp(l_0)));
  } else if (method == "poisson") {
    r2 = n_0 * (std::exp(l_star) - std::exp(l_0));
  }
  return r2;
}

inline double get_r(const double l_0, const double l_star, const double b_0,
                    const double Z_0, const double t_0, const double Y_0, 
                    const double n_0, const string& method) {
  const double r1 = Y_0 * (l_star - l_0);
  const double r2 = get_r2(l_0, l_star, n_0, method);
  const double delta = l_star - l_0;
  const double r3 = delta * (l_star + l_0 - 2.0 * (b_0 + Z_0));
  const double r = std::exp(r1 - r2 - 1.0 / (2.0 * t_0) * r3);
  return r;
}

//[[Rcpp::export]]
void update_lambda(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube lambda = sample["lambda"];
  const cube& Z = sample["Z"];
  const cube& beta = sample["beta"];
  mat tau2 = sample["tau2"];
  const List& data = RSTr_obj["data"];
  const cube& Y = data["Y"];
  const cube& n = data["n"];
  List priors = RSTr_obj["priors"];
  cube lambda_acpt = priors["lambda_acpt"];
  const cube& lambda_sd = priors["lambda_sd"];
  const List& sp_data = RSTr_obj["sp_data"];
  const uvec& isl_id = sp_data["isl_id"];
  const List& params = RSTr_obj["params"];
  const string method = Rcpp::as<string>(params["method"]);
  const uword n_region = Z.n_rows;
  const uword n_group = Z.n_cols;
  const uword n_time = Z.n_slices;
  if (tau2.n_cols == 1) tau2 = arma::repmat(tau2, 1, n_time);
  for (uword time = 0; time < n_time; ++time) {
    for (uword reg = 0; reg < n_region; ++reg) {
      for (uword grp = 0; grp < n_group; ++grp) {
        const double l_0 = lambda(reg, grp, time);
        const double l_star = std::clamp(R::rnorm(l_0, lambda_sd(reg, grp, time)), -100.0, 15.0);
        const double b_0 = beta(isl_id[reg], grp, time);
        const double Z_0 = Z(reg, grp, time);
        const double t_0 = tau2(grp, time);
        const double Y_0 = Y(reg, grp, time);
        const double n_0 = n(reg, grp, time);
        const double r = get_r(l_0, l_star, b_0, Z_0, t_0, Y_0, n_0, method);
        if (r >= R::runif(0, 1)) {
          lambda(reg, grp, time) = l_star;
          ++lambda_acpt(reg, grp, time);
        }
      }
    }
  }
  priors["lambda_acpt"] = lambda_acpt;
  RSTr_obj["priors"] = priors;
  sample["lambda"] = lambda;
  RSTr_obj["sample"] = sample;
}
