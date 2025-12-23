#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <cmath>
#include "helpers_indexing.h"
#include "helpers_car.h"
#include "helpers_prob.h"
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;
using std::string;

double get_adj_xprod(const vec& Zkt, const field<uvec>& adj) {
  const uword n_region = adj.n_elem;
  double adj_xprod = 0;
  for (uword reg = 0; reg < n_region; ++reg) {
    adj_xprod += Zkt[reg] * arma::sum(Zkt.elem(adj[reg]));
  }
  return adj_xprod;
}

mat get_scale_sig(const cube& Z, const field<uvec>& adjacency, const vec& n_adj,
                  const double sig_b) {
  const uword n_group = Z.n_cols;
  const uword n_time = Z.n_slices;
  mat scale_sig(n_group, n_time, arma::fill::none);
  for (uword grp = 0; grp < n_group; ++grp) {
    for (uword time = 0; time < n_time; ++time) {
      const vec& Zkt = get_row(Z, grp, time);
      const double adj_xprod = get_adj_xprod(Zkt, adjacency);
      const double dot_Zkt2 = arma::dot(arma::square(Zkt), n_adj);
      scale_sig(grp, time) = 1.0 / (0.5 * (dot_Zkt2 - adj_xprod) + sig_b);
    }
  }
  return scale_sig;
}

mat get_thres_sig(const cube& beta, const uvec& n_isl_region, 
  const uword n_region, const string method,
  const mat& A, const mat& tau2, const double m0) {
  const uword n_group = beta.n_cols;
  const uword n_time = beta.n_slices;
  mat thres_sig(n_group, n_time, arma::fill::none);
  if (method == "binomial") {
    const mat pi = get_pi_rcar(beta, n_isl_region, n_region);
    thres_sig = (1.0 / ((A + pi) % (1 - pi)) - tau2 * (1 + 1.0 / m0)) * m0;
  } else if (method == "poisson") {
    thres_sig = (arma::log(1.0 / A + 1) - tau2 * (1 + 1.0 / m0)) * m0;
  }
  thres_sig = 1.0 / arma::clamp(thres_sig, 0.0, arma::datum::inf);
  return thres_sig;
}

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
  const uword n_island = isl_region.n_elem;

  const double shape_sig = (n_region - n_island) / 2 + sig_a;
  const mat scale_sig = get_scale_sig(Z, adjacency, n_adj, sig_b);
  sig2 = irgamma_mat(shape_sig, scale_sig);
  
  sample["sig2"] = sig2;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_sig2_rcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat sig2 = sample["sig2"];
  const cube& Z = sample["Z"];
  const cube& beta = sample["beta"];
  const mat& tau2 = sample["tau2"];
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec>& adjacency = sp_data["adjacency"];
  const vec& n_adj = sp_data["n_adj"];
  const field<uvec>& isl_region = sp_data["isl_region"];
  const uvec& n_isl_region = sp_data["n_isl_region"];
  const List& params = RSTr_obj["params"];
  const mat& A = params["A"];
  const double m0 = params["m0"];
  const string method = Rcpp::as<string>(params["method"]);
  const List& priors = RSTr_obj["priors"];
  const double sig_a = priors["sig_a"];
  const double sig_b = priors["sig_b"];
  const uword n_region = Z.n_rows;
  const uword n_island = isl_region.n_elem;

  const double shape_sig = (n_region - n_island) / 2 + sig_a;
  const mat scale_sig = get_scale_sig(Z, adjacency, n_adj, sig_b);
  const mat thres_sig = get_thres_sig(beta, n_isl_region, n_region, method, A, tau2, m0);
  sig2 = irtgamma_mat(shape_sig, scale_sig, thres_sig);

  sample["sig2"] = sig2;
  RSTr_obj["sample"] = sample;
}
