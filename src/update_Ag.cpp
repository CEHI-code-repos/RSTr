#include <RcppArmadillo.h>
#include <RcppDist.h>
using arma::mat;
using arma::cube;
using arma::uword;
using Rcpp::List;

//[[Rcpp::export]]
void update_Ag(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat Ag = sample["Ag"];
  const cube& G = sample["G"];
  const List& priors = RSTr_obj["priors"];
  const mat& Ag_scale = priors["Ag_scale"];
  const double G_df = priors["G_df"];
  const double Ag_df = priors["Ag_df"];
  const uword n_time = G.n_slices;
  const uword n_group = G.n_rows;

  mat Ag_covar(n_group, n_group, arma::fill::zeros);
  Ag_covar += arma::inv_sympd(Ag_scale);
  for (uword time = 0; time < n_time; ++time) {
    Ag_covar += arma::inv_sympd(G.slice(time));
  }
  Ag = rwish(n_time * G_df + Ag_df, arma::inv_sympd(Ag_covar));
  
  sample["Ag"] = Ag;
  RSTr_obj["sample"] = sample;
}
