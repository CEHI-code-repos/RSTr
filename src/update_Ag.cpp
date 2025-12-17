#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_Ag(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  mat Ag = sample["Ag"];
  cube G = sample["G"];
  List priors = RSTr_obj["priors"];
  mat Ag_scale = priors["Ag_scale"];
  double G_df = priors["G_df"];
  double Ag_df = priors["Ag_df"];
  uword n_time  = G.n_slices;
  uword n_group = G.n_rows;
  mat Ag_covar(n_group, n_group, fill::zeros);
  Ag_covar += inv(Ag_scale);
  for (uword time = 0; time < n_time; time++) {
    Ag_covar += inv(G.slice(time));
  }
  Ag = rwish(n_time * G_df + Ag_df, inv(Ag_covar));
  sample["Ag"] = Ag;
  RSTr_obj["sample"] = sample;
}
