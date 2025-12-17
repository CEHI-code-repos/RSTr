#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_rho(List& RSTr_obj) {
  Rcpp::List sample = RSTr_obj["sample"];
  vec rho = sample["rho"];
  cube G = sample["G"];
  cube Z = sample["Z"];
  Rcpp::List priors = RSTr_obj["priors"];
  double rho_a = priors["rho_a"];
  double rho_b = priors["rho_b"];
  vec rho_sd = priors["rho_sd"];
  vec rho_acpt = priors["rho_acpt"];
  Rcpp::List sp_data = RSTr_obj["sp_data"];
  field<uvec> adjacency = sp_data["adjacency"];
  uword n_island = sp_data["n_island"];
  uword n_region = Z.n_rows;
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  vec logit_rho = log(rho / (1 - rho));
  vec rand = Rcpp::rnorm(n_group, 0, 1);
  vec expit_rho = rand % rho_sd + logit_rho;
  vec rho_star_0 = 1 / (1 + exp(-expit_rho));
  vec r(n_group, arma::fill::zeros);
  vec ra(n_group, arma::fill::zeros);
  vec rb(n_group, arma::fill::zeros);
  vec rc(n_group, arma::fill::zeros);
  cube Zm(n_region, n_group, n_time);
  for (uword reg = 0; reg < n_region; reg++) {
    Zm.row(reg) = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
  }
  for (uword grp = 0; grp < n_group; grp++) {
    vec rho_star = rho;
    rho_star[grp] = rho_star_0[grp];
    ra[grp] = (1 - pow(rho[grp], 2)) / (1 - pow(rho_star[grp], 2));
    field<mat> Sein_rho = Sig_eta_i(G, rho);
    field<mat> Sein_rho_star = Sig_eta_i(G, rho_star);
    field<mat> Sein_diff(n_time, n_time);    
    for (uword time1 = 0; time1 < n_time; time1++) {
      uword time1_l = (time1 == 0) ? 0 : time1 - 1;
      uword time2_u = (time1 == n_time - 1) ? n_time : time1 + 2;
      for (uword time2 = time1_l; time2 < time2_u; time2++) {
        Sein_diff(time2, time1) = Sein_rho_star(time2, time1) - Sein_rho(time2, time1);
      }
    }
    for (uword reg = 0; reg < n_region; reg++) {
      for (uword time1 = 0; time1 < n_time; time1++) {
        uword time2_l = (time1 == 0) ? 0 : time1 - 1;
        uword time2_u = (time1 == n_time - 1) ? n_time : time1 + 2;
        for (uword time2 = time2_l; time2 < time2_u; time2++) {
          mat rb_ik = get_grp(Z, reg, time2).t() * Sein_diff(time2, time1) * get_grp(Zm, reg, time1);
          rb[grp] += adjacency[reg].n_elem * rb_ik[0] / 2;
        }
      }
    }
    rc[grp] = pow(rho_star[grp] / rho[grp], rho_a) * pow((1 - rho_star[grp]) / (1 - rho[grp]), rho_b);
    r[grp] = exp((n_region - n_island) * (n_time - 1) / 2 * log(ra[grp]) - rb[grp]) * rc[grp];
    if (r[grp] >= R::runif(0, 1)) {
      rho[grp] = rho_star[grp];
      rho_acpt[grp]++;
    }
  }
  priors["rho_acpt"] = rho_acpt;
  RSTr_obj["priors"] = priors;
  sample["rho"] = rho;
  RSTr_obj["sample"] = sample;
}
