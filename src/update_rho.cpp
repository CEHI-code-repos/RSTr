#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <cmath>
#include "helpers_indexing.h"
#include "helpers_car.h"
using arma::vec;
using arma::rowvec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;

field<mat> get_Sein_diff(const cube& G, const vec& rho, const vec& rho_star) {
  const uword n_time = G.n_slices;
  const field<mat> Sein_rho = Sig_eta_i(G, rho);
  const field<mat> Sein_rho_star = Sig_eta_i(G, rho_star);
  field<mat> Sein_diff(n_time, n_time);    
  for (uword t1 = 0; t1 < n_time; ++t1) {
    const uword t1_l = (t1 == 0) ? 0 : t1 - 1;
    const uword t2_u = (t1 == n_time - 1) ? n_time : t1 + 2;
    for (uword t2 = t1_l; t2 < t2_u; ++t2) {
      Sein_diff(t2, t1) = Sein_rho_star(t2, t1) - Sein_rho(t2, t1);
    }
  }
  return Sein_diff;
}

inline double get_diff_quad(const cube& Z, const cube& Zm, 
                            const field<mat>& Sein_diff, const uword reg,
                            const uword t1, const uword t2) {
  const rowvec& Z_it2 = Z.slice(t2).row(reg);
  const mat& Se_t21 = Sein_diff(t2, t1);
  const vec& Zm_it1 = Zm.slice(t1).row(reg).t();
  return arma::as_scalar(Z_it2 * Se_t21 * Zm_it1);
}

double get_rb(const vec& rho, const vec& rho_star, const cube& G, const cube& Z,
              const cube& Zm, const vec& n_adj) {
  const uword n_time = G.n_slices;
  const uword n_region = Z.n_rows;
  const field<mat> Sein_diff = get_Sein_diff(G, rho, rho_star);
  double rb = 0;
  for (uword reg = 0; reg < n_region; ++reg) {
    for (uword t1 = 0; t1 < n_time; ++t1) {
      const uword t2_l = (t1 == 0) ? 0 : t1 - 1;
      const uword t2_u = (t1 == n_time - 1) ? n_time : t1 + 2;
      for (uword t2 = t2_l; t2 < t2_u; ++t2) {
        const double diff_quad = get_diff_quad(Z, Zm, Sein_diff, reg, t1, t2);
        rb += n_adj[reg] * diff_quad / 2;
      }
    }
  }
  return rb;
}

double get_r(const vec& rho, const vec& rho_star, const cube& G, const cube& Z,
             const cube& Zm, const vec& n_adj, const double rho_a,
             const double rho_b, const uword n_island, const uword grp) {
  const uword n_region = Z.n_rows;
  const uword n_time = Z.n_slices;
  const double r2 = rho[grp] * rho[grp];
  const double r_s2 = rho_star[grp] * rho_star[grp];
  const double ra = std::log((1 - r2) / (1 - r_s2));
  const double rb = get_rb(rho, rho_star, G, Z, Zm, n_adj);
  const double rsr = rho_star[grp] / rho[grp];
  const double rsr1m = (1 - rho_star[grp]) / (1 - rho[grp]);
  const double rc = std::pow(rsr, rho_a) * std::pow(rsr1m, rho_b);
  return (std::exp((n_region - n_island) * (n_time - 1) / 2 * ra - rb) * rc);
}

cube get_Zm(const cube& Z, const field<uvec>& adjacency) {
  const uword n_region = Z.n_rows;
  cube Zm(arma::size(Z), arma::fill::none);
  for (uword reg = 0; reg < n_region; ++reg) {
    Zm.row(reg) = Z.row(reg) - arma::mean(get_regs(Z, adjacency[reg]), 0);
  }
  return Zm;
}

vec get_rho_star_0(const vec& rho, const vec& rho_sd) {
  const uword n_group = rho.n_elem;
  const vec logit_rho = arma::log(rho / (1 - rho));
  const vec rand = Rcpp::rnorm(n_group, 0, 1);
  const vec mean_rho = rand % rho_sd + logit_rho;
  return (1.0 / (1 + arma::exp(-mean_rho)));
}

//[[Rcpp::export]]
void update_rho(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  auto rho = Rcpp::as<vec>(sample["rho"]);
  const auto G = Rcpp::as<cube>(sample["G"]);
  const auto Z = Rcpp::as<cube>(sample["Z"]);
  List priors = RSTr_obj["priors"];
  auto rho_acpt = Rcpp::as<vec>(priors["rho_acpt"]);
  const auto rho_sd = Rcpp::as<vec>(priors["rho_sd"]);
  const double rho_a = priors["rho_a"];
  const double rho_b = priors["rho_b"];
  const List& sp_data = RSTr_obj["sp_data"];
  const auto adjacency = Rcpp::as<field<uvec>>(sp_data["adjacency"]);
  const uword n_island = sp_data["n_island"];
  const auto n_adj = Rcpp::as<vec>(sp_data["n_adj"]);
  const uword n_group = Z.n_cols;

  const vec rho_star_0 = get_rho_star_0(rho, rho_sd);
  const cube Zm = get_Zm(Z, adjacency);
  for (uword grp = 0; grp < n_group; grp++) {
    vec rho_star = rho;
    rho_star[grp] = rho_star_0[grp];
    const double r = get_r(rho, rho_star, G, Z, Zm, n_adj, rho_a, rho_b, n_island, grp);
    if (r >= R::runif(0, 1)) {
      rho[grp] = rho_star[grp];
      rho_acpt[grp]++;
    }
  }

  priors["rho_acpt"] = rho_acpt;
  RSTr_obj["priors"] = priors;
  sample["rho"] = rho;
  RSTr_obj["sample"] = sample;
}
