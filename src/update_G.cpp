#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "helpers_indexing.h"
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;

mat get_scale_G(const mat& G_scale, const cube& Z, const field<uvec>& adjacency,
                const uword n_region, const uword time) {
  mat scale_G = G_scale;
  for (uword reg = 0; reg < n_region; ++reg) {
    const vec Zit = get_grp(Z, reg, time);
    const vec sum_adj = arma::sum(get_subgrp(Z, adjacency[reg], time), 0).t();
    scale_G += adjacency[reg].n_elem * Zit * Zit.t() - sum_adj * Zit.t();
  }
  return scale_G;
}

inline mat get_Zmikt(const cube& Z, const uword reg, const uvec& adj_regs,
                     const uword n_time) {
  mat Zmikt = Z.row(reg) - arma::mean(get_regs(Z, adj_regs), 0);
  if (n_time == 1) Zmikt = Zmikt.t();
  return Zmikt;
}

cube get_Ags(const mat& Ag, const cube& Z, const vec& rho,
             const field<uvec>& adjacency) {
  const uword n_region = Z.n_rows;
  const uword n_group = Z.n_cols;
  const uword n_time = Z.n_slices;
  cube Ags(n_group, n_group, n_time, arma::fill::zeros);
  Ags.each_slice() += Ag;
  const vec isr = 1.0 / arma::sqrt(1 - arma::square(rho));
  const vec rsr = rho / arma::sqrt(1 - arma::square(rho));
  for (uword reg = 0; reg < n_region; ++reg) {
    const uvec adj_regs = adjacency[reg];
    const double n_adj = adj_regs.n_elem;
    const mat Zmikt = get_Zmikt(Z, reg, adj_regs, n_time);
    const mat Zt0 = get_grp(Z, reg, 0).t();
    Ags.slice(0) += n_adj * Zmikt.col(0) * Zt0;
    for (uword time = 1; time < n_time; ++time) {
      const vec Zt = get_grp(Z, reg, time);
      const vec Zt1 = get_grp(Z, reg, time - 1);
      const vec Zt_diff = isr % Zt - rsr % Zt1;
      const vec Zm_diff = isr % Zmikt.col(time) - rsr % Zmikt.col(time - 1);
      Ags.slice(time) += n_adj * (Zm_diff * Zt_diff.t());
    }
  }
  //Rcpp::Rcout << Ags << "\n";
  return Ags;
}

//[[Rcpp::export]]
void update_G_default(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube G = sample["G"];
  const cube& Z = sample["Z"];
  const List& priors = RSTr_obj["priors"];
  const double G_df = priors["G_df"];
  const mat& G_scale = priors["G_scale"];
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec>& adjacency = sp_data["adjacency"];
  const uword n_island = sp_data["n_island"];
  const uword n_region = Z.n_rows;
  const uword n_time = Z.n_slices;

  const double df_G = n_region - n_island + G_df;
  for (uword time = 0; time < n_time; ++time) {
    const mat scale_G = get_scale_G(G_scale, Z, adjacency, n_region, time);
    G.slice(time) = riwish(df_G, scale_G);
  }

  sample["G"] = G;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_G_mstcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube G = sample["G"];
  const cube& Z = sample["Z"];
  const mat& Ag = sample["Ag"];
  const vec& rho = sample["rho"];
  const List& priors = RSTr_obj["priors"];
  const double G_df = priors["G_df"];
  const List& sp_data = RSTr_obj["sp_data"];
  const field<uvec>& adjacency = sp_data["adjacency"];
  const uword n_island = sp_data["n_island"];
  const uword n_region = Z.n_rows;
  const uword n_time = Z.n_slices;

  const cube Ags = get_Ags(Ag, Z, rho, adjacency);
  for (uword time = 0; time < n_time; ++time) {
    G.slice(time) = riwish((n_region - n_island) + G_df, Ags.slice(time));
  }
  
  sample["G"] = G;
  RSTr_obj["sample"] = sample;
}
