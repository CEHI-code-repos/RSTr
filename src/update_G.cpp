#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;
using Rcpp::List;

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
  for (uword time = 0; time < n_time; time++) {
    mat scale_G = G_scale;
    for (uword reg = 0; reg < n_region; reg++) {
      const vec Zit = get_grp(Z, reg, time);
      const vec sum_adj = sum(get_subgrp(Z, adjacency(reg), time), 0).t();
      scale_G += adjacency(reg).n_elem * Zit * Zit.t() - sum_adj * Zit.t();
    }
    G.slice(time) = riwish(df_G, scale_G);
  }
  sample["G"] = G;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_G_mstcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube G = sample["G"];
  cube Z = sample["Z"];
  mat Ag = sample["Ag"];
  vec rho = sample["rho"];
  List priors = RSTr_obj["priors"];
  double G_df = priors["G_df"];
  List sp_data = RSTr_obj["sp_data"];
  field<uvec> adjacency = sp_data["adjacency"];
  uword n_island = sp_data["n_island"];
  uword n_region = Z.n_rows;
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  cube Ags(n_group, n_group, n_time, arma::fill::zeros);
  Ags.each_slice() += Ag;
  vec r  = rho;
  vec sr = sqrt(1 - pow(rho, 2));

  for (uword reg = 0; reg < n_region; reg++) {
    double n_adj = adjacency[reg].n_elem;
    mat Zmikt = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
    if (n_time == 1) Zmikt = Zmikt.t();
    mat Zt = get_grp(Z, reg, 0).t();
    Ags.slice(0) += n_adj * Zmikt.col(0) * Zt;
    for (uword time = 1; time < n_time; time++) {
      vec Zt  = 1 / sr % get_grp(Z, reg, time);
      vec Ztl = r / sr % get_grp(Z, reg, time - 1);
      vec Zm  = 1 / sr % Zmikt.col(time);
      vec Zml = r / sr % Zmikt.col(time - 1);
      Ags.slice(time) += n_adj * ((Zm - Zml) * (Zt - Ztl).t());
    }
  }
  for (uword time = 0; time < n_time; time++) {
    G.slice(time) = riwish((n_region - n_island) + G_df, Ags.slice(time));
  }
  sample["G"] = G;
  RSTr_obj["sample"] = sample;
}
