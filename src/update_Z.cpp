#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_Z_ucar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube Z = sample["Z"];
  mat sig2 = sample["sig2"];
  cube lambda = sample["lambda"];
  cube beta = sample["beta"];
  mat tau2 = sample["tau2"];
  List sp_data = RSTr_obj["sp_data"];
  field<uvec> adjacency = sp_data["adjacency"];
  vec n_adj = sp_data["n_adj"];
  field<uvec> isl_region = sp_data["isl_region"];
  uvec isl_id = sp_data["isl_id"];
  uword n_region = Z.n_rows;
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  uword n_island = isl_region.n_elem;
  for (uword reg = 0; reg < n_region; reg++) {
    for (uword grp = 0; grp < n_group; grp++) {
      for (uword time = 0; time < n_time; time++) {
        double sum_adj = sum(get_subregs(Z, adjacency(reg), grp, time));
        double rate_diff = (lambda(reg, grp, time) - beta(isl_id(reg), grp, time));
        double var_Z  = 1 / (1 / tau2(grp, time) + n_adj(reg) / sig2(grp, time));
        double mean_Z = var_Z * (rate_diff / tau2(grp, time) + sum_adj / sig2(grp, time));
        Z(reg, grp, time) = R::rnorm(mean_Z, sqrt(var_Z));
      }
    }
  }
  cube Zkt(n_island, n_group, n_time);
  for (uword isl = 0; isl < n_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, isl_region(isl)), 0);
  }
  Z -= get_regs(Zkt, isl_id);
  sample["Z"] = Z;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_Z_mcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube Z = sample["Z"];
  cube G = sample["G"];
  cube lambda = sample["lambda"];
  cube beta = sample["beta"];
  mat tau2 = sample["tau2"];
  List sp_data = RSTr_obj["sp_data"];
  field<uvec> adjacency = sp_data["adjacency"];
  vec n_adj = sp_data["n_adj"];
  field<uvec> isl_region = sp_data["isl_region"];
  uvec isl_id = sp_data["isl_id"];
  uword n_region = Z.n_rows;
  uword n_group = Z.n_cols;
  uword n_time = Z.n_slices;
  uword n_island = isl_region.n_elem;
  field<mat> Z_cov(max(n_adj) + 1, n_time);
  field<mat> Z_coveig(max(n_adj) + 1, n_time);
  cube rate_diff = lambda - get_regs(beta, isl_id);
  vec unique_n_adj = unique(n_adj);
  for (uword time = 0; time < n_time; time++) {
    vec taut = tau2.col(time);
    mat Gt = G.slice(time);
    for (uword count : unique_n_adj) {
      Z_cov(count, time) = inv(diagmat(1 / taut) + count * Gt);
      Z_coveig(count, time) = geteig(Z_cov(count, time));
    }
    for (uword reg = 0; reg < n_region; reg++) {
      vec sum_adj = sum(get_subgrp(Z, adjacency(reg), time), 0).t();
      vec Z_mean = Z_cov(n_adj(reg), time) * (get_grp(rate_diff, reg, time) / taut + Gt * sum_adj);
      vec Z_new  = cpp_rmvnorm(Z_mean, Z_coveig(n_adj(reg), time));
      for (uword grp = 0; grp < n_group; grp++) {
        Z(reg, grp, time) = Z_new(grp);
      }
    }
  }
  cube Zkt(n_island, n_group, n_time);
  for (uword isl = 0; isl < n_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, isl_region[isl]), 0);
  }
  Z -= get_regs(Zkt, isl_id);
  sample["Z"] = Z;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_Z_mstcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  cube Z = sample["Z"];
  cube G = sample["G"];
  cube lambda = sample["lambda"];
  cube beta = sample["beta"];
  vec rho = sample["rho"];
  vec tau2 = sample["tau2"];
  List sp_data = RSTr_obj["sp_data"];
  field<uvec> adjacency = sp_data["adjacency"];
  vec n_adj = sp_data["n_adj"];
  field<uvec> isl_region = sp_data["isl_region"];
  uvec isl_id = sp_data["isl_id"];
  uword n_region  = Z.n_rows;
  uword n_group   = Z.n_cols;
  uword n_time    = Z.n_slices;
  uword n_island  = isl_region.n_elem;
  field<mat> Sein   = Sig_eta_i(G, rho);
  field<mat> SeSein = Sig_eta(Sein);
  field<mat> Z_cov(n_time, max(n_adj) + 1);
  field<mat> Z_coveig(n_time, max(n_adj) + 1);
  vec unique_n_adj = unique(n_adj);
  for (uword time = 0; time < n_time; time++) {
    for (uword count : unique_n_adj) {
      Z_cov   (time, count) = inv(diagmat(1 / tau2) + count * Sein(time, time));
      Z_coveig(time, count) = geteig(Z_cov(time, count));
    }
  }
  cube rate_diff = lambda - get_regs(beta, isl_id);
  for (uword reg = 0; reg < n_region; reg++) {
    mat nZm = mean(get_regs(Z, adjacency[reg]), 0);
    if (n_time == 1) nZm = nZm.t();
    for (uword time = 0; time < n_time; time++) {
      vec muZp = nZm.col(time);
      if (time > 0) {
        muZp += SeSein(time, time - 1) * (nZm.col(time - 1) - get_grp(Z, reg, time - 1));
      }
      if (time < n_time - 1) {
        muZp += SeSein(time, time + 1) * (nZm.col(time + 1) - get_grp(Z, reg, time + 1));
      }
      mat Z_mean = Z_cov(time, n_adj[reg]) * (get_grp(rate_diff, reg, time) / tau2 + (n_adj[reg] * Sein(time, time) * muZp));
      vec Z_new  = cpp_rmvnorm(Z_mean, Z_coveig(time, n_adj[reg]));
      for (uword grp = 0; grp < n_group; grp++) {
        Z(reg, grp, time) = Z_new(grp);
      }
    }
  }
  cube Zkt(n_island, n_group, n_time);
  for (uword isl = 0; isl < n_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, isl_region[isl]), 0);
  }
  Z -= get_regs(Zkt, isl_id);
  sample["Z"] = Z;
  RSTr_obj["sample"] = sample;
}
