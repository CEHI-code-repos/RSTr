#include <RcppArmadillo.h>
#include <RcppDist.h>
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

inline mat get_mean_Z_car(const cube& Z, const cube& lambda, const cube& beta, 
                   const mat& tau2, const mat& sig2, const field<uvec>& adjacency,
                   const uvec& isl_id, const mat& var_Z, const uword reg) {
  mat sum_adj = arma::sum(get_regs(Z, adjacency[reg]));
  mat lmb_i = lambda.row(reg) - beta.row(isl_id[reg]);
  if (Z.n_slices == 1) {
    sum_adj = sum_adj.t();
    lmb_i = lmb_i.t();
  }
  const mat mean_Z = var_Z % (lmb_i / tau2 + sum_adj / sig2);
  return mean_Z;
}

inline vec get_lmb_it(const cube& lambda, const cube& beta, const uvec& isl_id,
                      const uword time, const uword reg) {
  return (get_grp(lambda, reg, time) - get_grp(beta, isl_id[reg], time));
}

inline mat get_nZm(const cube& Z, const uvec& regs) {
  const uword n_time = Z.n_slices;
  mat nZm = arma::mean(get_regs(Z, regs), 0);
  if (n_time == 1) {
    return nZm.t();
  }
  return nZm;
}

inline vec get_muZp(const mat& nZm, const field<mat>& Se, const cube& Z,
                    const uword reg, const uword time) {
  const uword n_time = Z.n_slices;
  vec muZp = nZm.col(time);
  if (time > 0) {
   muZp += Se(time, time - 1) * (nZm.col(time - 1) - get_grp(Z, reg, time - 1));
  }
  if (time < n_time - 1) {
    muZp += Se(time, time + 1) * (nZm.col(time + 1) - get_grp(Z, reg, time + 1));
  }
  return muZp;
}

inline field<mat> get_cov_Z_mcar(const mat& tau2, const mat& G, const vec& n_adj) {
  const vec unique_n_adj = arma::unique(n_adj);
  field<mat> cov_Z(arma::max(n_adj) + 1);
  const mat it2_diag = arma::diagmat(1.0 / tau2);
  for (uword count : unique_n_adj) {
    cov_Z(count) = arma::inv_sympd(it2_diag + count * G);
  }
  return cov_Z;
}

inline field<mat> get_cov_Z_mstcar(const field<mat>& Sein, const vec& tau2, 
                            const vec& n_adj) {
  const uword n_time = Sein.n_rows;
  field<mat> cov_Z(n_time, arma::max(n_adj) + 1);
  const vec unique_n_adj = arma::unique(n_adj);
  const mat it2_diag = arma::diagmat(1.0 / tau2);
  for (uword time = 0; time < n_time; ++time) {
    for (uword count : unique_n_adj) {
      cov_Z(time, count) = arma::inv_sympd(it2_diag + count * Sein(time, time));
    }
  }
  return cov_Z;
}

inline mat geteig(const mat& covar) {
  vec eigval;
  mat eigvec;
  arma::eig_sym(eigval, eigvec, covar);
  eigvec *= eigvec.t() % arma::repmat(arma::sqrt(eigval), 1, covar.n_cols);
  return eigvec.t();
}

inline field<mat> get_coveig_Z_mcar(const field<mat>& cov_Z, const vec& n_adj) {
  const vec unique_n_adj = arma::unique(n_adj);
  field<mat> coveig_Z(arma::max(n_adj) + 1);
  for (uword count : unique_n_adj) {
    coveig_Z(count) = geteig(cov_Z(count));
  }
  return coveig_Z;
}

inline field<mat> get_coveig_Z_mstcar(const field<mat>& cov_Z, const vec& n_adj) {
  const uword n_time = cov_Z.n_rows;
  field<mat> coveig_Z(n_time, arma::max(n_adj) + 1);
  const vec unique_n_adj = arma::unique(n_adj);
  for (uword time = 0; time < n_time; ++time) {
    for (uword count : unique_n_adj) {
      coveig_Z(time, count) = geteig(cov_Z(time, count));
    }
  }
  return coveig_Z;
}

inline void demean_Z(cube& Z, const field<uvec>& isl_region, const uvec& isl_id) {
  const uword n_group = Z.n_cols;
  const uword n_time = Z.n_slices;
  const uword n_island = isl_region.n_elem;
  cube Zkt(n_island, n_group, n_time, arma::fill::none);
  for (uword isl = 0; isl < n_island; ++isl) {
    Zkt.row(isl) = arma::mean(get_regs(Z, isl_region[isl]), 0);
  }
  Z -= get_regs(Zkt, isl_id);
}

inline void replace_Zit(cube& Z, const vec& mean_Z, const mat& cov_Z, const uword reg,
                 const uword time) {
  const uword n_group = Z.n_cols;
  vec Z_new = rmvnorm_vec(mean_Z, cov_Z);
  for (uword grp = 0; grp < n_group; ++grp) {
    Z(reg, grp, time) = Z_new[grp];
  }
}

//[[Rcpp::export]]
void update_Z_car(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  auto Z = Rcpp::as<cube>(sample["Z"]);
  const auto sig2 = Rcpp::as<mat>(sample["sig2"]);
  const auto lambda = Rcpp::as<cube>(sample["lambda"]);
  const auto beta = Rcpp::as<cube>(sample["beta"]);
  const auto tau2 = Rcpp::as<mat>(sample["tau2"]);
  const List& sp_data = RSTr_obj["sp_data"];
  const auto adjacency = Rcpp::as<field<uvec>>(sp_data["adjacency"]);
  const auto n_adj = Rcpp::as<vec>(sp_data["n_adj"]);
  const auto isl_region = Rcpp::as<field<uvec>>(sp_data["isl_region"]);
  const auto isl_id = Rcpp::as<uvec>(sp_data["isl_id"]);
  const uword n_region = Z.n_rows;
  
  for (uword reg = 0; reg < n_region; ++reg) {
    const mat var_Z = 1.0 / (1.0 / tau2 + n_adj[reg] / sig2);
    const mat mean_Z = get_mean_Z_car(Z, lambda, beta, tau2, sig2, adjacency, 
                                      isl_id, var_Z, reg);
    Z.row(reg) = rnorm_mat(mean_Z, arma::sqrt(var_Z));
  }
  demean_Z(Z, isl_region, isl_id);

  sample["Z"] = Z;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_Z_mcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  auto Z = Rcpp::as<cube>(sample["Z"]);
  const auto G = Rcpp::as<cube>(sample["G"]);
  const auto lambda = Rcpp::as<cube>(sample["lambda"]);
  const auto beta = Rcpp::as<cube>(sample["beta"]);
  const auto tau2 = Rcpp::as<mat>(sample["tau2"]);
  const List& sp_data = RSTr_obj["sp_data"];
  const auto adjacency = Rcpp::as<field<uvec>>(sp_data["adjacency"]);
  const auto n_adj = Rcpp::as<vec>(sp_data["n_adj"]);
  const auto isl_region = Rcpp::as<field<uvec>>(sp_data["isl_region"]);
  const auto isl_id = Rcpp::as<uvec>(sp_data["isl_id"]);
  const uword n_region = Z.n_rows;
  const uword n_time = Z.n_slices;

  const cube& lmb = lambda - get_regs(beta, isl_id);
  for (uword time = 0; time < n_time; ++time) {
    const vec& tau2t = tau2.col(time);
    const mat& Gt = G.slice(time);
    const field<mat> cov_Z = get_cov_Z_mcar(tau2t, Gt, n_adj);
    const field<mat> coveig_Z = get_coveig_Z_mcar(cov_Z, n_adj);
    for (uword reg = 0; reg < n_region; ++reg) {
      const vec lmb_it = get_lmb_it(lambda, beta, isl_id, time, reg);
      const vec sum_adj = arma::sum(get_subgrp(Z, adjacency[reg], time), 0).t();
      const vec mean_Z = cov_Z(n_adj[reg]) * (lmb_it / tau2t + Gt * sum_adj);
      replace_Zit(Z, mean_Z, coveig_Z(n_adj[reg]), reg, time);
    }
  }
  demean_Z(Z, isl_region, isl_id);

  sample["Z"] = Z;
  RSTr_obj["sample"] = sample;
}

//[[Rcpp::export]]
void update_Z_mstcar(List& RSTr_obj) {
  List sample = RSTr_obj["sample"];
  auto Z = Rcpp::as<cube>(sample["Z"]);
  const auto G = Rcpp::as<cube>(sample["G"]);
  const auto lambda = Rcpp::as<cube>(sample["lambda"]);
  const auto beta = Rcpp::as<cube>(sample["beta"]);
  const auto rho = Rcpp::as<vec>(sample["rho"]);
  const auto tau2 = Rcpp::as<vec>(sample["tau2"]);
  const List& sp_data = RSTr_obj["sp_data"];
  const auto adjacency = Rcpp::as<field<uvec>>(sp_data["adjacency"]);
  const auto n_adj = Rcpp::as<vec>(sp_data["n_adj"]);
  const auto isl_region = Rcpp::as<field<uvec>>(sp_data["isl_region"]);
  const auto isl_id = Rcpp::as<uvec>(sp_data["isl_id"]);
  const uword n_region = Z.n_rows;
  const uword n_time = Z.n_slices;

  const field<mat> Sein = Sig_eta_i(G, rho);
  const field<mat> Se = Sig_eta(Sein);
  const field<mat> cov_Z = get_cov_Z_mstcar(Sein, tau2, n_adj);
  const field<mat> coveig_Z = get_coveig_Z_mstcar(cov_Z, n_adj);
  for (uword reg = 0; reg < n_region; ++reg) {
    const mat nZm = get_nZm(Z, adjacency[reg]);
    for (uword time = 0; time < n_time; ++time) {
      const vec lmb_it = get_lmb_it(lambda, beta, isl_id, time, reg);
      const vec muZp = get_muZp(nZm, Se, Z, reg, time);
      const vec agg_Z = lmb_it / tau2 + (n_adj[reg] * Sein(time, time) * muZp);
      const vec mean_Z = cov_Z(time, n_adj[reg]) * agg_Z;
      replace_Zit(Z, mean_Z, coveig_Z(time, n_adj[reg]), reg, time);
    }
  }
  demean_Z(Z, isl_region, isl_id);

  sample["Z"] = Z;
  RSTr_obj["sample"] = sample;
}
