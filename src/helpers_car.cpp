#include <RcppArmadillo.h>
#include "helpers_indexing.h"
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;

mat get_pi_rcar(const cube& beta, const uvec& n_isl_reg, const uword n_reg) {
  const uword n_group = beta.n_cols;
  const uword n_time = beta.n_slices;
  mat pi(n_group, n_time);
  for (uword grp = 0; grp < n_group; ++grp) {
    for (uword time = 0; time < n_time; ++time) {
      pi(grp, time) = arma::sum(get_row(beta, grp, time) % n_isl_reg / n_reg);
    }
  }
  mat exp_pi = arma::exp(pi);
  return (exp_pi / (1 + exp_pi));
}

field<mat> Sig_eta_i(const cube& G, const vec& rho) {
  const uword n_group = rho.n_elem;
  const uword n_time  = G.n_slices;
  const mat r  = arma::repmat(rho, 1, n_group);
  const mat sr = arma::sqrt(1 - pow(r, 2));
  field<mat> Sei(n_time, n_time);
  Sei(0, 0) = arma::inv_sympd(G.slice(0));
  for (uword time = 1; time < n_time; time++) {
    const mat Gi = arma::inv_sympd(G.slice(time));
    Sei(time - 1, time - 1) += ( r   / sr).t() % (r   / sr % Gi);
    Sei(time    , time    )  = ( 1.0 / sr).t() % (1.0 / sr % Gi);
    Sei(time - 1, time    )  = (-r   / sr)     % (1.0 / sr % Gi).t();
    Sei(time    , time - 1)  = (-r   / sr).t() % (1.0 / sr % Gi);
  }
  return Sei;
}

field<mat> Sig_eta(const field<mat>& Sein) {
  const uword n_time = Sein.n_rows;
  field<mat> Se(n_time, n_time);
  for (uword time = 0; time < n_time; time++) {
    const mat Sinv = arma::inv_sympd(Sein(time, time));
    if (time > 0) {
      Se(time, time - 1) = Sinv * Sein(time, time - 1);
    }
    if (time < n_time - 1) {
      Se(time, time + 1) = Sinv * Sein(time, time + 1);
    }
  }
  return Se;
}
