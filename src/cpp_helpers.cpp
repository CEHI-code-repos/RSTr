#include <RcppArmadillo.h>
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::field;
using arma::uvec;

cube get_regs(const cube& arr, const uvec& ind) {
  cube out(ind.n_elem, arr.n_cols, arr.n_slices);
  for (uword reg = 0; reg < arr.n_slices; ++reg) {
    out.slice(reg) = arr.slice(reg).rows(ind);
  }
  return out;
}

mat get_subgrp(const cube& arr, const uvec& ind, const uword& time) {
  return arr.slice(time).rows(ind);
}

vec get_subregs(const cube& arr, const uvec& ind, const uword& grp,
                const uword& time) {
  vec col = arr.slice(time).col(grp);
  return col.elem(ind);
}

vec get_grp(const cube& arr, const uword& reg, const uword& time) {
  return arr.slice(time).row(reg).t();
}

vec get_row(const cube& arr, const uword& grp,
                  const uword& time) {
  return arr.slice(time).col(grp);
}

field<mat> Sig_eta_i(const cube& G, const vec& rho) {
  uword n_group = rho.n_elem;
  uword n_time  = G.n_slices;
  mat r  = arma::repmat(rho, 1, n_group);
  mat sr = arma::sqrt(1 - pow(r, 2));
  field<mat> Sei(n_time, n_time);
  Sei(0, 0) = arma::inv_sympd(G.slice(0));
  for (uword time = 1; time < n_time; time++) {
    mat Gi = arma::inv_sympd(G.slice(time));
    Sei(time - 1, time - 1) += ( r / sr).t() % (r / sr % Gi);
    Sei(time    , time    )  = ( 1 / sr).t() % (1 / sr % Gi);
    Sei(time - 1, time    )  = (-r / sr)     % (1 / sr % Gi).t();
    Sei(time    , time - 1)  = (-r / sr).t() % (1 / sr % Gi);
  }
  return Sei;
}

field<mat> Sig_eta(const field<mat>& Sein) {
  uword n_time = Sein.n_rows;
  field<mat> Se(n_time, n_time);
  for (uword time = 0; time < n_time; time++) {
    if (time > 0) {
      Se(time, time - 1) = arma::inv_sympd(Sein(time, time)) * Sein(time, time - 1);
    }
    if (time < n_time - 1) {
      Se(time, time + 1) = arma::inv_sympd(Sein(time, time)) * Sein(time, time + 1);
    }
  }
  return Se;
}

mat cpp_rmvnorm(const vec& mean, const mat& covar) {
  vec out  = mean;
  vec rand = Rcpp::rnorm(covar.n_cols, 0, 1);
  out += covar * rand;
  return out;
}

mat geteig(const mat& covar) {
  vec eigval;
  mat eigvec;
  arma::eig_sym(eigval, eigvec, covar);
  eigvec *= eigvec.t() % arma::repmat(sqrt(eigval), 1, covar.n_cols);
  return eigvec.t();
}
