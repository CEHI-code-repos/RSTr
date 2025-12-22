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

mat get_subgrp(const cube& arr, const uvec& ind, uword time) {
  return arr.slice(time).rows(ind);
}

vec get_subregs(const cube& arr, const uvec& ind, uword grp,
                uword time) {
  vec col = arr.slice(time).col(grp);
  return col.elem(ind);
}

vec get_grp(const cube& arr, uword reg, uword time) {
  return arr.slice(time).row(reg).t();
}

vec get_row(const cube& arr, uword grp,
                  uword time) {
  return arr.slice(time).col(grp);
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

mat cpp_rmvnorm(const vec& mean, const mat& covar) {
  vec out  = mean;
  const vec rand = Rcpp::rnorm(covar.n_cols, 0, 1);
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

mat irgamma_mat(const double shape, const mat& scale) {
  const uword nr = scale.n_rows;
  const uword nc = scale.n_cols;
  mat x(nr, nc, arma::fill::none);
  for (uword r = 0; r < nr; r++) {
    for (uword c = 0; c < nc; c++) {
      x(r, c) = 1.0 / R::rgamma(shape, scale(r, c));
    }
  }
  return x;
}

mat irtgamma_mat(const double shape, const mat& scale, const mat& thres) {
  const uword nr = scale.n_rows;
  const uword nc = scale.n_cols;
  mat x(nr, nc, arma::fill::none);
  for (uword r = 0; r < nr; r++) {
    for (uword c = 0; c < nc; c++) {
      const double max = R::pgamma(thres(r, c), shape, scale(r, c), true, false);
      double u = R::runif(0, max);
      x(r, c) = 1.0 / R::qgamma(u, shape, scale(r, c), true, false);
    }
  }
  return x;
}

mat rtnorm_mat(const mat& mean, const mat& sd, const mat& thres) {
  const uword nr = mean.n_rows;
  const uword nc = mean.n_cols;
  mat x(nr, nc, arma::fill::none);
  for (uword r = 0; r < nr; ++r) {
    for (uword c = 0; c < nc; ++c) {
      double max = R::pnorm(thres(r, c), mean(r, c), sd(r, c), true, false);
      if (max > 0) {
        double u = R::runif(0, max);
        x(r, c) = R::qnorm(u, mean(r, c), sd(r, c), true, false);
      }
    }
  }
  return x;
}

mat rnorm_mat(const mat& mean, const mat& sd) {
  const uword nr = mean.n_rows;
  const uword nc = mean.n_cols;
  mat x(nr, nc, arma::fill::none);
  for (uword r = 0; r < nr; ++r) {
    for (uword c = 0; c < nc; ++c) {
      x(r, c) = R::rnorm(mean(r, c), sd(r, c));
    }
  }
  return x;
}
