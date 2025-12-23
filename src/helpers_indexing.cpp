#include <RcppArmadillo.h>
using arma::vec;
using arma::mat;
using arma::cube;
using arma::uword;
using arma::uvec;

cube get_regs(const cube& arr, const uvec& ind) {
  cube out(ind.n_elem, arr.n_cols, arr.n_slices);
  for (uword reg = 0; reg < arr.n_slices; ++reg) {
    out.slice(reg) = arr.slice(reg).rows(ind);
  }
  return out;
}

mat get_subgrp(const cube& arr, const uvec& ind, const uword time) {
  return arr.slice(time).rows(ind);
}

vec get_subregs(const cube& arr, const uvec& ind, const uword grp, const uword time) {
  vec col = arr.slice(time).col(grp);
  return col.elem(ind);
}

vec get_grp(const cube& arr, const uword reg, const uword time) {
  return arr.slice(time).row(reg).t();
}

vec get_row(const cube& arr, const uword grp, const uword time) {
  return arr.slice(time).col(grp);
}
