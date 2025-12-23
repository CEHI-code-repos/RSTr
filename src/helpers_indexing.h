// helpers_indexing.h
#ifndef HELPERS_INDEXING_H
#define HELPERS_INDEXING_H

#include <RcppArmadillo.h>

arma::cube get_regs(const arma::cube& arr, const arma::uvec& ind);
arma::mat get_subgrp(const arma::cube& arr, const arma::uvec& ind, const arma::uword time);
arma::vec get_subregs(const arma::cube& arr, const arma::uvec& ind, const arma::uword grp, const arma::uword time);
arma::vec get_grp(const arma::cube& arr, const arma::uword reg, const arma::uword time);
arma::vec get_row(const arma::cube& arr, const arma::uword grp, const arma::uword time);

#endif // HELPERS_INDEXING_H
