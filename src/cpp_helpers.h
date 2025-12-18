// cpp_helpers.h
#ifndef CPP_HELPERS_H
#define CPP_HELPERS_H

#include <RcppArmadillo.h>

arma::cube get_regs(const arma::cube& arr, const arma::uvec& ind);
arma::vec get_subregs(const arma::cube& arr, const arma::uvec& ind, arma::uword grp, arma::uword time);
arma::vec get_grp(const arma::cube& arr, arma::uword reg, arma::uword time);
arma::mat get_subgrp(const arma::cube& arr, const arma::uvec& ind, arma::uword time);
arma::vec get_row(const arma::cube& arr, arma::uword grp, arma::uword time);
arma::field<arma::mat> Sig_eta_i(const arma::cube& G, const arma::vec& rho);
arma::field<arma::mat> Sig_eta(const arma::field<arma::mat>& Sein);
arma::mat cpp_rmvnorm(const arma::vec& mean, const arma::mat& covar);
arma::mat geteig(const arma::mat& covar);

#endif // CPP_HELPERS_H
