// helpers_car.h
#ifndef HELPERS_CAR_H
#define HELPERS_CAR_H

#include <RcppArmadillo.h>

arma::mat get_pi_rcar(const arma::cube& beta, const arma::uvec& n_isl_reg, const arma::uword n_reg);
arma::field<arma::mat> Sig_eta_i(const arma::cube& G, const arma::vec& rho);
arma::field<arma::mat> Sig_eta(const arma::field<arma::mat>& Sein);

#endif // HELPERS_CAR_H
