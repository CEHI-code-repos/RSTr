# 09: Priors

## Overview

`priors` is a `list` specifying the priors for various hyperparameters
and auxiliary values in the model. By default, all of these values are
specified according to the literature, but `RSTr` allows the capability
of specifying your own priors. If you wish to provide priors, note that
you don’t have to specify `priors` for *all* parameters if you only want
to specify some of them - any undefined `priors` will be defined by the
default values. For example, you can specify only the priors for
`lambda_sd` and all other values will be generated on their own.
However, if one value is specified for a certain parameter in `priors`,
all values must be specified for that parameter in `priors`: you cannot,
for example, define priors for just one year of `lambda_sd`. Finally,
any values included in your `priors` list that aren’t aligned with the
above names will be ignored.

## Prior specifications

The models in `RSTr` share many `priors`, but a couple of models have
`initial_values` that are unique to them. All potential priors are
presented here.

### Priors for the MSTCAR model

The following are all priors used in the MSTCAR model:

- `Ag_scale` and `Ag_df`: These are the scale and degrees of freedom
  priors used with
  [Wishart-distributed](https://en.wikipedia.org/wiki/Wishart_distribution)
  random variable `Ag`. `Ag_scale` is a positive-definite symmetric
  matrix and `Ag_df` is a `double` of at least size `num_group`;

- `G_scale` and `G_df`: These are the scale and degrees of freedom
  priors used with [Inverse-Wishart
  distributed](https://en.wikipedia.org/wiki/Inverse-Wishart_distribution)
  matrix slices of random variable `G`. `G_scale` is a positive-definite
  symmetric matrix and `G_df` is a `double` of at least size
  `num_group`;

- `tau_a` and `tau_b`: These are the rate and scale priors used with
  [Inverse-Gamma
  distributed](https://en.wikipedia.org/wiki/Inverse-gamma_distribution)
  random variable `tau2`. `tau_a` and `tau_b` must both be positive real
  numbers;

- `rho_a` and `rho_b`: These are the shape priors used with
  [Beta-distributed](https://en.wikipedia.org/wiki/Beta_distribution)
  random variable `rho`. `rho_a` and `rho_b` must both be positive real
  numbers;

- `lambda_sd`: An array of positive real numbers describing the
  candidate standard deviation in the Metropolis update for the
  estimated rates `lambda`. These values will be adaptively updated at
  the start of each batch; and

- `rho_sd`: A vector of positive real numbers describing the candidate
  standard deviation in the Metropolis update for the temporal
  correlation `rho`. These values will be adaptively updated at the
  start of each batch. Note that this is only used if
  `update_rho = TRUE`.

### Priors for the MCAR model

The MCAR model shares all of the priors as the MSTCAR model, but does
not include the following: `Ag_scale`, `Ag_df`, `rho_a`, `rho_b`,
`rho_sd`.

### Priors for the UCAR/EUCAR model

The UCAR models include only the following from above: `lambda_sd`,
`tau_a`, and `tau_b`. UCAR models also take priors `sig_a` and `sig_b`,
which hold similar shape and restriction to `tau_a` and `tau_b`.
