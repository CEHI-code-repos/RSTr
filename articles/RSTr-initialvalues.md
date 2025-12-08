# 08: Initial Values

## Overview

`initial_values` is a `list` specifying the starting values for
parameters in the model. By default, all of these values are specified
according to the literature, but `RSTr` allows the capability of
specifying your own initial values. If you wish to provide initial
values, note that you don’t have to specify `initial_values` for *all*
parameters if you only want to specify some of them - any undefined
`initial_values` will be defined by the default values. For example, you
can specify only the initial values for `lambda` and all other values
will be generated on their own. However, if one value is specified for a
certain parameter in `initial_values`, all values must be specified for
that parameter in `initial_values`: you cannot, for example, define
initial values for just one year of `lambda`. Finally, any values
included in your `initial_values` list that aren’t aligned with the
above names will be ignored.

## Initial value specifications

The models in `RSTr` share many `initial_values`, but a couple of models
have `initial_values` that are unique to them. All potential initial
values are presented here.

### Initial values for the MSTCAR model

Here are the possible initial value parameters for the MSTCAR model:

- `lambda`: The estimated spatially smoothed rate for each
  region-group-time. `lambda` is an `array` of real numbers with
  dimensions `num_region x num_group x num_time`. Has support `(0, 1)`
  for `method = "binomial"` and support `(0, Inf)` for \`method =
  “poisson”;

- `beta`: The mean rate for each island-group-year on a logit- or
  log-transformed scale. Islands are sets of regions that exclusively
  share adjacency information. For example, in `miadj`, there are two
  islands that represent the counties of the Upper Peninsula and the
  Lower Peninsula. These islands don’t touch each other, and thus don’t
  share adjacency information. Each island is assigned its own `beta`.
  `beta` is an `array` of real numbers with dimensions
  `num_island x num_group x num_time`;

- `Z`: The spatiotemporal random effects. These are the parameters that
  induce smoothing on the counties, with the intensity of the smoothing
  dictated by the spatial covariance matrices `G`. `Z` is an `array` of
  real numbers with dimensions `num_region x num_group x num_time`;

- `G`: The spatial covariance matrices. This parameter determines the
  intensity of the spatial smoothing performed by `Z` and represents the
  strength of the relationship between each group in a given time
  period. `G` is an `array` of temporally-evolving positive-definite
  symmetric matrices with dimensions `num_group x num_group x num_time`;

- `rho`: The temporal correlation. This parameter decides the strength
  of the relationship between values in time period `t` to values in
  time period `t-1`. It is a `matrix` of size `num_group x 1` of real
  numbers with support `[0,1]`;

- `tau2`: The non-spatial variance. This parameter picks up any variance
  in values of `lambda` for each group. It is a `matrix` of size
  `num_group x 1` of positive real numbers; and

- `Ag`: The general spatial covariance matrix. This parameter describes
  the overall relationship between groups across the entire model and is
  used in the prior distribution for the matrices in `G`. `Ag` is a
  positive-definite symmetric matrix with dimensions
  `num_group x num_group`.

### Initial values for the MCAR model

The MCAR model utilizes a majority of the initial values of the MSTCAR
model. However, MCAR does not include `initial_values` for `rho` or
`Ag`. Note that specification for the MCAR model is slightly different
than that of the MSTCAR model. If an MCAR model is run with data
containing several time periods, `tau2` will require values for every
time period along with every group.

### Initial values for the UCAR/EUCAR model

The UCAR models have the smallest set of initial values, using only
`lambda`, `beta`, `Z`, and `tau2` from the MCAR model. Similar to the
MCAR, if a UCAR model is run with multiple groups and time periods,
`tau2` requires values for every group and time period present. The only
new initial value for the UCAR models is `sig2`, which takes the place
of `G` in the MCAR and MSTCAR models:

- `sig2` represents the spatial variance of a UCAR/EUCAR model. This
  parameter picks up any variance in values of `Z` for each group. It is
  a `matrix` of size `num_group x num_time` of positive real numbers.
