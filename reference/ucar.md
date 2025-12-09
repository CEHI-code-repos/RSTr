# Create CAR model

`*car()` generates an `RSTr` model object, samples, and estimates for
either an MSTCAR, MCAR, EUCAR, or UCAR model.

## Usage

``` r
ucar(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
  priors = NULL
)

eucar(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  A = NULL,
  m0 = NULL,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
  priors = NULL
)

mcar(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
  priors = NULL
)

mstcar(
  name,
  data,
  adjacency,
  dir = tempdir(),
  seed = NULL,
  perc_ci = 0.95,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  initial_values = NULL,
  priors = NULL,
  update_rho = FALSE
)
```

## Arguments

- name:

  Name of model and corresponding folder.

- data:

  Dataset including mortality (Y) and population (n) information.

- adjacency:

  Dataset including adjacency information.

- dir:

  Directory where model will live.

- seed:

  Set of random seeds to use for data replication.

- perc_ci:

  The percentage of the desired estimate credible interval. Defaults to
  95 percent (0.95).

- iterations:

  The number of iterations to run the model for.

- show_plots:

  If set to `FALSE`, suppresses traceplots.

- verbose:

  If set to `FALSE`, suppresses model progress messages.

- ignore_checks:

  If set to `TRUE`, skips model validation.

- method:

  Run model with either Binomial data or Poisson data.

- impute_lb:

  If counts are suppressed for privacy reasons, `impute_lb` is lower
  bound of suppression, typically 0 or 1.

- impute_ub:

  If counts are suppressed for privacy reasons, `impute_ub` is upper
  bound of suppression, typically 10.

- initial_values:

  Optional list of initial conditions for each parameter.

- priors:

  Optional list of priors for updates.

- A:

  For EUCAR models, describes maximum intensity of smoothing between
  regions.

- m0:

  For EUCAR models, baseline neighbor count by region.

- update_rho:

  For MSTCAR models, controls whether rho update is performed.

## Value

An `RSTr` model object.

## Examples

``` r
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
# MSTCAR model
mod_mst <- mstcar(
  name = "test",
  data = data_min,
  adjacency = adj_min,
  dir = tempdir(),
  show_plots = FALSE,
  verbose = FALSE
)
```
