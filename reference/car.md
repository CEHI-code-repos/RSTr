# Create CAR model

`*car()` generates an `RSTr` model object, samples, and estimates for
either an MSTCAR, MCAR, RCAR, or CAR model.

## Usage

``` r
car(
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
  impute_bounds = NULL,
  inits = NULL,
  priors = NULL
)

rcar(
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
  impute_bounds = NULL,
  inits = NULL,
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
  impute_bounds = NULL,
  inits = NULL,
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
  impute_bounds = NULL,
  inits = NULL,
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

- impute_bounds:

  If counts are suppressed for privacy reasons, `impute_bounds` is the
  lower/upper bound of suppression, typically 0 or 1 and 10,
  respectively.

- inits:

  Optional list of initial conditions for each parameter.

- priors:

  Optional list of priors for updates.

- A:

  For RCAR models, describes maximum intensity of smoothing between
  regions.

- m0:

  For RCAR models, baseline neighbor count by region.

- update_rho:

  For MSTCAR models, controls whether rho update is performed.

## Value

An `RSTr` model object.

## Examples

``` r
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
# MSTCAR model
on.exit(unlink(file.path(tempdir(), "test"), recursive = TRUE), add = TRUE)
mod_mst <- mstcar(
  name = "test",
  data = data_min,
  adjacency = adj_min,
  dir = tempdir(),
  show_plots = FALSE,
  verbose = FALSE
)
```
