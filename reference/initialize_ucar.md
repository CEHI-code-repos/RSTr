# Initialize CAR model

This function performs checks and prepares data for use with either an
MSTCAR, MCAR, or UCAR model. This function additionally specifies all of
the model parameters, such as model type, event data type, intensity of
smoothing in the UCAR model, and more.

## Usage

``` r
initialize_ucar(
  name,
  data,
  adjacency,
  dir = tempdir(),
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = 1234,
  inits = NULL,
  priors = NULL
)

initialize_ucar_restricted(
  name,
  data,
  adjacency,
  dir = tempdir(),
  A = NULL,
  m0 = NULL,
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = 1234,
  inits = NULL,
  priors = NULL
)

initialize_mcar(
  name,
  data,
  adjacency,
  dir = tempdir(),
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = 1234,
  inits = NULL,
  priors = NULL
)

initialize_mstcar(
  name,
  data,
  adjacency,
  dir = tempdir(),
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = 1,
  impute_ub = 9,
  seed = 1234,
  inits = NULL,
  priors = NULL,
  rho_up = FALSE
)
```

## Arguments

- name:

  Name of model and corresponding folder

- data:

  Dataset including mortality (Y) and population (n) information

- adjacency:

  Dataset including adjacency information

- dir:

  Directory where model will live

- show_plots:

  If set to `FALSE`, suppresses check plots generated for MSTCAR models

- ignore_checks:

  If set to `TRUE`, ignores data checks. Only use if you are certain
  that your input data is correct and you are encountering bugs during
  setup

- method:

  Run model with either Binomial data or Poisson data

- impute_lb:

  If counts are suppressed for privacy reasons, `impute_lb` is lower
  bound of suppression, typically 0 or 1

- impute_ub:

  If counts are suppressed for privacy reasons, `impute_ub` is upper
  bound of suppression, typically 10

- seed:

  Set of random seeds to use for data replication

- inits:

  Optional list of initial conditions for each parameter

- priors:

  Optional list of priors for updates

- A:

  For restricted UCAR models, describes intensity of smoothing between
  regions

- m0:

  For restricted UCAR models, baseline neighbor count by region

- rho_up:

  For MSTCAR models, controls whether rho update is performed for MSTCAR
  models

## Value

No output, only sets up model and saves files to directory

## Examples

``` r
# Initialize an MSTCAR model
initialize_mstcar(name = "test", data = miheart, adjacency = miadj, dir = tempdir())
#> Checking data...

#> Checking spatial data...
#> Checking inits...
#> The following objects were created using defaults in 'inits': beta theta Z G rho tau2 Ag
#> Checking priors...
#> The following objects were created using defaults in 'priors': theta_sd tau_a tau_b Ag_scale Ag_df G_df rho_a rho_b rho_sd
#> Model ready!
# Initialize an MCAR model
data_m <- lapply(miheart, \(x) x[, , "1979"])
initialize_mcar("test", data_m, miadj, tempdir())
#> Checking data...
#> Checking spatial data...
#> Checking inits...
#> The following objects were created using defaults in 'inits': beta theta Z G tau2
#> Checking priors...
#> The following objects were created using defaults in 'priors': theta_sd tau_a tau_b G_scale G_df
#> Model ready!
# Initialize an MCAR model with Poisson-distributed event data
initialize_mcar("test", data_m, miadj, tempdir(), method = "poisson")
#> Checking data...
#> Checking spatial data...
#> Checking inits...
#> The following objects were created using defaults in 'inits': beta theta Z G tau2
#> Checking priors...
#> The following objects were created using defaults in 'priors': theta_sd tau_a tau_b G_scale G_df
#> Model ready!
# Initialize a restricted UCAR model
data_u <- lapply(miheart, \(x) x[, "65-74", "1979"])
initialize_ucar_restricted("test", data_u, miadj, tempdir(), A = 6)
#> Checking data...
#> Checking spatial data...
#> Checking inits...
#> The following objects were created using defaults in 'inits': beta theta Z tau2 sig2
#> Checking priors...
#> Error in if (dim(theta_sd) != dim(data$Y)) {    errct <- errct + 1    errtxt <- paste(errct, ": theta_sd has different length than data. Ensure length(theta_sd) == length(Y) or use default value")    errout <- c(errout, errtxt)}: the condition has length > 1
```
