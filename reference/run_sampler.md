# Run Gibbs sampler

`run_sampler()` generates samples for model `name` in `dir`. The model
used to generate samples (e.g., MSTCAR, MCAR, UCAR) along with the
model's other parameters are specified in `initialize_*()`.

## Usage

``` r
run_sampler(
  name,
  dir = tempdir(),
  iterations = 6000,
  show_plots = TRUE,
  show_progress = TRUE,
  discard_burnin = FALSE
)
```

## Arguments

- name:

  Name of model and corresponding folder

- dir:

  Directory where model lives

- iterations:

  Specifies number of iterations to run

- show_plots:

  If set to `FALSE`, hides traceplots

- show_progress:

  If set to `FALSE`, hides progress bar

- discard_burnin:

  If set to `TRUE`, won't save burn-in samples

## Value

No output, saves sampler output to `dir`

## Examples

``` r
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
initialize_mstcar("test", data_min, adj_min, tempdir())
#> Checking data...

#> Checking spatial data...
#> Checking inits...
#> The following objects were created using defaults in 'inits': beta theta Z tau2 G rho Ag
#> Checking priors...
#> The following objects were created using defaults in 'priors': theta_sd tau_a tau_b G_df G_scale Ag_scale Ag_df rho_a rho_b rho_sd
#> Model ready!
run_sampler("test", show_plots = FALSE, show_progress = FALSE)
#> Starting sampler on Batch 1 at Wed Nov 19 14:56:35
#> Model finished at Wed Nov 19 14:56:37
```
