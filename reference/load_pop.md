# Load population array

Imports the population array from the model directory.

## Usage

``` r
load_pop(name, dir = tempdir())
```

## Arguments

- name:

  Name of model and corresponding folder

- dir:

  Directory where model will live

## Value

The population array for the specified model.

## Examples

``` r
initialize_mstcar("test", miheart, miadj, tempdir(), show_plots = FALSE)
#> Checking data...
#> Checking spatial data...
#> Checking inits...
#> The following objects were created using defaults in 'inits': beta theta Z G rho tau2 Ag
#> Checking priors...
#> The following objects were created using defaults in 'priors': theta_sd tau_a tau_b Ag_scale Ag_df G_df rho_a rho_b rho_sd
#> Model ready!
pop <- load_pop("test", tempdir())
```
