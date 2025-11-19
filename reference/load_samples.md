# Load MCMC samples

`load_samples()` gathers samples saved for model `name` in directory
`dir`. By default, loads the rate estimate samples `theta`, but any
model parameters can be loaded. Users can also specify a burn-in period.

## Usage

``` r
load_samples(name, dir = tempdir(), param = "theta", burn = 2000)
```

## Arguments

- name:

  Name of model

- dir:

  Directory where model lives

- param:

  Which parameter samples to load

- burn:

  Numer of burn-in samples to discard

## Value

An `array` of samples from model `name`

## Examples

``` r
# prepare truncated dataset
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
initialize_model("test", tempdir(), data_min, adj_min, show_plots = FALSE)
#> Error in initialize_model("test", tempdir(), data_min, adj_min, show_plots = FALSE): could not find function "initialize_model"
run_sampler("test", show_plots = FALSE, show_progress = FALSE)
#> Starting sampler on Batch 1 at Wed Nov 19 20:45:44
#> Model finished at Wed Nov 19 20:46:04
theta <- load_samples("test", tempdir()) * 1e5
```
