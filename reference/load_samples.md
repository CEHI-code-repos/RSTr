# Load MCMC samples

`load_samples()` gathers samples saved for model `RSTr_obj`. By default,
loads the rate estimate samples `lambda`, but any model parameters can
be loaded. Users can also specify a burn-in period.

## Usage

``` r
load_samples(RSTr_obj, param = "lambda", burn = 2000)
```

## Arguments

- RSTr_obj:

  `RSTr` model object to load in samples from.

- param:

  Which parameter samples to load.

- burn:

  Number of burn-in samples to discard.

## Value

An `array` of samples from model `RSTr_obj`.

## Examples

``` r
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
samples <- load_samples(mod_mst) * 1e5
```
