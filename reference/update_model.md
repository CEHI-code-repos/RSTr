# Update model

`update_model()` generates additional samples for model `RSTr_obj`.

## Usage

``` r
update_model(RSTr_obj, iterations = 6000, show_plots = TRUE, verbose = TRUE)
```

## Arguments

- RSTr_obj:

  The `RSTr` model object to generate samples for.

- iterations:

  Number of iterations to run.

- show_plots:

  If set to `FALSE`, hides traceplots.

- verbose:

  If set to `FALSE`, hides progress bar and other messages.

## Value

An `RSTr` model object.

## Examples

``` r
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
mod_mst <- update_model(mod_mst, iterations = 1000, show_plots = FALSE, verbose = FALSE)
```
