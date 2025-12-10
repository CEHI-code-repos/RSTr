# Suppress estimates based on reliability criteria

Generates suppressed estimates for an `RSTr` model object with a given
relative precision and population/event threshold.

## Usage

``` r
suppress_estimates(RSTr_obj, threshold = 0, type = c("population", "event"))
```

## Arguments

- RSTr_obj:

  An `RSTr` model object.

- threshold:

  The population/event suppression threshold.

- type:

  Determines whether suppression threshold is based on population counts
  or event counts.

## Value

An `RSTr` model object with suppressed estimates.

## Details

While the `threshold` argument is optional, population/event thresholds
are necessary for non-enhanced models. Population/event thresholds
should only be omitted for enhanced CAR models, such as the EUCAR.

## Examples

``` r
std_pop <- c(113154, 100640, 95799)
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
estimates_table <- get_estimates(mod_mst)
mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
estimates_table_as <- get_estimates(mod_mst)
```
