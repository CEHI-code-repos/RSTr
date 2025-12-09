# Extract estimates from RSTr model object

Gathers model and estimate information for an `RSTr` model object,
exported as a long table. Estimate rates and their respective credible
intervals are displayed by default in rates per 100,000.

## Usage

``` r
get_estimates(RSTr_obj, rates_per = 1e+05, standardized = TRUE)
```

## Arguments

- RSTr_obj:

  An `RSTr` model object.

- rates_per:

  The desired scaling for estimate rates.

- standardized:

  If `RSTr_obj` contains age-standardized rates, shows the
  age-standardized rates. If set to `FALSE`, always shows the
  non-age-standardized rates.

## Value

A long `table` containing region/group/time period names, estimates,
credible intervals, relative precisions, and the associated
event/population counts.

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
