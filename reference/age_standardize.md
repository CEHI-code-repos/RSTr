# Age-standardize model objects

Age-standardizes samples using a standard population for an `RSTr` model
object.

## Usage

``` r
age_standardize(RSTr_obj, std_pop, new_name, groups = NULL)
```

## Arguments

- RSTr_obj:

  An `RSTr` model object.

- std_pop:

  A vector of standard populations.

- new_name:

  The name to assign to the age-standardized group.

- groups:

  A vector of either indices for each group or a vector of strings for
  each group name. If set to `NULL`, will use all groups in the dataset.

## Value

An `RSTr` object with age-standardized estimates.

## Examples

``` r
std_pop <- c(113154, 100640, 95799)
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
# age-standardize by all age groups
mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
# Add onto age-standardized estimates. Age-standardize only by the first two age groups
mod_mst <- age_standardize(mod_mst, std_pop[1:2], "35-54", groups = 1:2)
```
