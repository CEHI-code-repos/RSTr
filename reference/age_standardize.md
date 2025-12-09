# Age-standardize samples

Age-standardizes samples using a standard population for either an
`RSTr` model object or an `array` of samples.

## Usage

``` r
age_standardize(
  object,
  std_pop,
  new_name = NULL,
  groups = NULL,
  margin = NULL,
  bind_new = FALSE
)
```

## Arguments

- object:

  An `RSTr` model object or an `array` of samples imported with
  [`load_samples()`](load_samples.md).

- std_pop:

  A vector of standard populations.

- new_name:

  The name to assign to the age-standardized group.

- groups:

  A vector of either indices for each group or a vector of strings for
  each group name. If set to `NULL`, will use all groups in the dataset.

- margin:

  For `array`s, The margin on which the groups of interest are
  stratified.

- bind_new:

  If set to `TRUE`, will bind an `array` to the original sample dataset.
  Otherwise, will generate a standalone array of samples.

## Value

An `RSTr` object or an `array` of age-standardized samples.

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

std_pop <- c(113154, 100640, 95799)
age_margin <- 2
# age-standardize by all age groups
samples_3564 <- age_standardize(minsample, std_pop, margin = age_margin)
# age-standardize only by the first two age groups
samples_3554 <- age_standardize(minsample, std_pop[1:2], groups = 1:2, margin = age_margin)
# bind age-standardized samples to original samples
samples_as <- age_standardize(
  minsample,
  std_pop,
  new_name = "35-64",
  margin = age_margin,
  bind_new = TRUE
)
```
