# Age-standardize samples

Age-standardizes samples using a standard population.

## Usage

``` r
standardize_samples(
  sample,
  std_pop,
  margin,
  groups = NULL,
  bind_new = FALSE,
  new_name = NULL
)
```

## Arguments

- sample:

  an `array` of samples imported with
  [`load_samples()`](load_samples.md)

- std_pop:

  A vector of standard populations.

- margin:

  For `array`s, The margin on which the groups of interest are
  stratified.

- groups:

  A vector of either indices for each group or a vector of strings for
  each group name. If set to `NULL`, will use all groups in the dataset.

- bind_new:

  If set to `TRUE`, will bind an `array` to the original sample dataset.
  Otherwise, will generate a standalone array of samples.

- new_name:

  The name to assign to the age-standardized group.

## Value

An `array` of age-standardized samples.

## Examples

``` r
std_pop <- c(113154, 100640, 95799)
age_margin <- 2
# age-standardize by all age groups
samples_3564 <- standardize_samples(minsample, std_pop, age_margin)
# age-standardize only by the first two age groups
samples_3554 <- standardize_samples(minsample, std_pop[1:2], age_margin, groups = 1:2)
# bind age-standardized samples to original samples
samples_as <- standardize_samples(
  minsample,
  std_pop,
  age_margin,
  bind_new = TRUE,
  new_name = "35-64"
)
```
