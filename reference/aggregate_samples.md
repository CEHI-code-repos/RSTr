# Aggregate samples by non-age group

Consolidates a set of samples over non-age groups using a population
array to create weighted-average samples.

## Usage

``` r
aggregate_samples(
  sample,
  pop,
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

- pop:

  The population array to be used for weighted averages.

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

An `array` of weighted-average samples.

## Details

`aggregate_samples()` is only meant for non-age group data, such as
spatial regions, time periods, or other sociodemographic groups (race,
sex, etc.). If you are interested in consolidating samples by age group,
use [`age_standardize()`](age_standardize.md) instead. Additionally, if
you plan on doing age-standardization along with aggregating by other
groups, always aggregate groups first before doing age-standardization
to ensure that the samples are properly standardized.

## Examples

``` r
pop <- miheart$n[1:2, 1:3, 1:3]
time_margin <- 3
# calculate prevalence by aggregating over time periods
samples_3564 <- aggregate_samples(minsample, pop, margin = time_margin)
# calculate prevalence of only the first two time periods
samples_3554 <- aggregate_samples(minsample, pop, time_margin, groups = 1:2)
# bind prevalence samples to original samples
samples_prev <- aggregate_samples(
  minsample,
  pop,
  time_margin,
  bind_new = TRUE,
  new_name = "1979-1981"
)
```
