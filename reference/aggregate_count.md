# Aggregate count arrays

Sums counts over event/population arrays. Useful when manually
generating group-aggregated/age-standardized estimates and a population
threshold is needed for suppression.

## Usage

``` r
aggregate_count(
  count,
  margin,
  groups = NULL,
  bind_new = FALSE,
  new_name = NULL
)
```

## Arguments

- count:

  The `array` to aggregate.

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

An `array` of aggregated count data.

## Examples

``` r
margin_time <- 3
# aggregate population from all years for each county-group
pop_7988 <- aggregate_count(miheart$n, margin_time)
# aggregate population from 1980-1984 for each county-group
pop_8084 <- aggregate_count(miheart$n, margin_time, groups = as.character(1980:1984))
# bind aggregated pop from all years to population data
pop_agg <- aggregate_count(miheart$n, margin_time, bind_new = TRUE, new_name = "1979-1988")
```
