# Aggregate population arrays

Aggregate population arrays

## Usage

``` r
aggregate_pop(pop, margin, groups = NULL, bind_new = FALSE, new_name = NULL)
```

## Arguments

- pop:

  The population array to aggregate

- margin:

  The margin on which the groups of interest are stratified

- groups:

  A vector of either indices for each group or a vector of strings for
  each group name. If set to `NULL`, will use all groups in the dataset.

- bind_new:

  If set to `TRUE`, will bind to the original sample dataset using
  [`abind::abind()`](https://rdrr.io/pkg/abind/man/abind.html).
  Otherwise, will generate a standalone array of samples.

- new_name:

  If binding the array to the original sample dataset, will assign a
  name to the new margin

## Value

An `array` of aggregated population data

## Examples

``` r
margin_time <- 3
# aggregate population from all years for each county-group
pop_7988 <- aggregate_pop(miheart$n, margin_time)
# aggregate population from 1980-1984 for each county-group
pop_8084 <- aggregate_pop(miheart$n, margin_time, groups = as.character(1980:1984))
# bind aggregated pop from all years to population data
pop_agg <- aggregate_pop(miheart$n, margin_time, bind_new = TRUE, new_name = "1979-1988")
```
