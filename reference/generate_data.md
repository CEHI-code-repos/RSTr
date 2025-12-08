# Generate count data for RSTr object

`generate_data()` converts a long `table` featuring event counts across
regions and other optional margins into a `list` that is readable by
`*car()`.

## Usage

``` r
generate_data(table, event, population, region, group = NULL, time = NULL)
```

## Arguments

- table:

  A `table` containing event and mortality counts stratified by
  group/region/time.

- event:

  The column containing event counts.

- population:

  The column containing population counts.

- region:

  The column containing region names.

- group:

  An optional column containing sociodemographic group names.

- time:

  An optional column containing time period names.

## Value

A `list` of mortality and population counts organized into
multi-dimensional arrays.

## Details

`generate_data()` will sum along any group/time stratifications that
aren't specified; for example, if your dataset contains time periods and
time is not specified in `generate_data()`, the output will be a sum of
all time periods. Filter data by desired groups and time periods before
running `generate_data()`.

## Examples

``` r
ma_data <- maexample[!is.na(maexample$Year), ]
# Generates data from 1979-1981 stratified by sex
ma_data_mst <- generate_data(ma_data, Deaths, Population, County.Code, Sex.Code, Year.Code)
ma_data_79 <- ma_data[ma_data$Year == 1979, ]
# Generates 1979 data stratified by sex
ma_data_m <- generate_data(ma_data_79, Deaths, Population, County.Code, Sex.Code)
# Generates 1979 data summarized for all sexes
ma_data_u <- generate_data(ma_data_79, Deaths, Population, County.Code)
```
