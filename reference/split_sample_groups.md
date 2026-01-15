# Split sample groups

Sequesters stratified sociodemographic group margin into individual
array margins.

## Usage

``` r
split_sample_groups(sample, new_groups, delimiter = "_")
```

## Arguments

- sample:

  an `array` of samples imported with
  [`load_samples()`](load_samples.md)

- new_groups:

  A string vector of names for each new group.

- delimiter:

  A character that specifies the break between group categories.

## Value

An `array` of samples with separate margins for stratified groups.

## Details

When using [`aggregate_samples()`](aggregate_samples.md) or
[`standardize_samples()`](standardize_samples.md), the group/age margin
must only feature groups of similar type, e.g., you cannot
age-standardize with groups that specify both age and race.
`split_sample_groups()` sequesters each category of group into its own
margin to allow group-aggregation and age-standardization of these
multiply-stratified groups. Ensure that the delimiter character is only
used to split groups. E.g., for an age-sex group named `35-64_m`, `"_"`
will split the margins with names `"35-64"` and `"m"`, whereas for a
group named `35_64_m`, `split_sample_group()` will fail.

## Examples

``` r
dimnames(minsplit)[2] # Can't age-standardize due to age-sex stratification
#> $group
#> [1] "35-44_f" "45-54_f" "55-64_f" "35-44_m" "45-54_m" "55-64_m"
#> 
new_groups = c("age", "sex")
delimiter = "_"
sample_split <- split_sample_groups(minsplit, new_groups, delimiter)
#> Error in split_sample_groups(minsplit, new_groups, delimiter): object 'samples_new_dims' not found
dimnames(sample_split)[2:3] # can now age-standardize
#> Error: object 'sample_split' not found
std_pop <- c(113154, 100640, 95799)
age_margin <- 2
sample_as <- standardize_samples(sample_split, std_pop, age_margin)
#> Error: object 'sample_split' not found
```
