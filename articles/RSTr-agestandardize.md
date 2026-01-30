# 04: Generating Estimates: Age-standardization

## Overview

In the previous vignette, we discussed the model setup process in-depth.
But how do we get our estimates once we’ve run our model? In this
vignette, we discuss extracting estimates from our model object with the
[`get_estimates()`](../reference/get_estimates.md) function, and how to
age-standardize those estimates with
[`age_standardize()`](../reference/age_standardize.md).

## The `get_estimates()` function

In the RSTr introductory vignette, we generated age-standardized
estimates for `lambda` based on our example Michigan dataset. To extract
rates from an `RSTr` object, we can simply run
[`get_estimates()`](../reference/get_estimates.md):

``` r
mod_mst <- mstcar(name = "my_test_model", data = miheart, adjacency = miadj)
```

![](RSTr-agestandardize_files/figure-html/unnamed-chunk-1-1.png)

``` r
estimates <- get_estimates(mod_mst, rates_per = 1e5)
head(estimates)
#>   county group year  medians ci_lower  ci_upper  rel_prec events population
#> 1  26001 35-44 1979 35.46173 24.16108  54.96574 1.1511806      1        964
#> 2  26003 35-44 1979 60.38549 43.35458 106.10903 0.9622502      1       1011
#> 3  26005 35-44 1979 19.97558 13.68447  29.12392 1.2938009      0       9110
#> 4  26007 35-44 1979 24.24089 16.48251  41.01521 0.9881057      0       3650
#> 5  26009 35-44 1979 28.19778 14.52846  42.56563 1.0057284      0       1763
#> 6  26011 35-44 1979 38.03138 27.70294  52.97160 1.5050813      0       1470
```

## The `age_standardization()` function

In many cases, we will want to age-standardize our estimates based on
some (or all) age groups in our dataset. In our Michigan dataset, we
have six ten-year age groups over which we can standardize; let’s
age-standardize from ages 35-64. For `RSTr` objects,
[`age_standardize()`](../reference/age_standardize.md) takes in four
arguments:

- `RSTr_obj`: The `RSTr` model object created with `*car()`;

- `std_pop`: A `vector` of standard populations associated with the age
  groups of interest. Since our Michigan data is from 1979-1988, we can
  use 1980 standard populations from
  [NIH](https://seer.cancer.gov/stdpopulations/stdpop.19ages.html). It
  is recommended that you use the standard population that is most
  closely associated with your dataset;

- `new_name`: The name of your new standard population group; and

- `groups`: A `vector` of names matching each group of interest. To
  age-standardize by all groups in a dataset, leave this argument blank.

Once we have our `std_pop` vector, we can age-standardize our estimates:

``` r
std_pop <- c(113154, 100640, 95799)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "35-64", groups = c("35-44", "45-54", "55-64"))
mod_mst
#> RSTr object:
#> 
#> Model name: my_test_model 
#> Model type: MSTCAR 
#> Data likelihood: binomial 
#> Estimate Credible Interval: 95% 
#> Number of geographic units: 83 
#> Number of samples: 6000 
#> Estimates age-standardized: Yes 
#> Age-standardized groups: 35-64 
#> Estimates suppressed: No
```

Notice now that the `mod_mst` object indicates we have age-standardized
our estimates and the names of our age-standardized group. We can also
add on to our list of age-standardized estimates by simply specifying a
different group:

``` r
std_pop <- c(68775, 34116, 9888)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "65up", groups = c("65-74", "75-84", "85+"))
mod_mst
#> RSTr object:
#> 
#> Model name: my_test_model 
#> Model type: MSTCAR 
#> Data likelihood: binomial 
#> Estimate Credible Interval: 95% 
#> Number of geographic units: 83 
#> Number of samples: 6000 
#> Estimates age-standardized: Yes 
#> Age-standardized groups: 35-64 65up 
#> Estimates suppressed: No
```

If we want to generate estimates for *all* groups, i.e. 35 and up, we
can omit the `groups` argument and expand `std_pop` to include all of
our populations:

``` r
std_pop <- c(113154, 100640, 95799, 68775, 34116, 9888)
mod_mst <- age_standardize(mod_mst, std_pop, new_name = "35up")
mod_mst
#> RSTr object:
#> 
#> Model name: my_test_model 
#> Model type: MSTCAR 
#> Data likelihood: binomial 
#> Estimate Credible Interval: 95% 
#> Number of geographic units: 83 
#> Number of samples: 6000 
#> Estimates age-standardized: Yes 
#> Age-standardized groups: 35-64 65up 35up 
#> Estimates suppressed: No
mst_estimates_as <- get_estimates(mod_mst)
head(mst_estimates_as)
#>   county group year  medians ci_lower ci_upper rel_prec events population
#> 1  26001 35-64 1979 177.1760 147.9660 209.6233 2.873563      7       3353
#> 2  26003 35-64 1979 264.3501 212.0934 323.4194 2.374558     12       3105
#> 3  26005 35-64 1979 126.8562 111.1458 148.1560 3.427600     27      23926
#> 4  26007 35-64 1979 144.9922 118.1114 179.2119 2.373010     15      10000
#> 5  26009 35-64 1979 144.1579 111.8566 185.7918 1.949788     11       5152
#> 6  26011 35-64 1979 196.5609 160.2600 235.3754 2.616786      8       4517
```

Now, `get_estimates(mod_mst)` shows the age-standardized estimates as
opposed to our non-standardized estimates. Should you want to see the
non-standardized estimates instead, you can set the argument
`standardized = FALSE`.

## Final thoughts

In this vignette, we explored the
[`get_estimates()`](../reference/get_estimates.md) function and
investigated age-standardization with the
[`age_standardize()`](../reference/age_standardize.md) function.
Age-standardization is one of the most important features of the RSTr
package; using just a few arguments, we can easily generate estimates
across our population groups.
