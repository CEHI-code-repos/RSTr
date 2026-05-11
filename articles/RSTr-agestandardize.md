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
#>   county group year  medians ci_lower  ci_upper rel_prec events population
#> 1  26001 35-44 1979 36.36755 24.48474  48.71868 1.500687      1        964
#> 2  26003 35-44 1979 88.28406 60.89756 154.96037 0.938565      1       1011
#> 3  26005 35-44 1979 22.40370 14.71798  32.14123 1.285851      0       9110
#> 4  26007 35-44 1979 31.01046 21.02774  43.53680 1.377688      0       3650
#> 5  26009 35-44 1979 32.54242 19.32274  42.02634 1.433359      0       1763
#> 6  26011 35-44 1979 46.32233 31.63450  63.83114 1.438732      0       1470
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
#> 1  26001 35-64 1979 168.8406 140.0636 198.8595 2.871638      7       3353
#> 2  26003 35-64 1979 281.6352 237.2297 373.8780 2.061021     12       3105
#> 3  26005 35-64 1979 122.7245 103.8202 146.9175 2.847616     27      23926
#> 4  26007 35-64 1979 152.3309 118.8884 182.3210 2.401460     15      10000
#> 5  26009 35-64 1979 157.3300 113.3391 184.9485 2.197056     11       5152
#> 6  26011 35-64 1979 211.4462 169.1618 260.8167 2.306983      8       4517
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
