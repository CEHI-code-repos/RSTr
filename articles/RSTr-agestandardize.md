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
#> Starting sampler on Batch 1 at Thu Dec 18 20:47:05
```

![](RSTr-agestandardize_files/figure-html/unnamed-chunk-2-1.png)

    #> Generating estimates...
    #> Model finished at Thu Dec 18 20:47:31

``` r
estimates <- get_estimates(mod_mst, rates_per = 1e5)
head(estimates)
#>   county group year  medians ci_lower ci_upper rel_prec events population
#> 1  26001 35-44 1979 45.07149 31.48676 59.17320 1.627926      1        964
#> 2  26003 35-44 1979 56.47176 35.01048 85.83469 1.111119      1       1011
#> 3  26005 35-44 1979 24.35092 17.47599 30.61988 1.852641      0       9110
#> 4  26007 35-44 1979 35.13386 18.28414 48.08389 1.178999      0       3650
#> 5  26009 35-44 1979 34.24613 19.40691 45.97550 1.288970      0       1763
#> 6  26011 35-44 1979 46.13104 28.16870 65.56977 1.233415      0       1470
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
#> 1  26001 35-64 1979 203.3206 157.7320 241.9712 2.413612      7       3353
#> 2  26003 35-64 1979 290.8983 223.9136 384.6817 1.809428     12       3105
#> 3  26005 35-64 1979 124.6186 103.7021 142.7059 3.195035     27      23926
#> 4  26007 35-64 1979 174.7536 138.9969 204.4177 2.671222     15      10000
#> 5  26009 35-64 1979 171.6994 134.1164 209.4447 2.279347     11       5152
#> 6  26011 35-64 1979 223.7704 183.4417 276.1411 2.413933      8       4517
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
