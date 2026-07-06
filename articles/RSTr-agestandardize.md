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
#>   county group year  medians ci_lower ci_upper rel_prec events population
#> 1  26001 35-44 1979 36.36755 24.48474 48.71868 1.500687      1        964
#> 2  26001 35-44 1980 33.96586 25.86092 44.63596 1.809096      1        995
#> 3  26001 35-44 1981 30.30507 23.45830 40.78845 1.748690      0        988
#> 4  26001 35-44 1982 29.32680 20.88678 40.84538 1.469382      1       1018
#> 5  26001 35-44 1983 26.57939 18.66534 35.69001 1.561228      2       1039
#> 6  26001 35-44 1984 21.74369 17.00168 31.67577 1.481775      0       1058
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
#> 2  26001 35-64 1980 158.6862 135.2334 184.0058 3.253610     12       3421
#> 3  26001 35-64 1981 143.1623 123.4511 168.0858 3.207417      7       3431
#> 4  26001 35-64 1982 142.3206 121.1283 170.4376 2.886282      3       3455
#> 5  26001 35-64 1983 138.2588 118.3071 164.8520 2.970443     11       3478
#> 6  26001 35-64 1984 119.4551 104.3687 142.4879 3.133723      5       3519
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
