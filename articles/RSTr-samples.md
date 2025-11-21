# 04: Gathering and Loading Samples

## Overview

After initializing the model with `initialize_*()`, you can gather and
load samples using the [`run_sampler()`](../reference/run_sampler.md)
and [`load_samples()`](../reference/load_samples.md) functions. Finally,
using some basic plotting, we can look at our results and make sure they
make sense.

## The `run_sampler()` function

[`run_sampler()`](../reference/run_sampler.md) takes in four arguments:
`name`, `dir`, `show_plots`, and `discard_burnin`. `name` and `dir` are
associated with the name of the folder that contains the model content
and the directory where that folder lives, respectively. `show_plots`
allows you to hide the traceplots generated while the model runs and
`discard_burnin` prevents samples from being saved before the
2000-iteration burn-in period in case the dataset is particularly large.
The information needed to load the model is already saved in `name`, and
so [`run_sampler()`](../reference/run_sampler.md) pulls this information
in to generate parameter samples. When running
[`run_sampler()`](../reference/run_sampler.md), the R console will
generate outputs describing the batch number, the total iteration
number, and the time when the last batch started. By default,
[`run_sampler()`](../reference/run_sampler.md) will generate 6000
samples, but other values can be specified as desired, rounded down to
the nearest 100. For example, we can begin a simple model and generate
samples:

``` r
library(RSTr)
initialize_mstcar(name = "my_test_model", data = miheart, adjacency = miadj)
#> Checking data...
```

![](RSTr-samples_files/figure-html/unnamed-chunk-1-1.png)

    #> Warning in initialize_model(name = name, data = data, adjacency = adjacency, :
    #> Seed is not set using `seed` arg in `initialize_*()`; samples may not be
    #> replicable.
    #> Checking spatial data...
    #> Checking initial_values...
    #> The following objects were created using defaults in 'initial_values': beta lambda Z tau2 G rho Ag
    #> Checking priors...
    #> The following objects were created using defaults in 'priors': lambda_sd tau_a tau_b G_df G_scale Ag_scale Ag_df rho_a rho_b rho_sd
    #> Model ready!

``` r
run_sampler("my_test_model", dir = tempdir())
#> Starting sampler on Batch 1 at Fri Nov 21 15:15:57
```

![](RSTr-samples_files/figure-html/unnamed-chunk-2-1.png)

    #> Model finished at Fri Nov 21 15:16:22

## The `load_samples()` function

Once [`run_sampler()`](../reference/run_sampler.md) tells you that the
model is finished running, you can import samples into R using
[`load_samples()`](../reference/load_samples.md).
[`load_samples()`](../reference/load_samples.md) takes in four
arguments:

- `name`: Name of the model;

- `dir`: Directory of the model;

- `param`: The parameter to import samples for; and

- `burn`: Specifies a burn-in period for samples. This allows the model
  time to stabilize before using samples to generate estimates. By
  default, has a burn-in period of 2000 samples.

Any `dimnames` that were saved to `data` will be applied to the samples
as appropriate. Here, we pull in the `lambda` samples for our test
Michigan dataset with a 2000-sample burn-in period (as specified by the
default arguments). We also multiply by 100,000 as it is common to
display mortality rates per 100,000 individuals:

``` r
samples <- load_samples("my_test_model") * 1e5
```

### Group aggregation

In many cases, we will want to aggregate our data across non-age groups,
such as when looking at prevalence estimates or to simply consolidate
other sociodemographic groups. In our Michigan dataset, we have 10 years
of data over which we can consolidate to look at prevalence. In these
cases, we need to pull in the population array as a weight. First, we
need to check which margin contains our year information using the
[`dim()`](https://rdrr.io/r/base/dim.html) function:

``` r
dim(samples)
#> [1]  83   6  10 400
```

Our `samples` array has four margins with dimensions `83`, `6`, `10`,
and `400`, representing the spatial regions, age groups, time periods,
and iterations, respectively. Let’s set a variable `margin_time` to
represent our time period margin and aggregate our `samples` estimates
across 1979-1988 using the
[`aggregate_groups()`](../reference/aggregate_groups.md) function. The
population data needed to weight our samples can be pulled from the
model folder using the [`load_pop()`](../reference/load_pop.md)
function:

``` r
margin_time <- 3
pop <- load_pop("my_test_model")
samples_7988 <- aggregate_groups(samples, pop, margin_time)
```

Now, we have a standalone sample array for our 1979-1988 samples. But
what if we are interested in both the individual year data *and* the
prevalence data? We can alternatively bind these new samples to our main
`samples` array by adding in values for the `bind_new` and `new_name`
arguments:

``` r
samples <- aggregate_groups(samples, pop, margin_time, bind_new = TRUE, new_name = "1979-1988")
```

### Age-standardization

The process of age-standardization is similar to that of
group-aggregation, but requires a bit more nuance in its use. In our
Michigan dataset, we have six 10-year age groups that start at age 35
years. We can age-standardize these into a new 35-64 group. Since we are
using data from 1979-1988, we can use 1980 standard populations from
[NIH](https://seer.cancer.gov/stdpopulations/stdpop.19ages.html) to
generate a `std_pop` vector:

``` r
age <- c("35-44", "45-54", "55-64")
std_pop <- c(113154, 100640, 95799)
names(std_pop) <- age
```

With `std_pop` generated, we need to then check which margin contains
our age group information using the
[`dim()`](https://rdrr.io/r/base/dim.html) function:

``` r
dim(samples)
#> [1]  83   6  11 400
```

Our age groups lay along the second margin. Let’s set a variable
`margin_age` and standardize our `samples` estimates across ages 35-64
using the [`age_standardize()`](../reference/age_standardize.md)
function:

``` r
margin_age <- 2
samples_3564 <- age_standardize(samples, std_pop, margin_age, groups = c("35-44", "45-54", "55-64"))
```

Note that there may be times where you have groups stratified by both
age and other sociodemographic groups. In these cases, you’ll have to
refactor your sample array so that the age groups are separated from
your other groups before doing age-standardization, such as using
[`as.data.frame.table()`](https://rdrr.io/r/base/table.html) in
conjunction with [`xtabs()`](https://rdrr.io/r/stats/xtabs.html).

Now, we have a standalone array for our age-standardized 35-64 age
group. Similarly to
[`aggregate_groups()`](../reference/aggregate_groups.md), we can
alternatively consolidate this into our main `samples` array by adding
in values for the `bind_new` and `new_name` arguments:

``` r
samples <- age_standardize(
    samples,
    std_pop,
    margin_age,
    groups = c("35-44", "45-54", "55-64"),
    bind_new = TRUE,
    new_name = "35-64"
)
```

Now, our samples for `samples` are aggregated by year and
age-standardized, and we have matching array values for `pop`. Note that
if you plan on doing a mix of non-age aggregation and
age-standardization, do age-standardization *after* aggregation, as
doing age-standardization first will alter the results of any
aggregation done afterward.

### Traceplots

Before we generate the median estimates for our samples, let’s quickly
check some high- and low-population counties to see the stability of the
samples:

``` r
counties <- c(
  which.max(miheart$n[, "55-64", "1979"]),
  which.min(miheart$n[, "55-64", "1979"])
)
x <- dimnames(samples)[[4]]
yhigh <- samples[counties[1], "35-64", "1979", ]
ylow <- samples[counties[2], "35-64", "1979", ]
raw_high <- sum(miheart$Y[counties[1], 1:3, "1979"]) /
  sum(miheart$n[counties[1], 1:3, "1979"]) *
  1e5
raw_low <- sum(miheart$Y[counties[2], 1:3, "1979"]) /
  sum(miheart$n[counties[2], 1:3, "1979"]) *
  1e5
par(mfrow = c(1, 2))
title_cty1 <- paste0("Smoothed rates, 1979, 35-64,\nFIPS ", names(counties[1]))
plot(x, yhigh, type = "l", main = title_cty1)
abline(h = raw_high, col = "red")
title_cty2 <- paste0("Smoothed rates, 1979, 35-64,\nFIPS ", names(counties[2]))
plot(x, ylow, type = "l", main = title_cty2)
abline(h = raw_low, col = "red")
```

![](RSTr-samples_files/figure-html/unnamed-chunk-12-1.png)

This code generates two traceplots representing our sample values over
time for a given county-group-year. On the left, we can see a traceplot
for the highest-population county, and on the right is a traceplot for
the lowest-population county. Each traceplot also includes a red line
representing the mean event rate for that county-group-year. Note that
in both traceplots, the mean line is nearby the plots, but higher,
indicating that the rates in both of these counties were attenuated
thanks to spatial smoothing.

The traceplot on the left is exactly what we want to see: it is clearly
fluctuating around a certain value and will give us reliable estimates
for that county. The right traceplot, however, is a bit less favorable:
the value doesn’t seem to want to stabilize and jumps between values
over the course of the model. However, the rate itself is naturally
high, so the intensity of the fluctuation isn’t shocking. Smaller
counties like the one shown here demonstrate the limits of `RSTr`: while
the samples do hover around a single value for some iterations, the
estimated values we would get for the county-group-year on the right
will not be as reliable as estimates on the left due to the large
variability of the samples gathered. To learn more about reliability
measures, read
[`vignette("RSTr-medians")`](../articles/RSTr-medians.md).

## Closing Thoughts

In this vignette, we discussed generating samples with
[`run_sampler()`](../reference/run_sampler.md), importing those samples
into R, age-standardization using
[`load_samples()`](../reference/load_samples.md), generating traceplots
to get a gut-check on our dataset, along with the implications of
sporadic traceplots.
