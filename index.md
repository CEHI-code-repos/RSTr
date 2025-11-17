# RSTr

The Rate Stabilizing Tool R Package (RSTr): Gibbs Samplers for Bayesian
Spatiotemporal CAR Models

The RSTr package is a tool that uses Bayesian spatiotemporal models in
conjunction with C++ to help you quickly and easily generate spatially
smoothed estimates for your spatiotemporal small area data. For the
ArcGIS Python Toolbox, visit the [Rate Stabilizing Toolbox
(RSTbx)](https://github.com/CEHI-code-repos/RateStabilizingToolbox).

## Installation

Windows users must first install
[RTools](https://cran.r-project.org/bin/windows/Rtools/) before
installing any package from source. You may find it useful to use the
`installr` package to install RTools like below:

``` r
install.packages("installr")
installr::install.Rtools()
```

You can install the development version of `RSTr`.

``` r
devtools::install_github("CEHI-code-repos/RSTr", build_vignettes = TRUE)
```

If this is your first time using RSTr, check out
[`vignette("RSTr")`](articles/RSTr.md) to learn how to use the package.

## Thanks!

Thanks for checking out the RSTr package. I hope you find it useful in
your spatiotemporal modeling! Please feel free to give feedback on bugs,
ways to make the package more accessible to use, and features you’d like
to see added to the package.
