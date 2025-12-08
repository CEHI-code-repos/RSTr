# Appendix B: Troubleshooting

## Overview

The `*car()` functions include checks to ensure that all data is set up
in a way that the model can use. If the checks notice something
incompatible with the model, an error will be displayed in the console.
Here, we list each possible error for the MSTCAR model, including its
reasoning and ways to check and fix your inputs.

## Generic errors and warnings

All checks that are run on the data check for general information, such
as extraneous list objects which cannot be used by the model. In these
cases, `RSTr` will not use these parts of the dataset. It is recommended
to remove these datasets from the model to prevent unnecessary clutter.

## Errors and warnings with list `data`

- `One or more objects missing from list 'data'`: The `data` object
  requires a list with two objects: an array `Y` for mortality counts
  and an array `n` for population counts. These are needed to perform
  parameter updates inside of the MSTCAR model. Check `names(data)` to
  see a vector of names that exists inside of your array and fix names
  accordingly.

- `Data not same dimensions. Ensure dim(Y) == dim(n)`: This error occurs
  when the arrays fed into the `data` list are not of the same size. Use
  the [`dim()`](https://rdrr.io/r/base/dim.html) function to check that
  the dimensions of your two arrays do line up. You can use the
  [`aperm()`](https://rdrr.io/r/base/aperm.html) function to rearrange
  arrays if they have the same set of dimensions on different margins.

- `Invalid Y/n values. Check that all Y's/n's are at least 0 and finite`:
  While `NA` and `NULL` are accepted values in `Y`, both `Y` and `n`
  must be comprised of positive finite values. Check your data with
  [`summary()`](https://rdrr.io/r/base/summary.html) for a quick look at
  your minimum values.

## Errors and warnings with list `spatial_data`

- `Adjacency different length than data. Ensure length(adjacency) == dim(Y)[1]`:
  Adjacency information must be present for all regions. Look at the
  length of your adjacency information: it should line up with the
  number of rows in your `Y` and `n` arrays.

- `Some regions have no neighbors. Ensure all regions have at least one neighbor`:
  If a region does not have a neighbor, the model will not run. If
  working with `list` adjacency information, you can first check which
  regions have no neighbors by running
  `which(lapply(adjacency, \(x) all(x == 0)))` and then fixing the
  regions accordingly. For assistance with fixing the adjacency
  information, read
  [`vignette("RSTr-adjacency")`](../articles/RSTr-adjacency.md).

## `initial_values` and `priors`

For troubleshooting inputs with `initial_values` and `priors`, reference
[`vignette("RSTr-initialvalues")`](../articles/RSTr-initialvalues.md)
and [`vignette("RSTr-priors")`](../articles/RSTr-priors.md),
respectively.

## `ignore_checks`

If you have input data that you are certain is correct, yet are
receiving error messages anyway, set the `ignore_checks` argument to
`TRUE` when setting up your model to skip data checks. This is *not*
recommended as it is easy to break `RSTr` when data isn’t properly
checked, and should only be used as a last resort when using `RSTr`.
