# Load model

`load_model()` imports an `RSTr` object with name `name` in directory
`dir`.

## Usage

``` r
load_model(name, dir = tempdir())
```

## Arguments

- name:

  The name of the model to load.

- dir:

  The directory in which the model lives.

## Value

An `RSTr` model object.

## Examples

``` r
data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
adj_min <- list(2, 1)
mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE, verbose = FALSE)
mod_mst <- load_model(name = "test", dir = tempdir())
```
