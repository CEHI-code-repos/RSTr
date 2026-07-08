# Changelog

## RSTr 1.2.0

- Added limits to A
- Refactored estimate generation process
- Updated introductory vignette

## RSTr 1.1.4

CRAN release: 2026-01-31

- Increased redundancy of requireNamespace() checks on examples and
  vignettes that use sf and spdep
- Translated update_lambda() to R from C++ to fix memory safety issues

## RSTr 1.1.3

- Issues with testing examples that use package ‘sf’ fixed with
  requireNamespace()
- Fixed valgrind memory leaks directly associated with package cpp
  functions

## RSTr 1.1.2

- Fixed segmentation faults caused by sf::st_as_sf()
- Modified temporary file management in example functions

## RSTr 1.1.1

- Fixed memory safety errors exposed by clang-ASAN, gcc-ASAN, and
  valgrind

## RSTr 1.1.0

- a new function has been introduced: split_sample_groups()

- The following functions have been tweaked for optimization, though no
  arguments or output has changed: standardize_samples(),
  aggregate_samples(), aggregate_count(), get_medians(), and
  get_credible_interval()

## RSTr 1.0.0

CRAN release: 2026-01-17

- Initial CRAN submission.
