# RSTr 1.1.2

- Fixed segmentation faults caused by sf::st_as_sf()
- Modified temporary file management in example functions

# RSTr 1.1.1

- Fixed memory safety errors exposed by clang-ASAN, gcc-ASAN, and valgrind

# RSTr 1.1.0

- a new function has been introduced: split_sample_groups()

- The following functions have been tweaked for optimization, though no arguments or output has changed: standardize_samples(), aggregate_samples(), aggregate_count(), get_medians(), and get_credible_interval()

# RSTr 1.0.0

* Initial CRAN submission.
