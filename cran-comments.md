<!-- version 1.1.3 -->

RSTr has been incremented to verison 1.1.3. Two bugs have been fixed:
- Issues with testing examples that use package 'sf' fixed with requireNamespace()
- After testing on valgrind with rhub, zero 'definitely lost' memory leaks found. I could not find any function calls that trace back to functions in my package, so these are likely from sources outside of the package.

<!-- version 1.1.2 -->
## R CMD check results

0 errors | 0 warnings | 0 notes

RSTr has been incremented to verison 1.1.2. One bug and some small tweaks have been made:
- Fixed segmentation faults caused by sf::st_as_sf()
- Modified temporary file management in example functions
- Missed one potential instance of memory-unsafe calling of Lists in Rcpp code

<!-- version 1.1.1 -->
## R CMD check results

0 errors | 0 warnings | 0 notes

RSTr has been incremented to version 1.1.1. One new feature and a few bugs have been patched to address CRAN issues:

- clang-ASAN, gcc-ASAN, and valgrind exposed memory safety errors from Armadillo objects. I have done my best to rectify potential sources of error in all Rcpp parameter updates by converting calls such as `armatype& object = list["object"];` to `auto object = Rcpp::as<armatype>(list["object"]);`

- a new function has been introduced: split_sample_groups()

- The following functions have been tweaked for optimization, though no arguments or output has changed: standardize_samples(), aggregate_samples(), aggregate_count(), get_medians(), and get_credible_interval()

- Some of the CRAN-only outputs to ensure vignettes don't run model functions had `echo = FALSE` omitted. I've added these in where necessary.

- The tests also pointed out a NOTE that some file directories were particularly large. I have further shrunk the file size of the example models to as small as possible without impeding functionality

<!-- version 1.0.0 -->
## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The previous CRAN submission on 2026/01/08 failed due to errors with vignettes and testing on Windows. I have reduced the scope of testing to help these tests run faster and have also modified the vignettes so that the MCMC models don't run when generating the vignettes for testing with CRAN.
