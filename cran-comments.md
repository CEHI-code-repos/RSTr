R CMD check --as-cran results:
0 errors | 0 warnings | 0 notes

This version addresses memory safety errors and segfaults in my previous CRAN submission. C++ code that was associated with safety errors has been entirely rewritten in R and taken out of the package.

The package was checked using R CMD check . --as-cran along with different rhub sanitizer builds (clang-asan, clang-ubsan, gcc-asan). All checks completed without errors.

Examples and vignettes that depend on sf and spdep are now conditionally guarded with requireNamespace().