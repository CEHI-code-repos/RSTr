## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The previous CRAN submission on 2026/01/08 failed due to errors with vignettes and testing on Windows. I have reduced the scope of testing to help these tests run faster and have also modified the vignettes so that the MCMC models don't run when generating the vignettes for testing with CRAN.