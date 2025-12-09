# 04: Model Informativeness

## Overview

- RSTr’s models help to smooth estimates with low event counts
- What if events are too low?
- Standard BYM models oversmooth in these regions
- Need to restrict the model informativeness
- Discuss how informativeness works
- How do estimates differ between unrestricted/restricted models
- When should you use one model over the other
- Developments in restricted models

## Section

Vignette in development.

Informativeness:

The goal of informativeness is to limit/prevent a situation in which the
model oversmooths the data. Oversmoothing can lead to masking/hiding
real geographic disparities. A map with data that is oversmoothed has
overly distinct geographic patterns \[show an example\].

Informativeness refers to the amount/degree of spatial smoothing that
the model allows. As a default, we suggest an informativeness value of
6. This is based upon.

However, there are situations in which you may want to use a different
level of informativeness, for example:

## Final thoughts

Vignette in development.
