# Rate Stabilizing Tool: Gibbs Samplers for Bayesian Spatiotemporal CAR Models

Takes Poisson or Binomial discrete spatial data and runs a Gibbs sampler
for a variety of Spatiotemporal Conditional Autoregressive (CAR) models.
Includes measures to prevent estimate over-smoothing through a
restriction of model informativeness for select models. Also provides
tools to load output and get median estimates. Implements methods from
Besag, York, and Mollié (1991) "Bayesian image restoration, with two
applications in spatial statistics" \doi{10.1007/BF00116466}, Gelfand
and Vounatsou (2003) "Proper multivariate conditional autoregressive
models for spatial data analysis" doi:10.1093/biostatistics/4.1.11,
Quick et al. (2017) "Multivariate spatiotemporal modeling of
age-specific stroke mortality" doi:10.1214/17-AOAS1068, and Quick et al.
(2021) "Evaluating the informativeness of the Besag-York-Mollié CAR
model" doi:10.1016/j.sste.2021.100420.

## Details

The RSTr package uses Bayesian spatiotemporal modeling to spatially
smooths discrete small-area event rates using information from
neighboring spatial regions. See \`browseVignettes("RSTr")\` for a
series of tutorials on basic usage of the RSTr functions.

## Author

David DeLara \[aut, cre\] (ORCID:
\<https://orcid.org/0000-0003-0485-7549\>), Centers for Disease Control
and Prevention \[aut\] (https://ror.org/042twtr12)

Maintainer: David DeLara \<sfq1@cdc.gov\>

## References

Besag, J., York, J., and Mollié, A. (1991). Bayesian Image Restoration
with Two Applications in Spatial Statistics (with Discussion). Annals of
the Institute of Statistical Mathematics, 43, 1–59.
[doi:10.1007/BF00116466](https://doi.org/10.1007/BF00116466)

Gelfand, A. E., & Vounatsou, P. (2003). Proper multivariate conditional
autoregressive models for spatial data analysis. Biostatistics, 4(1),
11–25.
[doi:10.1093/biostatistics/4.1.11](https://doi.org/10.1093/biostatistics/4.1.11)

Quick, et al. (2017). Multivariate spatiotemporal modeling of
age-specific stroke mortality. Annals of Applied Statistics, 11(4),
2165–2177.
[doi:10.1214/17-AOAS1068](https://doi.org/10.1214/17-AOAS1068)

Quick, et al. (2021). Evaluating the informativeness of the
Besag-York-Mollié CAR model. Spatial and Spatio-temporal Epidemiology,
37, 100420.
[doi:10.1016/j.sste.2021.100420](https://doi.org/10.1016/j.sste.2021.100420)
