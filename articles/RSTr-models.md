# Appendix A: The CAR Hierarchical Models

## Overview

In this vignette, we outline the hierarchical models used in the `RSTr`
package, along with the full-conditional distributions used for each
update.

## The UCAR Hierarchical Model

The UCAR model used by `RSTr` is based on the model developed by [Besag,
York, and Mollié (1991)](https://doi.org/10.1007/BF00116466) with
modifications using inverse transform sampling for restricted
informativeness based on [Quick, et
al. (2021)](https://doi.org/10.1016/j.sste.2021.100420):

For models using `method = "binomial"`,

\\ \begin{split} Y\_{i} &\sim \text{Binomial}(n\_{i}, \lambda\_{i}) \\
\theta\_{i} &= \text{Logit}(\lambda\_{i}) \\ \end{split} \\

For models using `method = "poisson"`, \\ \begin{split} Y\_{i} &\sim
\text{Poisson}(n\_{i} \lambda\_{i}) \\ \theta\_{i} &=
\text{Log}(\lambda\_{i}) \\ \end{split} \\

For both models,

\\ \begin{split} \theta\_{i} &\sim \text{Normal}(\beta\_{j} + Z\_{i},
\tau^2), \\ i &=\\1,...,N\_{s}\\,\\ j =\\1,...,N\_{is}\\ \\
p(\beta\_{j}) &\propto 1 \\ Z &\sim \text{CAR}(\sigma^2) \\ \sigma^2
&\sim \text{InvGamma}(a\_\sigma,b\_\sigma) \\ \tau^2 &\sim
\text{InvGamma}(a\_\tau,b\_\tau) \end{split} \\

## The MCAR Hierarchical Model

The MCAR model used by `RSTr` is based on the model developed by
[Gelfand and Vounatsou
(2003)](https://doi.org/10.1093/biostatistics/4.1.11):

For models using `method = "binomial"`,

\\ \begin{split} Y\_{ik} &\sim \text{Binomial}(n\_{ik}, \lambda\_{ik})
\\ \theta\_{ik} &= \text{Logit}(\lambda\_{ik}) \\ \end{split} \\

For models using `method = "poisson"`,

\\ \begin{split} Y\_{ik} &\sim \text{Poisson}(n\_{ik}, \lambda\_{ik}) \\
\theta\_{ik} &= \text{Log}(\lambda\_{ik}) \\ \end{split} \\

For both models,

\\ \begin{split} \theta\_{ik} &\sim \text{Normal}(\beta\_{jk} + Z\_{ik},
\tau_k^2), \\ i &=\\1,...,N_s\\, k =\\1,...,N\_{g}\\,
j=\\1,...,N\_{is}\\ \\ p(\beta\_{jk}) &\propto 1 \\ Z &\sim
\text{CAR}(G) \\ G &\sim \text{InvWishart}(\nu,G_0) \\ \tau^2 &\sim
\text{InvGamma}(a\_\tau,b\_\tau) \end{split} \\

## The MSTCAR Hierarchical Model

The MSTCAR model used by `RSTr` is based on the model developed by
[Quick, et al. (2017)](https://doi.org/10.1214/17-AOAS1068):

For models using `method = "binomial"`,

\\ \begin{split} Y\_{ikt} &\sim \text{Binomial}(n\_{ikt},
\lambda\_{ikt}) \\ \theta\_{ikt} &= \text{Logit}(\lambda\_{ikt}) \\
\end{split} \\

For models using `method = "poisson"`,

\\ \begin{split} Y\_{ikt} &\sim \text{Poisson}(n\_{ikt} \lambda\_{ikt})
\\ \theta\_{ikt} &= \text{Log}(\lambda\_{ikt}) \\ \end{split} \\

For both models,

\\ \begin{split} \theta\_{ikt} &\sim \text{Normal}(\beta\_{jkt} +
Z\_{ikt}, \tau_k^2), \\ i &=\\1,...,N_s\\,\\ k =\\1,...,N_g\\,\\
t=\\1,...,N_t\\,\\ j=\\1,...,N\_{is}\\ \\ p(\beta\_{j}) &\propto 1 \\ Z
&\sim \text{MSTCAR}(\mathcal{G}, \mathcal{R}), \\
\mathcal{G}=\\G_1,...,G\_{N_t}\\, \\ \mathcal{R}=\\R_1,...,R\_{N_g}\\ \\
G_t &\sim \text{InvWishart}(A_G, \nu) \\ A_G &\sim
\text{Wishart}(A\_{G_0}, \nu_0) \\ R_k &= \text{AR}(1,\rho_k) \\ \rho_k
&\sim \text{Beta}(a\_{\rho}, b\_{\rho}) \\ \tau_k^2 &\sim
\text{InvGamma}(a\_\tau,b\_\tau) \end{split} \\
