
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cobin: R package for cobin and micobin regression models

<!-- badges: start -->
<!-- badges: end -->

Cobin and micobin regression models are scalable and robust alternative
to beta regression model for continuous proportional data. See the
following paper for more details:

> Lee, C. J., Dahl, B. K., Ovaskainen, O., Dunson, D. B. (2025).
> Scalable and robust regression models for continuous proportional
> data. arXiv preprint arXIV:2504.15269
> <https://arxiv.org/abs/2504.15269>

A dedicated Github repository for reproducing the analysis in the paper
is available at <https://github.com/changwoo-lee/cobin-reproduce>. This
R package repository contains the functions for the cobin and micobin regression
models, as well as sampler for Kolmogorov-Gamma random variables.

Install the package from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("changwoo-lee/cobin")
```

Glossaries: GLM: generalized linear model; GLMM: generalized linear
mixed model; GP: Gaussian process; NNGP: nearest neighbor Gaussian
process; cobin: continuous binomial; micobin: mixture of continuous
binomial;

### vignette

To be added soon!

## Code structure

### Basic functions

- cobin.R:
  - `dcobin(x, theta, lambda)`: Density of
    $\mathrm{cobin}(\theta, \lambda^{-1})$ at x  
  - `rcobin(n, theta, lambda)`: Random variate generation from
    $\mathrm{cobin}(\theta, \lambda^{-1})$
- micobin.R:
  - `dmicobin(x, theta, psi)`: Density of
    $\mathrm{micobin}(\theta, \psi)$ at x  
  - `rmicobin(n, theta, psi)`: Random variate generation from
    $\mathrm{micobin}(\theta, \psi)$

### Cobin / Micobin regression (Bayesian, cobit link)

- cobinreg.R:
  - `cobinreg()`: fit Bayesian cobin GLM or GLMM
  - └── fit_cobin_fixedeffect.R: backend function for cobin GLM
  - └── fit_cobin_mixedeffect.R: backend function for cobin GLMM
- micobinreg.R:
  - `micobinreg()`: fit Bayesian micobin GLM or GLMM
  - └── fit_micobin_fixedeffect.R: backend function for micobin GLM
  - └── fit_micobin_mixedeffect.R: backend function for micobin GLMM

### Spatial cobin / micobin regression (Bayesian, cobit link)

- spcobinreg.R:
  - `cobinreg()`: fit spatial cobin regresson
  - └── fit_cobin_spatial.R: backend function with GP random effect
  - └── fit_cobin_spatial_NNGP.R: backend function with NNGP random
    effect
- spmicobinreg.R:
  - `micobinreg()`: fit Bayesian cobin GLM or GLMM
  - └── fit_micobin_spatial.R: backend function with GP random effect
  - └── fit_micobin_spatial_NNGP.R: backend function with NNGP random
    effect

### Helper functions

- CB.R:
  - `qcb()`, `rcb()`: quantile and random variate generation of
    continuous Bernoulli
- IH.R:
  - `dIH()`: density of Irwin-Hall distribution
- varfunctions.R: collection of functions related to variance function
  of cobin, with numerically stable computation
  - `bft()`: $B(x) = \log((\exp(x)-1)/x)$, cumulant (log partition)
    function
  - `bftprime()`: $B'(x) = 1/(1-\exp(-x))-1/x$, corresponding to inverse
    of cobit link function
  - `bftprimeprime()`, `bftprimeprimeprime()`: $B''(x)$ and $B'''(x)$
  - `bftprimeinv()`: inverse of $B'(x)$, corresponding to cobit link
    function
  - `Vft()`: $B''((B')^{-1}(\mu))$, variance function of cobin
- cobinfamily.R:
  - `cobinfamily()`: a list of functions and expressions needed to fit
    cobin GLM

### cobin regression (non-Bayesian)

- glm.cobin.R:
  - `glm.cobin()`: fit cobin GLM using iteratively reweighted least
    squares (`stats::glm.fit`). Supports link functions “cobit”,
    “logit”, “probit”, “cloglog”, “cauchit”.
