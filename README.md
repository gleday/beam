
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/beam)](https://CRAN.R-project.org/package=beam)
[![CRAN checks](https://badges.cranchecks.info/summary/beam.svg)](https://cran.r-project.org/web/checks/check_results_beam.html)
[![Coverage status](https://codecov.io/gh/gleday/beam/branch/master/graph/badge.svg)](https://app.codecov.io/github/gleday/beam?branch=master)
[![](http://cranlogs.r-pkg.org/badges/grand-total/beam?color=#1F65CC)](https://cran.r-project.org/package=beam)
[![](http://cranlogs.r-pkg.org/badges/last-month/beam?color=#4197D9)](https://cran.r-project.org/package=beam)
<!-- badges: end -->

# beam

This R package implements the method described in

Leday, G.G.R. and Richardson, S. (2019).
[Fast Bayesian inference in large Gaussian graphical models](https://doi.org/10.1111/biom.13064).
*Biometrics.* 75(4), 1288--1298.

## Description

beam:

* enables the reconstruction of marginal and conditional independence
graphs from high-dimensional data
* is computationally very efficient, addressing problems with
thousands of variables in just a few seconds.
* outperforms popular Bayesian and non-Bayesian methods

## Installation

If you wish to install **beam** from R:

```R
# Install/load R package devtools
install.packages("devtools")
library(devtools)

# Install/load R package beam from github
install_github("gleday/beam")
library(beam)
```

