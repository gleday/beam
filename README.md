<!-- header = image + name -->
<h1 align="center">
  <br>
  <img src="https://github.com/gleday/gleday/blob/main/images/inferred_graph.png?raw=true" alt="beam" width="200">
  <br>
  beam
  <br>
</h1>

<!-- headline -->
<h4 align="center">Fast Bayesian inference of network structures.</h4>

<!-- badges: start -->
<p align="center">
  <a href="https://CRAN.R-project.org/package=beam"><img src="https://www.r-pkg.org/badges/version/beam"></a>
  <a href="https://cran.r-project.org/web/checks/check_results_beam.html"><img src="https://badges.cranchecks.info/summary/beam.svg"></a>
  <a href="https://app.codecov.io/github/gleday/beam?branch=master"><img src="https://codecov.io/gh/gleday/beam/branch/master/graph/badge.svg"></a>
  <a href="https://cran.r-project.org/package=beam"><img src="http://cranlogs.r-pkg.org/badges/grand-total/beam?color=#1F65CC"></a>
  <a href="https://cran.r-project.org/package=beam"><img src="http://cranlogs.r-pkg.org/badges/last-month/beam?color=#4197D9"></a>
</p>
<!-- badges: end -->

## Features

* inference of conditional independence structures
* inference of marginal independence structures
* computationally efficient (no MCMC) 
* memory-efficient
* able to address problems with thousands of variables on
a standard laptop in just a few seconds
* outperforms popular Bayesian and non-Bayesian methods

## Installation

To install **beam** from R:

```R
# Install/load R package devtools
install.packages("devtools")
library(devtools)

# Install/load R package beam from github
install_github("gleday/beam")
library(beam)
```

## Citation

This R package implements the method described in

Leday, G.G.R. and Richardson, S. (2019).
[Fast Bayesian inference in large Gaussian graphical models](https://doi.org/10.1111/biom.13064).
*Biometrics.* 75(4), 1288--1298.
