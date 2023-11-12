# geomTS

The goal of geomTS is to analysis geometrical time series data. For time series of full rank covariance (or correlation) matrices offering a systematic Riemannian framework for adapting statistical analysis to non-linear geometrical settings, this package defines different Riemannian metrics, including Euclidean metric in symmetric matrix space, affine invariant metric and Log-Cholesky metric in symmetric positive definite (SPD) matrix space, and quotient metric in quotient geometry consisting of full-rank correlation matrices. Under different manifold frameworks, we develop manifold-adapted time series models for matrix-valued data and explore how modelling results are influenced by manifold choices.

## Installation

You can install the development version of geomTS from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TaoDing2/geomTS")
library(geomTS)
```
