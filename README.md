# geomTS

The goal of geomTS is to analysis geometrical time series data. For time series of full rank covariance (or correlation) matrices offering a systematic Riemannian framework for adapting statistical analysis to non-linear geometrical settings, this package defines different Riemannian metrics, including Euclidean metric in symmetric matrix space, affine invariant metric and Log-Cholesky metric in symmetric positive definite (SPD) matrix space, and quotient metric in quotient geometry consisting of full-rank correlation matrices. Under different manifold frameworks, we develop manifold-adapted time series models for matrix-valued data and explore how modelling results are influenced by manifold choices.

# Paper
Readers can get access to the arXiv paper (https://doi.org/10.48550/arXiv.2402.06410). There is the abstract of this paper:


We propose a model for time series taking values on a Riemannian manifold and fit it to time series of covariance matrices derived from EEG data for patients suffering from epilepsy. The aim of the study is two-fold: to develop a model with interpretable parameters for different possible modes of EEG dynamics, and to explore the extent to which modelling results are affected by the choice of manifold and its associated geometry. The model specifies a distribution for the tangent direction vector at any time point, combining an autoregressive term, a mean reverting term and a form of Gaussian noise. Parameter inference is carried out by maximum likelihood estimation, and we compare modelling results obtained using the standard Euclidean geometry on covariance matrices and the affine invariant geometry. Results distinguish between epileptic seizures and interictal periods between seizures in patients: between seizures the dynamics have a strong mean reverting component and the autoregressive component is missing, while for the majority of seizures there is a significant autoregressive component and the mean reverting effect is weak. The fitted models are also used to compare seizures within and between patients. The affine invariant geometry is advantageous and it provides a better fit to the data.


## Installation

You can install the development version of geomTS from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TaoDing2/geomTS")
library(geomTS)
```
