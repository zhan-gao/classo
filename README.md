# classo

<!-- badges: start -->
<!-- badges: end -->

This is the package implements Classifier-Lasso developed in the following paper.

[Su, L., Shi, Z., & Phillips, P. C. (2016). Identifying latent structures in panel data. *Econometrica*, *84*(6), 2215-2264.](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA12560)

An accompanying note focused on the computation details: [Two Examples of Convex-Programming-Based High-Dimensional Econometric Estimators](https://arxiv.org/abs/1806.10423) and illustration example code hosted on https://github.com/zhan-gao/convex_prog_in_econometrics .



The package is still under active development...

## Installation

This package is dependent on [Rmosek](https://cran.r-project.org/web/packages/Rmosek/index.html). An installation gist can be found at https://gist.github.com/mikelove/67ea44d5be5a053e599257fe357483dc . Please make sure `Rmosek` is successfully installed and activated before install this package.

In later versions, we will relax this restriction to allow users use an open source solver via [CVXR](https://github.com/anqif/CVXR) and make `Rmosek` as an option.

You can install the development version of classo from [Github](https://CRAN.R-project.org) with:

``` r
library(devtools)
devtools::install_github("zhan-gao/classo", INSTALL_opts=c("--no-multiarch"))
```

