
<!-- README.md is generated from README.Rmd. Please edit that file -->

# matrans

<!-- badges: start -->
<!-- badges: end -->

This package provides prediction tools under multi-source transfer
learning framework based on frequentist model averaging strategy. It is
primarily built on statistical model frameworks, including linear
regression models, partially linear models. Unlike existing approaches,
the proposed methods can adaptively integrate the auxiliary information
from different sources and possess asymptotic optimality for prediction
on the target model. It is worth noting that this package is the first
open-source software package for transfer learning based on the optimal
model averaging methods, providing convenient and efficient
computational tools for practitioners in multi-source data modeling. For
specific details, please refer to the following literature:

\[1\] Hu, X., & Zhang, X. (2023). [Optimal Parameter-Transfer Learning
by Semiparametric Model
Averaging](https://jmlr.org/papers/v24/23-0030.html). Journal of Machine
Learning Research, 24(358), 1-53.

Any questions or comments, please don’t hesitate to contact with me any
time.

## Installation

You can install the development version of the package like so:

``` r
library("devtools")
devtools::install_github("XnhuUcas/matrans")
```

## Maintainer

Xiaonan Hu (<xiaonanhu@cnu.edu.cn>)

## Usage

This is a simple example which shows users how to use the provided
functions to estimate weights and make predictions.

First, we generate simulation datasets (M=3) under the corrected target
model and homogeneous dimension settings.

``` r
library(matrans)

## generate simulation datasets (M=3)
size <- c(150, 200, 200, 150)
coeff0 <- cbind(
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3)),
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.02),
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.3),
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3))
)
px <- 6
err.sigma <- 0.5
rho <- 0.5
size.test <- 500

whole.data <- simdata.gen(
  px = px, num.source = 4, size = size, coeff0 = coeff0, coeff.mis = as.matrix(c(coeff0[, 2], 1.8)),
  err.sigma = err.sigma, rho = rho, size.test = size.test, sim.set = "homo", tar.spec = "cor",
  if.heter = FALSE
)
data.train <- whole.data$data.train
data.test <- whole.data$data.test
```

Then, we apply the functions to implement weights estimation and
out-of-sample predictions.

``` r
## optimal weights estimation
bs.para <- list(bs.df = rep(3, 3), bs.degree = rep(3, 3))
data.train$data.x[[2]] <- data.train$data.x[[2]][, -7]
fit.transsmap <- trans.smap(train.data = data.train, nfold = 5, bs.para = bs.para)
ma.weights <- fit.transsmap$weight.est
time.transsmap <- fit.transsmap$time.transsmap

## out-of-sample prediction results
pred.res <- pred.transsmap(object = fit.transsmap, newdata = data.test, bs.para = bs.para)
pred.val <- pred.res$predict.val
```
