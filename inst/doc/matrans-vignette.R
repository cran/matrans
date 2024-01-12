## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = FALSE-------------------------------------------------------------
library(formatR)

## ----eval=FALSE---------------------------------------------------------------
#  library("devtools")
#  devtools::install_github("XnhuUcas/matrans")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("matrans")

## -----------------------------------------------------------------------------
library(matrans)

## ----tidy=TRUE, tidy.opts=list(width.cutoff=70)-------------------------------
set.seed(1)
## sample size
size <- c(150, 200, 200, 150)
## shared coefficient vectors for different models
coeff0 <- cbind(
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3)),
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.02),
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.3),
  as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3))
)
## dimension of parametric component for all models
px <- 6
## standard deviation for random errors
err.sigma <- 0.5
## the correlation coefficient for covariates
rho <- 0.5
## sample size for testing data
size.test <- 500

whole.data <- simdata.gen(
  px = px, num.source = 4, size = size, coeff0 = coeff0, coeff.mis = as.matrix(c(coeff0[, 2], 1.8)),
  err.sigma = err.sigma, rho = rho, size.test = size.test, sim.set = "homo", tar.spec = "cor",
  if.heter = FALSE
)
## multi-source training datasets
data.train <- whole.data$data.train
## testing target dataset
data.test <- whole.data$data.test

## ----tidy=TRUE, tidy.opts=list(width.cutoff=70), collapse=FALSE, comment='', warning=FALSE----
## hyperparameters for B-splines
bs.para <- list(bs.df = rep(3, 3), bs.degree = rep(3, 3))
## the second model is misspecified
data.train$data.x[[2]] <- data.train$data.x[[2]][, -7]
## fitting the Trans-SMAP
fit.transsmap <- trans.smap(train.data = data.train, nfold = 5, bs.para = bs.para)
## weight estimates
fit.transsmap$weight.est
## computational time of algorithm (sec)
fit.transsmap$time.transsmap

## ----tidy=TRUE, tidy.opts=list(width.cutoff=70), collapse=FALSE, comment='', warning=FALSE----
## prediction using testing data
pred.res <- pred.transsmap(object = fit.transsmap, newdata = data.test, bs.para = bs.para)
## predicted values for the new observations of predictors
pred.val <- pred.res$predict.val
## mean squared prediction risk for Trans-SMAP
sum((pred.val - data.test$data.x %*% data.test$beta.true - data.test$gz.te)^2) / 500

