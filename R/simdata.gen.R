#' Generate multi-source data from partially linear models.
#'
#' @description Generate simulation datasets containing training data and testing data from partially linear models under various settings.
#'
#' @param px the dimension of the shared parametric component for all models. Should be an integer smaller than sample size.
#' @param num.source the number of datasets. Should be the value 4 or 7.
#' @param size the sample size of different datasets. Should be a vector of \code{num.source}.
#' @param coeff0 a px * num.source matrix of the shared coefficient vector for all models.
#' @param coeff.mis the shared coefficient vector for the misspecified model. If tar.spec = 'cor', it should be a parameter vector of length px + 1 for the second misspecified source model. If tar.spec = 'mis', it should be a (px+1) * 2 matrix, in which the first column is the parameter vector for the misspecified target model and the second column is for the second misspecified source model. The last component of predictors for the misspecified model will be omitted in the estimation.
#' @param err.sigma the standard deviations of the normal random errors in regression models.
#' @param rho the correlation coefficient in the multivariate normal distribution of the parametric variables.
#' @param size.test the sample size of the testing target data.
#' @param sim.set the type of the nonparametric settings. Can be "heter" or "homo", which represents the heterogeneous and homogeneous dimension settings, respectively.
#' @param tar.spec the type of the target model specification. Can be "cor" or "mis", which represents the corrected and misspecified target model, respectively.
#' @param if.heter the logical variable, whether to allow a heteroscedastic setup. Default is False.
#'
#' @return a list of the training data and testing data, including the response, parametric predictors, nonparametric predictors, nonparametric values, coefficient vector.
#' @references Hu, X., & Zhang, X. (2023). Optimal Parameter-Transfer Learning by Semiparametric Model Averaging. Journal of Machine Learning Research, 24(358), 1-53.
#' @export
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom MASS mvrnorm
#'
#' @examples
#' coeff0 <- cbind(
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3)),
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.02),
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.3),
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3))
#' )
#' # correct target model setting
#' whole.data <- simdata.gen(
#'   px = 6, num.source = 4, size = c(150, 200, 200, 150), coeff0 = coeff0,
#'   coeff.mis = as.matrix(c(coeff0[, 2], 1.8)), err.sigma = 0.5, rho = 0.5, size.test = 500,
#'   sim.set = "homo", tar.spec = "cor", if.heter = FALSE
#' )
#'
#' # misspecified target model setting
#' coeff.mis <- matrix(c(c(coeff0[, 1], 0.1), c(coeff0[, 2], 1.8)), ncol = 2)
#' whole.data <- simdata.gen(
#'   px = 6, num.source = 4, size = c(150, 200, 200, 150), coeff0 = coeff0,
#'   coeff.mis = coeff.mis, err.sigma = 0.5, rho = 0.5, size.test = 500,
#'   sim.set = "homo", tar.spec = "mis", if.heter = FALSE
#' )
simdata.gen <- function(px, num.source = 4, size, coeff0, coeff.mis, err.sigma, rho, size.test, sim.set = c("heter", "homo"), tar.spec = c("cor", "mis"), if.heter = FALSE) {
  sim.set <- match.arg(arg = sim.set, choices = c("heter", "homo"))
  tar.spec <- match.arg(arg = tar.spec, choices = c("cor", "mis"))
  respon <- vector(mode = "list", length = num.source)
  beta.true <- vector(mode = "list", length = num.source)
  datax <- vector(mode = "list", length = num.source)
  dataz <- vector(mode = "list", length = num.source)
  gz <- vector(mode = "list", length = num.source)

  ## heterogeneous settings
  if (sim.set == "heter") {
    if (length(size) == 4) {
      qz <- c(3, 2, 2, 1)
      smooth.fun.true1 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true2 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 1]) + u[, 2])
      }
      smooth.fun.true3 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 2]) + u[, 1])
      }
      smooth.fun.true4 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 1]) + u[, 1])
      }
    } else if (length(size) == 7) {
      qz <- c(3, 2, 2, 1, 3, 2, 2)
      smooth.fun.true1 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true2 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 1]) + u[, 2])
      }
      smooth.fun.true3 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 2]) + u[, 1])
      }
      smooth.fun.true4 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 1]) + u[, 1])
      }
      smooth.fun.true5 <- function(u) {
        return((u[, 1] + 0.5)^3 + cos(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true6 <- function(u) {
        return((u[, 1] + 0.5)^3 + cos(pi * u[, 1]) + u[, 2])
      }
      smooth.fun.true7 <- function(u) {
        return((u[, 1] + 0.5)^3 + cos(pi * u[, 2]) + u[, 1])
      }
    }
  }

  ## homogeneous settings
  else if (sim.set == "homo") {
    qz <- rep(3, num.source)
    if (length(size) == 4) {
      smooth.fun.true1 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true2 <- function(u) {
        return(2 * (u[, 1] + 0.5)^3 + cos(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true3 <- function(u) {
        return(2.5 * (u[, 1] + 0.3)^3 + sin(pi * u[, 2]) + 1.5 * u[, 3])
      }
      smooth.fun.true4 <- function(u) {
        return((1.8 * u[, 1] + 0.3)^3 + cos(pi * u[, 2]) + u[, 3])
      }
    } else if (length(size) == 7) {
      smooth.fun.true1 <- function(u) {
        return(2 * (u[, 1] - 0.5)^3 + sin(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true2 <- function(u) {
        return(2 * (u[, 1] + 0.5)^3 + cos(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true3 <- function(u) {
        return(2.5 * (u[, 1] + 0.3)^3 + sin(pi * u[, 2]) + 1.5 * u[, 3])
      }
      smooth.fun.true4 <- function(u) {
        return((1.8 * u[, 1] + 0.3)^3 + cos(pi * u[, 2]) + u[, 3])
      }
      smooth.fun.true5 <- function(u) {
        return(1.5 * (u[, 1] + 0.5)^3 + cos(2 * pi * u[, 2]) + u[, 3]^2)
      }
      smooth.fun.true6 <- function(u) {
        return((u[, 1] + 0.6)^2 + cos(pi * u[, 2]) + 1.3 * u[, 3]^2)
      }
      smooth.fun.true7 <- function(u) {
        return(1.3 * (u[, 1] + 0.5)^2 + cos(2 * pi * u[, 2]) + 1.6 * u[, 3])
      }
    }
  }

  if (tar.spec == "cor") {
    pnew <- length(coeff.mis)
    # train data
    for (i in 1:num.source) {
      for (j in 1:qz[i]) {
        dataz[[i]] <- cbind(dataz[[i]], as.matrix(runif(size[i], 0, 1)))
      }
      fun.name <- get(paste("smooth.fun.true", i, sep = ""))
      gz[[i]] <- as.matrix(fun.name(dataz[[i]]))
    }
    for (k in 1:num.source) {
      if (k == 2) {
        beta.true[[k]] <- coeff.mis
        xmat <- rho^abs(outer(1:pnew, 1:pnew, "-"))
        datax[[k]] <- mvrnorm(size[k], rep(0, pnew), xmat)
        if (if.heter) {
          err <- c()
          for (ii in 1:size[k]) {
            err[ii] <- rnorm(1, 0, err.sigma * (datax[[k]][ii, 1])^2)
          }
        } else {
          err <- rnorm(size[k], 0, err.sigma)
        }
        respon[[k]] <- datax[[k]] %*% beta.true[[k]] + gz[[k]] + err
      } else {
        beta.true[[k]] <- coeff0[, k]
        xmat <- rho^abs(outer(1:px, 1:px, "-"))
        datax[[k]] <- mvrnorm(size[k], rep(0, px), xmat)
        if (if.heter) {
          err <- c()
          for (ii in 1:size[k]) {
            err[ii] <- rnorm(1, 0, err.sigma * (datax[[k]][ii, 1])^2)
          }
        } else {
          err <- rnorm(size[k], 0, err.sigma)
        }
        respon[[k]] <- datax[[k]] %*% beta.true[[k]] + gz[[k]] + err
      }
    }
    data.train <- list(data.y = respon, beta.true = beta.true, data.x = datax, data.z = dataz, gz = gz)

    # test data
    beta.true.te <- coeff0[, 1]
    xmat.te <- rho^abs(outer(1:px, 1:px, "-"))
    datax.te <- mvrnorm(size.test, rep(0, px), xmat.te)
    if (if.heter) {
      err.te <- c()
      for (ii in 1:size.test) {
        err.te[ii] <- rnorm(1, 0, err.sigma * (datax.te[ii, 1])^2)
      }
    } else {
      err.te <- rnorm(size.test, 0, err.sigma)
    }
    dataz.te <- NULL
    for (j in 1:qz[1]) {
      dataz.te <- cbind(dataz.te, as.matrix(runif(size.test, 0, 1)))
    }
    gz.te <- as.matrix(smooth.fun.true1(dataz.te))
    respon.te <- datax.te %*% beta.true.te + gz.te + err.te
    data.test <- list(data.y = respon.te, data.x = datax.te, data.z = dataz.te, gz.te = gz.te, beta.true = beta.true.te)
  } else if (tar.spec == "mis") {
    pnew <- nrow(coeff.mis)
    para.mis.tar <- coeff.mis[, 1]
    para.mis.src <- coeff.mis[, 2]

    # train data
    for (i in 1:num.source) {
      for (j in 1:qz[i]) {
        dataz[[i]] <- cbind(dataz[[i]], as.matrix(runif(size[i], 0, 1)))
      }
      fun.name <- get(paste("smooth.fun.true", i, sep = ""))
      gz[[i]] <- as.matrix(fun.name(dataz[[i]]))
    }
    beta.true[[1]] <- para.mis.tar
    xmat <- rho^abs(outer(1:pnew, 1:pnew, "-"))
    datax[[1]] <- mvrnorm(size[1], rep(0, pnew), xmat)
    if (if.heter) {
      err <- c()
      for (ii in 1:size[1]) {
        err[ii] <- rnorm(1, 0, err.sigma * (datax[[1]][ii, 1])^2)
      }
    } else {
      err <- rnorm(size[1], 0, err.sigma)
    }
    respon[[1]] <- datax[[1]] %*% beta.true[[1]] + gz[[1]] + err

    for (k in 2:num.source) {
      if (k == 2) {
        beta.true[[k]] <- para.mis.src
        xmat <- rho^abs(outer(1:pnew, 1:pnew, "-"))
        datax[[k]] <- mvrnorm(size[k], rep(0, pnew), xmat)
        if (if.heter) {
          err <- c()
          for (ii in 1:size[k]) {
            err[ii] <- rnorm(1, 0, err.sigma * (datax[[k]][ii, 1])^2)
          }
        } else {
          err <- rnorm(size[k], 0, err.sigma)
        }
        respon[[k]] <- datax[[k]] %*% beta.true[[k]] + gz[[k]] + err
      } else {
        beta.true[[k]] <- coeff0[, k]
        xmat <- rho^abs(outer(1:px, 1:px, "-"))
        datax[[k]] <- mvrnorm(size[k], rep(0, px), xmat)
        if (if.heter) {
          err <- c()
          for (ii in 1:size[k]) {
            err[ii] <- rnorm(1, 0, err.sigma * (datax[[k]][ii, 1])^2)
          }
        } else {
          err <- rnorm(size[k], 0, err.sigma)
        }
        respon[[k]] <- datax[[k]] %*% beta.true[[k]] + gz[[k]] + err
      }
    }
    data.train <- list(data.y = respon, beta.true = beta.true, data.x = datax, data.z = dataz, gz = gz)

    # test data
    beta.true.te <- para.mis.tar
    xmat.te <- rho^abs(outer(1:pnew, 1:pnew, "-"))
    datax.te <- mvrnorm(size.test, rep(0, pnew), xmat.te)
    if (if.heter) {
      err.te <- c()
      for (ii in 1:size.test) {
        err.te[ii] <- rnorm(1, 0, err.sigma * (datax.te[ii, 1])^2)
      }
    } else {
      err.te <- rnorm(size.test, 0, err.sigma)
    }
    dataz.te <- NULL
    for (j in 1:qz[1]) {
      dataz.te <- cbind(dataz.te, as.matrix(runif(size.test, 0, 1)))
    }
    gz.te <- as.matrix(smooth.fun.true1(dataz.te))
    respon.te <- datax.te %*% beta.true.te + gz.te + err.te
    data.test <- list(data.y = respon.te, data.x = datax.te, data.z = dataz.te, gz.te = gz.te, beta.true = beta.true.te)
  }

  return(list(data.train = data.train, data.test = data.test))
}
