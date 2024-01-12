#' Parameter-transfer learning for partially linear models based on semiparametric model averaging.
#'
#' @description Obtain optimal weights and estimated coefficients based on Trans-SMAP.
#'
#' @param train.data a list containing the observations of predictors and response for fitting models. Should be a list with elements "data.y", "data.x" and "data.z", where "data.y" indicates a response list for all data sources, "data.x" indicates a parametric predictor list for all data sources, and "data.z" indicates a nonparametric predictor list for all data sources. Each element in "data.x" and "data.z" is a matrix with each row as an observation and each column as a variable. By default, the first element in "data.y", "data.x" and "data.z" is target data, and others are source data.
#' @param nfold the number of folds for the cross-validation weight criterion. Default is NULL (leave-one-out).
#' @param bs.para a list containing the parameters for B-spline construction in function \code{bs}. Should be a list with elements "bs.df" and "bs.degree", each component of which is a vector with the same length as the number of nonparametric variables. For example, bs.para = list(bs.df=c(3,3,3), bs.degree=c(3,3,3)).
#' \itemize{
#' \item "bs.df": degrees of freedom for each nonparametric component; The details can be referred to the arguments in function \code{bs}.
#' \item "bs.degree": degree of the piecewise polynomial for each nonparametric component; The default is 3 for cubic splines.
#' }
#' @param lm.set the vector of indices for the linear regression models, which means the corresponding models are constructed by ordinary linear models instead of partially linear models. Default is NULL.
#'
#' @return a result list containing the estimated weight vector, the execution time of solving the optimal weights and the summarized results of fitting models.
#' @references Hu, X., & Zhang, X. (2023). Optimal Parameter-Transfer Learning by Semiparametric Model Averaging. Journal of Machine Learning Research, 24(358), 1-53.
#' @export
#' @importFrom splines bs
#' @importFrom quadprog solve.QP
#' @importFrom stats lm
#' @importFrom stats as.formula
#' @importFrom caret createFolds
#'
#' @examples
#' ## correct target model setting
#'
#' # generate simulation dataset
#' coeff0 <- cbind(
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3)),
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.02),
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3) + 0.3),
#'   as.matrix(c(1.4, -1.2, 1, -0.8, 0.65, 0.3))
#' )
#' whole.data <- simdata.gen(
#'   px = 6, num.source = 4, size = c(150, 200, 200, 150), coeff0 = coeff0,
#'   coeff.mis = as.matrix(c(coeff0[, 2], 1.8)), err.sigma = 0.5, rho = 0.5, size.test = 500,
#'   sim.set = "homo", tar.spec = "cor", if.heter = FALSE
#' )
#' data.train <- whole.data$data.train
#' data.test <- whole.data$data.test
#'
#' # running Trans-SMAP and obtain the optimal weight vector
#' data.train$data.x[[2]] <- data.train$data.x[[2]][, -7]
#' fit.transsmap <- trans.smap(
#'   train.data = data.train, nfold = 5,
#'   bs.para = list(bs.df = rep(3, 3), bs.degree = rep(3, 3))
#' )
#' ma.weights <- fit.transsmap$weight.est
#'
#' \donttest{
#' ## misspecified target model setting
#'
#' # generate simulation dataset
#' coeff.mis <- matrix(c(c(coeff0[, 1], 0.1), c(coeff0[, 2], 1.8)), ncol = 2)
#' whole.data <- simdata.gen(
#'   px = 6, num.source = 4, size = c(150, 200, 200, 150), coeff0 = coeff0,
#'   coeff.mis = coeff.mis, err.sigma = 0.5, rho = 0.5, size.test = 500,
#'   sim.set = "homo", tar.spec = "mis", if.heter = FALSE
#' )
#' data.train <- whole.data$data.train
#' data.test <- whole.data$data.test
#'
#' # running Trans-SMAP and obtain the optimal weight vector
#' data.train$data.x[[1]] <- data.train$data.x[[1]][, -7]
#' data.train$data.x[[2]] <- data.train$data.x[[2]][, -7]
#' fit.transsmap <- trans.smap(
#'   train.data = data.train, nfold = 5,
#'   bs.para = list(bs.df = rep(3, 3), bs.degree = rep(3, 3))
#' )
#' ma.weights <- fit.transsmap$weight.est
#' }
trans.smap <- function(train.data, nfold = NULL, bs.para, lm.set = NULL) {
  s.cnt <- length(train.data$data.y)
  p <- ncol(train.data$data.x[[1]])
  size.train <- sapply(train.data$data.x, function(x) nrow(x))
  bs.df <- bs.para$bs.df
  bs.degree <- bs.para$bs.degree
  xstring <- paste(paste("x", 1:p, sep = ""), collapse = "+")
  reg.res <- vector(mode = "list", length = s.cnt)
  for (k in 1:s.cnt) {
    data.merge <- cbind(train.data$data.y[[k]], train.data$data.x[[k]], train.data$data.z[[k]])
    data.train.frame <- as.data.frame(data.merge)
    colnames(data.train.frame) <- c("respon", paste("x", 1:p, sep = ""), paste("z", 1:ncol(train.data$data.z[[k]]), sep = ""))
    if (k %in% lm.set) {
      reg.res[[k]] <- lm(respon ~ ., data = data.train.frame)
    } else {
      zstring <- paste(paste(paste("bs(", paste(paste("z", 1:ncol(train.data$data.z[[k]]), sep = ""),
        paste("df = bs.df[", 1:ncol(train.data$data.z[[k]]), "]", sep = ""),
        paste("degree = bs.degree[", 1:ncol(train.data$data.z[[k]]), "]", sep = ""),
        sep = ","
      ), sep = ""), ")", sep = ""), collapse = "+")
      reg.res[[k]] <- lm(as.formula(paste("respon~", paste(xstring, zstring, sep = "+"), sep = "")), data = data.train.frame)
    }
  }

  if (is.null(lm.set)) {
    ## jackknife
    if (is.null(nfold)) {
      timestart <- Sys.time()
      proj.mat <- matrix(0, s.cnt, s.cnt)
      for (i in 1:size.train[1]) {
        train.data.cv <- train.data
        est.beta <- matrix(NA, nrow = s.cnt, ncol = p)
        lm.tr <- vector(mode = "list", length = s.cnt)
        for (j in 1:s.cnt) {
          if (j == 1) {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
            data.merge.cv <- data.merge.cv[-i, ]
          } else {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
          }
          qz <- ncol(train.data.cv$data.z[[j]])
          datalm <- as.data.frame(data.merge.cv)
          colnames(datalm) <- c("respon", paste("x", 1:p, sep = ""), paste("z", 1:qz, sep = ""))
          xstring <- paste(paste("x", 1:p, sep = ""), collapse = "+")
          zstring <- paste(paste(paste("bs(", paste(paste("z", 1:qz, sep = ""),
            paste("df = bs.df[", 1:qz, "]", sep = ""),
            paste("degree = bs.degree[", 1:qz, "]", sep = ""),
            sep = ","
          ), sep = ""), ")", sep = ""), collapse = "+")
          lm.tr[[j]] <- lm(as.formula(paste("respon~", paste(xstring, zstring, sep = "+"), sep = "")), data = datalm)$coefficients
          est.beta[j, ] <- lm.tr[[j]][2:(p + 1)]
        }
        beta.est.train.mat <- NULL
        for (j in 1:s.cnt) {
          beta.est.train.mat <- cbind(beta.est.train.mat, as.matrix(c(lm.tr[[1]][1], est.beta[j, ], lm.tr[[1]][(p + 2):length(lm.tr[[1]])])))
        }
        q <- ncol(train.data$data.z[[1]])
        bsz.tar <- NULL
        for (j in 1:q) {
          bsz.tar <- cbind(bsz.tar, bs(train.data$data.z[[1]][, j]))
        }
        data.merge.new <- cbind(train.data$data.x[[1]], bsz.tar)
        pred.y <- t(c(1, data.merge.new[i, ]) %*% beta.est.train.mat)
        proj.mat <- proj.mat + (pred.y - train.data$data.y[[1]][i] * matrix(1, s.cnt, 1)) %*% t(pred.y - train.data$data.y[[1]][i] * matrix(1, s.cnt, 1))
      }

      Dmat <- proj.mat / size.train[1]
      d <- rep(0, s.cnt)
      Amat <- t(rbind(matrix(1, nrow = 1, ncol = s.cnt), diag(s.cnt), -diag(s.cnt)))
      bvec <- rbind(1, matrix(0, nrow = s.cnt, ncol = 1), matrix(-1, nrow = s.cnt, ncol = 1))
      solve.qr <- try(
        {
          solve.QP(Dmat, d, Amat, bvec, meq = 1)
        },
        silent = TRUE
      )
      time.transsmap <- as.double(difftime(Sys.time(), timestart, units = "secs"))
      if ("try-error" %in% class(solve.qr)) {
        return(list(weight.est = NA, time.transsmap = time.transsmap, reg.res = reg.res))
      } else {
        return(list(weight.est = solve.qr$solution, time.transsmap = time.transsmap, reg.res = reg.res))
      }
    }

    ## k-fold
    else {
      timestart <- Sys.time()
      proj.mat <- matrix(0, s.cnt, s.cnt)
      split <- createFolds(train.data$data.y[[1]], k = nfold)

      for (i in 1:nfold) {
        train.data.cv <- train.data
        est.beta <- matrix(NA, nrow = s.cnt, ncol = p)
        lm.tr <- vector(mode = "list", length = s.cnt)
        for (j in 1:s.cnt) {
          if (j == 1) {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
            data.merge.cv <- data.merge.cv[-split[[i]], ]
          } else {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
          }
          qz <- ncol(train.data.cv$data.z[[j]])
          datalm <- as.data.frame(data.merge.cv)
          colnames(datalm) <- c("respon", paste("x", 1:p, sep = ""), paste("z", 1:qz, sep = ""))
          xstring <- paste(paste("x", 1:p, sep = ""), collapse = "+")
          zstring <- paste(paste(paste("bs(", paste(paste("z", 1:qz, sep = ""),
            paste("df = bs.df[", 1:qz, "]", sep = ""),
            paste("degree = bs.degree[", 1:qz, "]", sep = ""),
            sep = ","
          ), sep = ""), ")", sep = ""), collapse = "+")
          lm.tr[[j]] <- lm(as.formula(paste("respon~", paste(xstring, zstring, sep = "+"), sep = "")), data = datalm)$coefficients
          est.beta[j, ] <- lm.tr[[j]][2:(p + 1)]
        }
        beta.est.train.mat <- NULL
        for (j in 1:s.cnt) {
          beta.est.train.mat <- cbind(beta.est.train.mat, as.matrix(c(lm.tr[[1]][1], est.beta[j, ], lm.tr[[1]][(p + 2):length(lm.tr[[1]])])))
        }
        q <- ncol(train.data$data.z[[1]])
        bsz.tar <- NULL
        for (j in 1:q) {
          bsz.tar <- cbind(bsz.tar, bs(train.data$data.z[[1]][, j]))
        }
        data.merge.new <- cbind(as.matrix(rep(1, size.train[1])), train.data$data.x[[1]], bsz.tar)
        pred.y <- data.merge.new[split[[i]], ] %*% beta.est.train.mat
        err.fold <- t(pred.y - train.data$data.y[[1]][split[[i]], ] %*% matrix(1, 1, s.cnt)) %*% (pred.y - train.data$data.y[[1]][split[[i]], ] %*% matrix(1, 1, s.cnt))
        proj.mat <- proj.mat + err.fold
      }

      Dmat <- proj.mat / size.train[1]
      d <- rep(0, s.cnt)
      Amat <- t(rbind(matrix(1, nrow = 1, ncol = s.cnt), diag(s.cnt), -diag(s.cnt)))
      bvec <- rbind(1, matrix(0, nrow = s.cnt, ncol = 1), matrix(-1, nrow = s.cnt, ncol = 1))
      solve.qr <- try(
        {
          solve.QP(Dmat, d, Amat, bvec, meq = 1)
        },
        silent = TRUE
      )
      time.transsmap <- as.double(difftime(Sys.time(), timestart, units = "secs"))
      if ("try-error" %in% class(solve.qr)) {
        return(list(weight.est = NA, time.transsmap = time.transsmap, reg.res = reg.res))
      } else {
        return(list(weight.est = solve.qr$solution, time.transsmap = time.transsmap, reg.res = reg.res))
      }
    }
  } else {
    ## jackknife
    if (is.null(nfold)) {
      timestart <- Sys.time()
      proj.mat <- matrix(0, s.cnt, s.cnt)
      for (i in 1:size.train[1]) {
        train.data.cv <- train.data
        est.beta <- matrix(NA, nrow = s.cnt, ncol = p)
        lm.tr <- vector(mode = "list", length = s.cnt)
        for (j in 1:s.cnt) {
          if (j == 1) {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
            data.merge.cv <- data.merge.cv[-i, ]
          } else {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
          }
          qz <- ncol(train.data.cv$data.z[[j]])
          datalm <- as.data.frame(data.merge.cv)
          colnames(datalm) <- c("respon", paste("x", 1:p, sep = ""), paste("z", 1:qz, sep = ""))
          if (j %in% lm.set) {
            lm.tr[[j]] <- lm(respon ~ ., data = datalm)$coefficients
          } else {
            xstring <- paste(paste("x", 1:p, sep = ""), collapse = "+")
            zstring <- paste(paste(paste("bs(", paste(paste("z", 1:qz, sep = ""),
              paste("df = bs.df[", 1:qz, "]", sep = ""),
              paste("degree = bs.degree[", 1:qz, "]", sep = ""),
              sep = ","
            ), sep = ""), ")", sep = ""), collapse = "+")
            lm.tr[[j]] <- lm(as.formula(paste("respon~", paste(xstring, zstring, sep = "+"), sep = "")), data = datalm)$coefficients
          }
          est.beta[j, ] <- lm.tr[[j]][2:(p + 1)]
        }
        beta.est.train.mat <- NULL
        for (j in 1:s.cnt) {
          beta.est.train.mat <- cbind(beta.est.train.mat, as.matrix(c(lm.tr[[1]][1], est.beta[j, ], lm.tr[[1]][(p + 2):length(lm.tr[[1]])])))
        }
        q <- ncol(train.data$data.z[[1]])
        if (1 %in% lm.set) {
          data.merge.new <- cbind(train.data$data.x[[1]], train.data$data.z[[1]])
        } else {
          bsz.tar <- NULL
          for (j in 1:q) {
            bsz.tar <- cbind(bsz.tar, bs(train.data$data.z[[1]][, j]))
          }
          data.merge.new <- cbind(train.data$data.x[[1]], bsz.tar)
        }
        pred.y <- t(c(1, data.merge.new[i, ]) %*% beta.est.train.mat)
        proj.mat <- proj.mat + (pred.y - train.data$data.y[[1]][i] * matrix(1, s.cnt, 1)) %*% t(pred.y - train.data$data.y[[1]][i] * matrix(1, s.cnt, 1))
      }

      Dmat <- proj.mat / size.train[1]
      d <- rep(0, s.cnt)
      Amat <- t(rbind(matrix(1, nrow = 1, ncol = s.cnt), diag(s.cnt), -diag(s.cnt)))
      bvec <- rbind(1, matrix(0, nrow = s.cnt, ncol = 1), matrix(-1, nrow = s.cnt, ncol = 1))
      solve.qr <- try(
        {
          solve.QP(Dmat, d, Amat, bvec, meq = 1)
        },
        silent = TRUE
      )
      time.transsmap <- as.double(difftime(Sys.time(), timestart, units = "secs"))
      if ("try-error" %in% class(solve.qr)) {
        return(list(weight.est = NA, time.transsmap = time.transsmap, reg.res = reg.res))
      } else {
        return(list(weight.est = solve.qr$solution, time.transsmap = time.transsmap, reg.res = reg.res))
      }
    }

    ## k-fold
    else {
      timestart <- Sys.time()
      proj.mat <- matrix(0, s.cnt, s.cnt)
      split <- createFolds(train.data$data.y[[1]], k = nfold)

      for (i in 1:nfold) {
        train.data.cv <- train.data
        est.beta <- matrix(NA, nrow = s.cnt, ncol = p)
        lm.tr <- vector(mode = "list", length = s.cnt)
        for (j in 1:s.cnt) {
          if (j == 1) {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
            data.merge.cv <- data.merge.cv[-split[[i]], ]
          } else {
            data.merge.cv <- cbind(train.data$data.y[[j]], train.data$data.x[[j]], train.data$data.z[[j]])
          }
          qz <- ncol(train.data.cv$data.z[[j]])
          datalm <- as.data.frame(data.merge.cv)
          colnames(datalm) <- c("respon", paste("x", 1:p, sep = ""), paste("z", 1:qz, sep = ""))
          if (j %in% lm.set) {
            lm.tr[[j]] <- lm(respon ~ ., data = datalm)$coefficients
          } else {
            xstring <- paste(paste("x", 1:p, sep = ""), collapse = "+")
            zstring <- paste(paste(paste("bs(", paste(paste("z", 1:qz, sep = ""),
              paste("df = bs.df[", 1:qz, "]", sep = ""),
              paste("degree = bs.degree[", 1:qz, "]", sep = ""),
              sep = ","
            ), sep = ""), ")", sep = ""), collapse = "+")
            lm.tr[[j]] <- lm(as.formula(paste("respon~", paste(xstring, zstring, sep = "+"), sep = "")), data = datalm)$coefficients
          }
          est.beta[j, ] <- lm.tr[[j]][2:(p + 1)]
        }
        beta.est.train.mat <- NULL
        for (j in 1:s.cnt) {
          beta.est.train.mat <- cbind(beta.est.train.mat, as.matrix(c(lm.tr[[1]][1], est.beta[j, ], lm.tr[[1]][(p + 2):length(lm.tr[[1]])])))
        }
        q <- ncol(train.data$data.z[[1]])
        if (1 %in% lm.set) {
          data.merge.new <- cbind(1, train.data$data.x[[1]], train.data$data.z[[1]])
        } else {
          bsz.tar <- NULL
          for (j in 1:q) {
            bsz.tar <- cbind(bsz.tar, bs(train.data$data.z[[1]][, j]))
          }
          data.merge.new <- cbind(as.matrix(rep(1, size.train[1])), train.data$data.x[[1]], bsz.tar)
        }
        pred.y <- data.merge.new[split[[i]], ] %*% beta.est.train.mat
        err.fold <- t(pred.y - train.data$data.y[[1]][split[[i]], ] %*% matrix(1, 1, s.cnt)) %*% (pred.y - train.data$data.y[[1]][split[[i]], ] %*% matrix(1, 1, s.cnt))
        proj.mat <- proj.mat + err.fold
      }

      Dmat <- proj.mat / size.train[1]
      d <- rep(0, s.cnt)
      Amat <- t(rbind(matrix(1, nrow = 1, ncol = s.cnt), diag(s.cnt), -diag(s.cnt)))
      bvec <- rbind(1, matrix(0, nrow = s.cnt, ncol = 1), matrix(-1, nrow = s.cnt, ncol = 1))
      solve.qr <- try(
        {
          solve.QP(Dmat, d, Amat, bvec, meq = 1)
        },
        silent = TRUE
      )
      time.transsmap <- as.double(difftime(Sys.time(), timestart, units = "secs"))
      if ("try-error" %in% class(solve.qr)) {
        return(list(weight.est = NA, time.transsmap = time.transsmap, reg.res = reg.res))
      } else {
        return(list(weight.est = solve.qr$solution, time.transsmap = time.transsmap, reg.res = reg.res))
      }
    }
  }
}
