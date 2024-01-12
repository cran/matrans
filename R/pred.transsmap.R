#' Prediction for new data based on Trans-SMAP.
#'
#' @description Obtain predictions from a "trans.smap" object based on new samples.
#'
#' @param object the output of function \code{trans.smap}.
#' @param newdata a list containing the new observations of predictors for prediction, the components of which is named as "data.x" for parametric variables and "data.z" for nonparametric variables. Should be in accordance with the data for training \code{object}.
#' @param bs.para a list containing the parameters for B-spline construction in function \code{bs}. Should be a vector with names "bs.df" and "bs.degree", each component of which is a vector with the same length as the number of nonparametric variables. For example, bs.para = list(bs.df=c(3,3,3), bs.degree=c(3,3,3)).
#' \itemize{
#' \item "bs.df": degrees of freedom for each nonparametric component; The details can be referred to the arguments in function \code{bs}.
#' \item "bs.degree": degree of the piecewise polynomial for each nonparametric component; The default is 3 for cubic splines.
#' }
#' @param if.lm the logical variable, whether to set the target model as ordinary linear model. Default is False.
#'
#' @return a result list containing the predicted values on new data and the estimated coefficient vector.
#' @seealso \code{\link{trans.smap}}.
#' @references Hu, X., & Zhang, X. (2023). Optimal Parameter-Transfer Learning by Semiparametric Model Averaging. Journal of Machine Learning Research, 24(358), 1-53.
#' @export
#' @importFrom splines bs
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
#' # predict for new data
#' pred.res <- pred.transsmap(
#'   object = fit.transsmap, newdata = data.test,
#'   bs.para = list(bs.df = rep(3, 3), bs.degree = rep(3, 3))
#' )
#' pred.val <- pred.res$predict.val
#' predict.risk <- sum((pred.val - data.test$data.x %*% data.test$beta.true - data.test$gz.te)^2) / 500
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
#'
#' # predict for new data
#' data.test.mis <- data.test
#' data.test.mis$data.x <- data.test.mis$data.x[, -7]
#' pred.res <- pred.transsmap(
#'   object = fit.transsmap, newdata = data.test.mis,
#'   bs.para = list(bs.df = rep(3, 3), bs.degree = rep(3, 3))
#' )
#' pred.val <- pred.res$predict.val
#' predict.risk <- sum((pred.val - data.test$data.x %*% data.test$beta.true - data.test$gz.te)^2) / 500
#' }
pred.transsmap <- function(object, newdata, bs.para, if.lm = FALSE) {
  q <- ncol(newdata$data.z)
  p <- ncol(newdata$data.x)
  size.test <- nrow(newdata$data.x)
  bs.df <- bs.para$bs.df
  bs.degree <- bs.para$bs.degree
  reg.res <- object$reg.res
  ma.weights <- object$weight.est

  beta.est.train.mat <- matrix(NA, p, length(reg.res))
  for (k in 1:length(reg.res)) {
    beta.est.train.mat[, k] <- reg.res[[k]]$coefficients[2:(p + 1)]
  }
  beta.ma <- beta.est.train.mat %*% as.matrix(ma.weights)

  if (if.lm) {
    nonpara.est <- newdata$data.z %*% reg.res[[1]]$coefficients[(p + 2):length(reg.res[[1]]$coefficients)]
  } else {
    bsz.tar <- NULL
    for (j in 1:q) {
      bsz.tar <- cbind(bsz.tar, bs(newdata$data.z[, j], df = bs.df[j], degree = bs.degree[j]))
    }
    nonpara.est <- bsz.tar %*% reg.res[[1]]$coefficients[(p + 2):length(reg.res[[1]]$coefficients)]
  }

  if (all(is.na(ma.weights))) {
    return(list(predict.val = NA, beta.ma = NA))
  } else {
    predict.val <- (newdata$data.x %*% beta.est.train.mat + matrix(rep(reg.res[[1]]$coefficients[1] + nonpara.est, length(reg.res)), size.test, length(reg.res))) %*% as.matrix(ma.weights)
    return(list(
      predict.val = predict.val,
      beta.ma = beta.ma
    ))
  }
}
