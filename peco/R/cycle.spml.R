#' @title Estimate gene weights for cell time
#'
#' @param Y gene by sample expression matrix
#' @param theta sample cell time vector
#'
#' @export
cycle.spml.trainmodel <- function(Y, theta) {

  library(Rfast)
  library(assertthat)
  fit <- spml.reg(theta, t(Y), seb=TRUE)
  return(fit)
}


#' @title Estimate gene weights for cell time
#'
#' @param Y_test gene by testing samples
#' @param theta_test gene by training samples
#' @param theta_train cell times for training samples
#' @param theta_test cell times for test samples
#'
#' @export
cycle.spml.testmodel <- function(Y_test, Y_train, theta_test, theta_train) {

  library(Rfast)
  library(assertthat)
  assert_that(is.matrix(Y_test))
  assert_that(dim(Y_test)[2]==length(theta_test),
              msg = "dimension of testing expression matrix doesn't match length of cell time vector")
  assert_that(is.matrix(Y_train))
  assert_that(dim(Y_train)[2]==length(theta_train),
              msg = "dimension of training expression matrix doesn't match length of cell time vector")

  fit_train <- cycle.spml.trainmodel(Y_train, theta_train)

  pred_cart <- cbind(1,t(Y_test))%*%fit_train$be
  pred_polar <- atan( pred_cart[, 2] / pred_cart[, 1] ) + pi * I(pred_cart[, 1] < 0)

  rho_test <- rFLIndTestRand(pred_polar, theta_test, 9999)
  boot_ci <- rhoFLCIBoot(pred_polar, theta_test, 95, 9999)

  return(list(betahat=fit_train$be,
              theta_pred=pred_polar,
              theta_test=theta_test,
              rho=rho_test[1],
              boot_95ci_low=boot_ci[1],
              boot_95ci_high=boot_ci[2],
              pval=rho_test[2]))
}
