#' @title Circle fit
#'
#' @description Compute circle fit. Output projected positions.
#'
#' @param XY matrix object of two columns of measurements.
#'

circle.fit <- function(XY, IterMAX=50) {
  library(matrixStats)
  library(coincfit)

  parini <- colMeans(as.matrix(XY))

  cfit <- LMcircleFit(as.matrix(XY), ParIni=parini, IterMAX=IterMAX)
  cpred <- calculateCircle(cfit[1],cfit[2],cfit[3], 360)

  cpos <- do.call(rbind, lapply(1:nrow(XY), function(j) {
    obs <- unlist(XY[j,])
    diff <- sweep(cpred, 2, obs, FUN="-")
    rr <- rowSds(diff^2)
    pos <- which.min(rr)
    pos
  }))
  return(cpred[cpos,])
}



#' @title Circle fit sample deleted residuals
#'
#' @description Compute residuals of the circle fit after deleting each sample point
#'
#' @param XY matrix object of two columns of measurements.

circle.fit.resid.delete <- function(XY, IterMAX=50) {
  library(matrixStats)
  library(coincfit)

  resid.delete <- do.call(rbind, lapply(1:nrow(XY), function(i) {
    XY.del <- XY[-i,]
    parini <- colMeans(XY.del)

    cfit <- LMcircleFit(as.matrix(XY.del), ParIni=parini, IterMAX=IterMAX)
    cpred <- calculateCircle(cfit[1],cfit[2],cfit[3], 360)

    cpos <- do.call(rbind, lapply(1:nrow(XY.del), function(j) {
      obs <- unlist(XY.del[j,])
      diff <- sweep(cpred, 2, obs, FUN="-")
      rr <- rowSds(diff^2)
      pos <- which.min(rr)
      pos
    }))

    cresid <- sum(XY.del - cpred[cpos,])^2
    return(cresid)
  }) )

  return(resid.delete)
}



