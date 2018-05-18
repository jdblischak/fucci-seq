#' @title compute Fisher-Lee correlation coefficient
#'
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#'
#' @references Pewsey et al. Circular statistics in R
#'
#' @export
rFLCorrCoeff <- function(lcdat1, lcdat2) {
  A <- sum(cos(lcdat1)*cos(lcdat2))
  B <- sum(sin(lcdat1)*sin(lcdat2))
  C <- sum(cos(lcdat1)*sin(lcdat2))
  D <- sum(sin(lcdat1)*cos(lcdat2))
  E <- sum(cos(2*lcdat1)); F <- sum(sin(2*lcdat1))
  G <- sum(cos(2*lcdat2)); H <- sum(sin(2*lcdat2))
  n <- length(lcdat1)
  denom <- sqrt(((n*n)-(E*E)-(F*F))*((n*n)-(G*G)-(H*H)))
  rFLVal <- 4*((A*B)-(C*D))/denom
  return(rFLVal)
}
