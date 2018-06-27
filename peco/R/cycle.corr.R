#' @title Compute correlations between rank orders of the cell times
#'
#' @param y1 numeric vector of length N with values betwen 0 to 2pi
#' @param y2 numeric vector of length N with values betwen 0 to 2pi
#'
#' @export
rotation <- function(y1, y2, ...) {
  theta_grid <- y2
#  theta_grid <- theta_grid[1:2]
  diff_time_sum <- c()
  for(which_cutpoint in 1:length(theta_grid)) {
    cutpoint <- theta_grid[which_cutpoint]
    y2shift <- rep(0, length(y2))
    for (c in 1:length(y2)) {
      if (y2[c] >= cutpoint) {y2shift[c] <- y2[c]-2*pi}
      else {
        y2shift[c] <- y2[c]
      }
    }
    if (min(y2shift) < 0) {y2shift <- y2shift-min(y2shift)}
    diff_time_sum <- c(diff_time_sum,
                       sum(pmin(abs(y1-y2shift), abs(y1-(2*pi-y2shift)))^2))
  }

  cutpoint_choose <- theta_grid[which.min(diff_time_sum)]
  y2shift <- rep(0, length(y2))
  for (c in 1:length(y2)) {
    if (y2[c] >= cutpoint_choose) {y2shift[c] <- y2[c]-2*pi}
    else {
      y2shift[c] <- y2[c]
    }
  }
  if (min(y2shift) < 0) {y2shift <- y2shift-min(y2shift)}
  names(y2shift) <- names(y2)

  return(list(y1=y1,
              y2=y2,
              y2shift=y2shift,
              cutpoint_choose=cutpoint_choose))
}


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

#' @title compute statistical significance of Fisher-Lee correlation coefficient
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#' @param NR number of permutations
#'
#' @references Pewsey et al. Circular statistics in R
#'
#' @export
rFLIndTestRand <- function(lcdat1, lcdat2, NR) {
  rFLObs <- rFLCorrCoeff(lcdat1, lcdat2)
  nxtrm <- 1
  for (r in 1:NR) {
    lcdat1Rand <- sample(lcdat1)
    rFLRand <- rFLCorrCoeff(lcdat1Rand, lcdat2)
    if (abs(rFLRand) >= abs(rFLObs)) {nxtrm <- nxtrm + 1} }
  pval <- nxtrm/(NR+1); return(c(rFLObs,pval))
}
