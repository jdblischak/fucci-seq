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
