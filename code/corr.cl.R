#' @title Correlation between angular and linear variable
#'
#' @lvar Linear variable
#' @cvar angular variable
R2xtCorrCoeff <- function(lvar, cvar) {
  rxc <- cor(lvar, cos(cvar))
  rxs <- cor(lvar, sin(cvar))
  rcs <- cor(cos(cvar), sin(cvar))
  R2xtVal <- ((rxc*rxc)+rxs*rxs- (2*rxc*rxs*rcs))/(1-rcs*rcs)
  return(R2xtVal)
}

#' @title Test for signfiicance of dorrelation between angular and linear variable
#'
#' @lvar Linear variable
#' @cvar angular variable
R2xtIndTestRand <- function(lvar, cvar,NR) {
  R2xtObs <- R2xtCorrCoeff(lvar, cvar)
  nxtrm <- 1
  for (r in 1:NR) {
    lvarRand <- sample(lvar)
    R2xtRand <- R2xtCorrCoeff(lvarRand, cvar)
    if (R2xtRand >= R2xtObs) { nxtrm <- nxtrm + 1} }
  pval <- nxtrm/(NR+1)
  return(c(R2xtObs, pval))
}
