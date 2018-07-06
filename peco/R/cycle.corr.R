#' @title pairwise distance between two circular variables
#'
circ_dist <- function(y1,y2){
  pmin(abs(y2-y1), abs(2*pi-(abs(y2-y1))))
}


#' @title rotate circular variable y2 to minimize distance between y1 and y2
#'
#' @export
rotation <- function(ref_var,shift_var){

  df <- data.frame(flip=rep(c(1,-1), each=length(shift_var)),
                   shift = rep(shift_var, 2))

  for (i in 1:nrow(df)) {
    shift_var_tmp <- df$flip[i]*((shift_var+df$shift[i])%%(2*pi))
    df$dis[i] <- mean(circ_dist(ref_var, shift_var_tmp))
  }
  which_cutoff <- which.min(df$dis)

  shift_var_new <- df$flip[which_cutoff]*((shift_var+df$shift[which_cutoff])%%(2*pi))

  return(shift_var_new)
}



#' #' @title Circular data rotation
#' #'
#' #' @description For any pair of circular data (y1,y2), shift y2 until
#' #'   the distance between y1 and y2 is minimized
#' #'
#' #' @param y1 numeric vector of length N with values betwen 0 to 2pi
#' #' @param y2 numeric vector of length N with values betwen 0 to 2pi
#' #'
#' #' @export
#' rotation <- function(y1, y2, ...) {
#'   theta_grid <- y2
#'   #  theta_grid <- theta_grid[1:2]
#'   diff_time_sum <- c()
#'   for(which_cutpoint in 1:length(theta_grid)) {
#'     cutpoint <- theta_grid[which_cutpoint]
#'     y2shift <- rep(0, length(y2))
#'     for (c in 1:length(y2)) {
#'       if (y2[c] >= cutpoint) {y2shift[c] <- y2[c]-2*pi}
#'       else {
#'         y2shift[c] <- y2[c]
#'       }
#'     }
#'     if (min(y2shift) < 0) {y2shift <- y2shift-min(y2shift)}
#'     diff_time_sum <- c(diff_time_sum,
#'                        sum((pmin(abs(y1-y2shift), abs(y1-(2*pi-y2shift))))^2))
#'   }
#'
#'   cutpoint_choose <- theta_grid[which.min(diff_time_sum)]
#'   y2shift <- rep(0, length(y2))
#'   for (c in 1:length(y2)) {
#'     if (y2[c] >= cutpoint_choose) {y2shift[c] <- y2[c]-2*pi}
#'     else {
#'       y2shift[c] <- y2[c]
#'     }
#'   }
#'   if (min(y2shift) < 0) {y2shift <- y2shift-min(y2shift)}
#'   names(y2shift) <- names(y2)
#'
#'   return(list(y1=y1,
#'               y2=y2,
#'               y2shift=y2shift,
#'               cutpoint_choose=cutpoint_choose))
#' }


#' @title Circular-circular rank correlation coefficient
#'
#' @references Fisher NJ and Lee AJ (1983). A correlation coefficient for circular data. Biometrika, 70, 327-332.
#'
#' @export

rFLRank.CorrCoeff <- function(lcdat1, lcdat2) {

  # convert data to linear ranks
  n <- length(lcdat1)
  crank1 <- rank(lcdat1)*2*pi/n
  crank2 <- rank(lcdat2)*2*pi/n

  A <- sum(cos(crank1)*cos(crank2))
  B <- sum(sin(crank1)*sin(crank2))
  C <- sum(cos(crank1)*sin(crank2))
  D <- sum(sin(crank1)*cos(crank2))

  val <- ((A*B)-(C*D))*(4/(n^2))

  return(val)
}



#' @title Significance test for circular rank correlation, permutation-based
#'
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#' @param NR number of permutations
#'
#' @references Pewsey et al. Circular statistics in R
#'
#' @export
rFLRank.IndTestRand <- function(lcdat1, lcdat2, NR=10000) {

  n <- length(lcdat1)
  crank1 <- rank(lcdat1)*2*pi/n
  crank2 <- rank(lcdat2)*2*pi/n

  rFLObs <- rFLRank.CorrCoeff(lcdat1, lcdat2)
  nxtrm <- 1
  for (r in 1:NR) {
    crank1Rand <- sample(crank1)
    rFLRand <- rFLRank.CorrCoeff(crank1Rand, crank2)
    if (abs(rFLRand) >= abs(rFLObs)) {nxtrm <- nxtrm + 1} }
  pval <- nxtrm/(NR+1); return(c(rFLObs,pval))
}



#' @title Fisher-Lee correlation coefficient
#'
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#'
#' @references
#'    Fisher NJ and Lee AJ (1983). A correlation coefficient for circular data. Biometrika, 70, 327-332.
#'    Pewsey et al. Circular statistics in R
#'
#' @export
rFL.CorrCoeff <- function(lcdat1, lcdat2) {
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

#' @title statistical significance of Fisher-Lee correlation coefficient, permutation-based
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#' @param NR number of permutations
#'
#' @references
#'    Fisher NJ and Lee AJ (1983). A correlation coefficient for circular data. Biometrika, 70, 327-332.
#'    Pewsey et al. Circular statistics in R
#'
#' @export
rFL.IndTestRand <- function(lcdat1, lcdat2, NR=10000) {
  rFLObs <- rFL.CorrCoeff(lcdat1, lcdat2)
  nxtrm <- 1
  for (r in 1:NR) {
    lcdat1Rand <- sample(lcdat1)
    rFLRand <- rFL.CorrCoeff(lcdat1Rand, lcdat2)
    if (abs(rFLRand) >= abs(rFLObs)) {nxtrm <- nxtrm + 1} }
  pval <- nxtrm/(NR+1); return(c(rFLObs,pval))
}

