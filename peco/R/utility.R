#' @title compute proportation of variable explained
#'
#' @export

get.pve <- function(yy, method="trendfilter", polyorder=2) {
  if (method=="trendfilter") {
    fit <- fit.trendfilter.generic(yy, polyorder=polyorder)
    trend.yy <- fit$trend.yy
    ss_total <- (length(yy)-1)*var(yy)
    ss_error <- (length(yy)-1)*var(yy-trend.yy)

    pve <- 1-ss_error/ss_total
    pval <- var.test(yy, yy-trend.yy, alternative="two.sided")$p.value
  }
  return(list(pve=pve, pval=pval))
}


#' @title compute variance explained for discrete cell time
#'
#' @param yy continous measurement
#' @param xx categorical variable
#'
#' @export
get.aov <- function(yy, xx) {
  nclass <- length(unique(xx))
  xx <- as.numeric(xx)

  ss_total <- sum((yy-mean(yy))^2)
  ss_within <- sum(sapply(1:nclass, function(i) {
    rr <- yy[xx==i];
    sum((rr-mean(rr))^2) }))

  pve <- 1-(ss_within/ss_total)
  ss_between <- ss_total-ss_within

  df_total <- length(xx)-1
  df_between <- nclass-1
  f.obs <- (ss_between/df_between)/(ss_total/df_total)
  pval <- pf(f.obs, df1=df_between, df2=df_total, lower.tail=F)

  return(list(pve=pve, pval=pval))
}



#' #' @title difference in circular rank
#' #'
#' #' @export
#' get.diffRank <- function(pred_time, ref_time) {
#'
#'   n <- length(ref_time)
#'   crank1 <- rank(ref_time)*2*pi/n
#'   crank2 <- rank(pred_time)*2*pi/n
#'
#'   pmin(abs(crank2-crank1),
#'        abs(crank2-(2*pi-crank2)))
#'
#'     fit <- fit.trendfilter.generic(yy, polyorder=3)
#'     trend.yy <- fit$trend.yy
#'     pve <- 1-var(yy-trend.yy)/var(yy)
#'   }
#'   return(pve)
#' }
