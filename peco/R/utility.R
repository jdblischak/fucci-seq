#' @title compute proportation of variable explained
#'
#' @export

get.pve <- function(yy, method="trendfilter") {
  if (method=="trendfilter") {
    fit <- fit.trendfilter.generic(yy, polyorder=3)
    trend.yy <- fit$trend.yy
    pve <- 1-var(yy-trend.yy)/var(yy)
  }
  return(pve)
}
