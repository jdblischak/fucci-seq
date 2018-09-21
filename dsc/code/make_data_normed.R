#' @title Transform raw counts to log2cpm
#'
#' @param eset eset objected of final, filtered data containing counts
#'
#' @import parallel
#'
#' @export
#'
make_data_normed < - function(log2cpm, ncores=1) {

  log2cpm_normed <- do.call(rbind, mclapply(1:nrow(log2cpm), function(g) {
    yy <- log2cpm[g,]
    is.zero <- which(yy == 0)
    qq.map <- qqnorm(yy)
    yy.qq <- qq.map$x
    yy.qq[is.zero] <- sample(qq.map$x[is.zero])
    return(yy.qq)
  }, mc.cores=ncores) )
  rownames(log2cpm_normed) <- rownames(log2cpm)

  return(log2cpm_normed)
}
