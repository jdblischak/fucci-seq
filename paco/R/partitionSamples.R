#' @title partition samples by specified sizes of training and testing data
#'
#' @param y data vector
#' @param runs number of times that data vector will be partitioned
#' @param nsize_each number of samples in each partition
#'
#' @export
partitionSamples <- function(y, runs, nsize.each = NULL) {
  dat <- data.frame(y)
  dat$index <- c(1:length(y))
  nseeds <- runs*length(nsize.each)
  seeds <- matrix(primes(nseeds), nrow=runs, ncol=length(nsize.each))

  out <- lapply(1:runs, function(r) {
    lapply(1:length(nsize.each), function(i) {
      set.seed(seeds[r,i])
      indices <- sample(dat$index, size=nsize.each[i])
      indices <- sort(indices)
      dat$y[indices]
    })
  })
  return(list(partitions=out,
              seeds=seeds))
}
