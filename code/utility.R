
#' @title Print prime numbers
#' 
#' @param nprimes number of prime numbers needed (default starting from 1)
#' @param prime.start start search for prime numbers from this value
#' 
#' @export
primes <- function(nprimes, prime.start=1) {
  is.prime <- function(n) n == 2L || all(n %% 2L:max(2,floor(sqrt(n))) != 0)
  out <- c()
  num <- prime.start
  
  while(length(out)<nprimes) {
    if (is.prime(num) ) { 
      out <- c(out, num) } else {
      out <- out
      }
    num <- num+1
  }
  return(out)  
}

#' @title partition samples by specified sizes of training and testing data
#' 
#' @param y data vector
#' @param runs number of times that data vector will be partitioned
#' @param nsize_each number of samples in each partition
#' 
#' @export
partitionSamples.n <- function(y, runs, nsize.each = NULL) {
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





