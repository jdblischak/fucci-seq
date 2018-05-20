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
