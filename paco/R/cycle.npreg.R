#' @title Initialize cell time estimates for the nonparametric approach
#'
#' @param Y gene by sample matrix
#' @param nbins Y is grouped to nbins between 0 to 2pi.
#' @export
initialize_cell_times <- function(Y, nbins=NULL) {
  library(circular)
  pc_res <- prcomp(t(Y), scale = TRUE)
  theta <- coord2rad(cbind(pc_res$x[,1], pc_res$x[,2]))
  theta <- as.numeric(theta)
  names(theta) <- colnames(Y)

  if (is.null(nbins)) {
    return(theta)
  } else {

    bin.count <- c(1:nbins)
    bin.length <- (2 * pi)/nbins
    theta.assign <- c(1:length(theta))
    mids <- seq(bin.length/2, 2 * pi - pi/nbins, length = nbins)
    for (i in 1:nbins) {
      ii.bin <- theta <= i * bin.length & theta > (i - 1) * bin.length
      bin.count[i] <- sum(ii.bin)
      theta.assign[ii.bin] <- mids[i]
    }
    return(theta.assign)
  }
}


#' @title Compute estimates in one iteration
#'
#' @param Y gene by sample log2CPM gene expression matrix
#' @param cell_times cell time estimate to be used in estimating cyclical patterns
#'    of gene exprssion levels for each gene
#' @param ncores number of cores to be used
#'
#' @export
cycle.npreg.iter <- function(Y, cell_times, ncores, ...) {

  G <- nrow(Y)
  N <- ncol(Y)
  cell_times_ordered <- cell_times[order(cell_times)]
  Y_ordered <- Y[,order(cell_times)]

  # estimate expected mean expression and variance
  # given cell times for every gene
  fit <- mclapply(1:G, function(g) {
    y_g <- Y_ordered[g,]
    fit_g <- fit.trendfilter.generic(yy=y_g)
    mu_g <- fit_g$trend.yy
    sigma_g <- sqrt(sum((y_g-mu_g)^2)/N)
    lambda_g <- fit_g$lambda

    list(y_g =y_g,
         mu_g=mu_g,
         sigma_g=sigma_g,
         lambda_g=lambda_g)
  }, mc.cores = ncores)

  # for every cell, compute the data likelihood
  # summing log-likelihood at given cell times across genes
  # compute expected cell times for each bin
  cell_times_ordered_bins <- unique(cell_times_ordered)
  cell_times_ordered_bins <- cell_times_ordered_bins[order(cell_times_ordered_bins)]
  nbins <- length(cell_times_ordered_bins)
  mu_g_bins <- do.call(rbind, lapply(fit, function(x) {
    tmp <- aggregate(x$mu_g~factor(cell_times_ordered), FUN=mean)
    return(tmp[,2])
  }))
  colnames(mu_g_bins) <- cell_times_ordered_bins

  # compute the log-likelihood for every cell
  loglik_per_cell_by_celltimes <- matrix(0, N, nbins)

  for (n in 1:N) {
    for (g in 1:G)
      # compute the log-likelihood for every cell and gene g
      loglik_per_cell <- do.call(rbind, lapply(1:G, function(g) {
        with(fit[[g]], dnorm(y_g[n], mu_g_bins, sigma_g, log = TRUE))
      }))
    loglik_per_cell_by_celltimes[n,] <- colSums(loglik_per_cell)
  }

  # compute likelihood after initialization
  loglik_max_per_cell <- sapply(1:N, function(n) {
    ll <- loglik_per_cell_by_celltimes[n,]
    ll[which.max(ll)]
  })

  loglik_est <- sum(loglik_max_per_cell)
  cell_times_est <- sapply(1:nrow(loglik_per_cell_by_celltimes), function(n) {
    ll_n <- loglik_per_cell_by_celltimes[n,]
    cell_times_ordered_bins[which.max(ll_n)] })
  cell_times_est <- as.numeric(cell_times_est)
  names(cell_times_est) <- colnames(Y_ordered)

  lambda <- sapply(fit, "[[", "lambda_g")
  names(lambda) <- rownames(Y_ordered)

  sigma_est <- sapply(fit, "[[", "sigma_g")
  names(sigma_est) <- rownames(Y_ordered)

  mu_est <- do.call(rbind, lapply(fit, "[[", "mu_g"))
  colnames(mu_est) <- colnames(Y_ordered)
  rownames(mu_est) <- rownames(Y_ordered)

  return(list(Y = Y_ordered,
              mu_est = mu_est,
              sigma_est = sigma_est,
              loglik_est = loglik_est,
              cell_times_est = cell_times_est,
              lambda = lambda))
}





#' @title Estimate cell cycle ordering
cycle.npreg <- function(Y, nbins, tol = .001, ncores=12, verbose = TRUE, maxiter = 20,...) {

  cell_times_initial <- initialize_cell_times(Y, nbins=100)

  # order data by initial cell times
  G <- nrow(Y)
  N <- ncol(Y)
  Y_ordered <- Y[,order(cell_times_initial)]
  cell_times_previous <- cell_times_initial[order(cell_times_initial)]

  # initialize converge check
  loglik_previous <- .Machine$double.xmax
  eps <- tol+1
  iter <- 1

  # save loglik
  loglik <- c()
  lambda <- rep(0, nrow(Y_ordered))

  while (TRUE) {

    estimates <- cycle.npreg.iter(Y=Y_ordered, cell_times = cell_times_previous,
                                  ncores = ncores)

    loglik_iter <- estimates$loglik_est
    cell_times_iter <- estimates$cell_times_est
    lambda_iter <- estimates$lambda

    loglik <- c(loglik, loglik_iter)
    lambda <- cbind(lambda, lambda_iter)

    eps <- abs(loglik_iter - loglik_previous)/abs(loglik_previous)
    if (verbose) {
      cat("After iter", iter, "; log-likelihood = ", loglik_iter,
          "; eps = ", eps, "\n") }

    # loop out if converged
    if (!(eps > tol & iter < maxiter)) break
    iter <- iter + 1

    loglik_previous <- loglik_iter
    cell_times_previous <- cell_times_iter
  }

  lambda <- lambda[,-1]
  colnames(lambda) <- paste0("lambda.",c(1:iter))

  out <- list(Y_ordered=Y_ordered,
              mu_est=estimates$mu_est,
              sigma_est=estimates$sigma_est,
              loglik_est=loglik_iter,
              cell_times_est=cell_times_iter,
              loglik=loglik,
              iter=iter,
              lambda=lambda)
}
