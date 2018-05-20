#' @title Initialize cell time estimates for the nonparametric approach
#'
#' @description Estimate cell times using PC1 and PC2. If \code{nbins} is specified,
#'   then group the estimated cell times to \code{nbins} groups.
#'
#' @param Y gene by sample matrix
#' @param nbins Y is grouped to nbins between 0 to 2pi.
#'
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


#' @title log-likelihood of nonparametric smoothing
#'
#' @param Y gene by sample expression matrix
#' @param mu_est gene by sample matrix of expected mean
#' @param sigma_est vector of standard errors for each gene
#'
#' @export
cycle.npreg.loglik <- function(Y, theta, mu_est, sigma_est) {

  nbins <- ncol(mu_est)
  N <- ncol(Y)
  G <- nrow(Y)
  loglik_per_cell_by_celltimes <- matrix(0, N, nbins)

  for (n in 1:N) {
    # for each cell, sum up the loglikelihood for each gene
    # at the observed cell times
    loglik_per_cell <- do.call(rbind, lapply(1:G, function(g) {
      dnorm(Y[g,n], mu_est[g,], sigma_est[g], log = TRUE)
    }))
    loglik_per_cell <- colSums(loglik_per_cell)

    loglik_per_cell_by_celltimes[n,] <- loglik_per_cell
  }

  # use multinomial sampling to assign samples
  prob_per_cell_by_celltimes <- sweep(exp(loglik_per_cell_by_celltimes), 1,
                                      FUN = "/",
                                      STATS = rowMax(exp(loglik_per_cell_by_celltimes)))

  cell_times_samp_ind <- sapply(1:N, function(n) {
    which(rmultinom(1,1,prob_per_cell_by_celltimes[n,])==1)
  })

  cell_times_est <- sapply(1:N, function(n) {
    theta[cell_times_samp_ind[n]]
  })
  names(cell_times_est) <- colnames(Y)

  # compute likelihood based on the selected cell times
  loglik_max_per_cell <- sapply(1:N, function(n) {
    ll <- loglik_per_cell_by_celltimes[n,]
    ll[cell_times_samp_ind[n]]
  })
  loglik_est <- sum(loglik_max_per_cell)

  return(list(loglik_est=loglik_est,
              cell_times_est=cell_times_est))
}




#' @title Estimate parameters for the cyclial ordering using nonparametric smoothing
#'
#' @description Conditioned on observed cell times, estimate the cyclical trend
#'   for each gene, and compute the likelihood of the observed cell times.
#'
#' @param Y gene by sample expression matrix (log2CPM).
#' @param theta observed cellt times

cycle.npreg.mstep <- function(Y, theta, ncores=12, ...) {

  G <- nrow(Y)
  N <- ncol(Y)
  Y_ordered <- Y[,order(theta)]
  theta_ordered <- theta[order(theta)]

  # for each gene, estimate the cyclical pattern of gene expression
  # conditioned on the given cell times
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

  lambda <- sapply(fit, "[[", "lambda_g")
  names(lambda) <- rownames(Y_ordered)

  sigma_est <- sapply(fit, "[[", "sigma_g")
  names(sigma_est) <- rownames(Y_ordered)

  mu_est <- do.call(rbind, lapply(fit, "[[", "mu_g"))
  colnames(mu_est) <- colnames(Y_ordered)
  rownames(mu_est) <- rownames(Y_ordered)

  return(list(Y = Y_ordered,
              theta = theta_ordered,
              mu_est = mu_est,
              sigma_est = sigma_est,
              lambda = lambda))
}





#' @title Compute estimates in one iteration
#'
#' @description Given the estimated cell times, gene mean and gene variation,
#'   (possibly from m-step), compute the expected cell times at midpoint of defined bins
#'
#' @param Y gene by sample log2CPM gene expression matrix
#' @param theta cell time estimate to be used in estimating cyclical patterns
#'    of gene exprssion levels for each gene
#' @param ncores number of cores to be used
#'
#' @export
cycle.npreg.estep <- function(Y, theta, mu_est, sigma_est,
                              nbins=100, ncores, ...) {

  G <- nrow(Y)
  N <- ncol(Y)
  theta_ordered <- theta[order(theta)]
  Y_ordered <- Y[,order(theta)]

  # split cell_times to nbins equaly sized bins
  while (TRUE) {
    breaks <- quantile(theta_ordered, prob=seq(0,1,1/nbins))

    if (sum(duplicated(breaks)) == 0) break

    nbins <- nbins-sum(duplicated(breaks))
    breaks <- quantile(theta_ordered, prob=seq(0,1,1/nbins))
  }

  theta_ordered_bins <- data.frame(theta_ordered=theta_ordered,
                                        bins=cut(theta_ordered,
                                                 breaks=breaks,
                                                 include.lowest = TRUE))
  # mean cell time in each bin
  theta_ordered_bins_mean <- aggregate(theta_ordered ~ bins,
                                            data=theta_ordered_bins,
                                            FUN=mean)

  # expected mean gene expression in each bin
  mu_g_bins <- do.call(rbind, lapply(1:G, function(g) {
    tmp <- aggregate(mu_est[g,]~theta_ordered_bins$bins,
                     FUN=mean)
    return(tmp[,2])
  }))
  colnames(mu_g_bins) <- as.character(theta_ordered_bins_mean$bins)
  rownames(mu_g_bins) <- rownames(mu_est)

  loglik_max_est <- cycle.npreg.loglik(Y = Y_ordered,
                                       theta = theta_ordered,
                                       mu_est = mu_g_bins,
                                       sigma_est = sigma_est)

  Y_ordered_update <- Y_ordered[,order(loglik_max_est$cell_times_est)]
  cell_times_est_update <- with(loglik_max_est, cell_times_est[order(cell_times_est)])

  return(list(Y = Y_ordered_update,
              cell_times_est = cell_times_est_update,
              nbins=length(breaks)-1))
}



#' @title Estimate cell cycle ordering in the current sample
#'
#' @export
cycle.npreg.insample <- function(Y, theta, nbins=200, ncores=12,...) {

  # order data by initial cell times
  G <- nrow(Y)
  N <- ncol(Y)
  theta_ordered_initial <- theta[order(theta)]
  Y_ordered <- Y[,order(theta)]

  # initialize mu and sigma
  initial <- cycle.npreg.mstep(Y = Y_ordered,
                              theta = theta_ordered_initial, ncores = 10)
  # compute expected cell time under initial mu and sigma
  estep <- cycle.npreg.estep(Y = initial$Y,
                             theta = initial$theta,
                             mu_est = initial$mu_est,
                             sigma_est = initial$sigma_est, nbins = nbins)
  # compute estimated mu and sigma given expected cell time
  mstep <- cycle.npreg.mstep(Y = estep$Y,
                             theta = estep$cell_times_est, ncores = 10)

  mstep.loglik <- cycle.npreg.loglik(Y = mstep$Y,
                                     theta = mstep$theta,
                                     mu_est=mstep$mu_est,
                                     sigma_est=mstep$sigma_est)$loglik_est
  out <- list(Y_ordered=estep$Y,
              cell_times_est=estep$cell_times_est,
              loglik_est=mstep.loglik,
              mu_est=mstep$mu_est,
              sigma_est=mstep$sigma_est)
}





#' @title Predict cell cycle ordering in the test samples
#'
#' @export
cycle.npreg.outsample <- function(Y_test, theta_est,
                                  mu_est, sigma_est, ncores=12,...) {

  # order data by initial cell times
  G <- nrow(Y_test)
  N <- ncol(Y_test)

  # compute expected cell time for the test samples
  # under mu and sigma estimated from the training samples
  pred <- cycle.npreg.loglik(Y = Y_test, theta = theta_est,
                             mu_est = mu_est, sigma_est = sigma_est)

  # compute estimated mu and sigma given expected cell time
  mstep <- cycle.npreg.mstep(Y = Y_test,
                             theta = pred$cell_times_est, ncores = 10)

  mstep.loglik <- cycle.npreg.loglik(Y = mstep$Y,
                                     theta = mstep$theta,
                                     mu_est=mstep$mu_est,
                                     sigma_est=mstep$sigma_est)$loglik_est
  out <- list(Y_ordered=mstep$Y,
              cell_times_est=mstep$cell_times_est,
              loglik_est=mstep.loglik,
              mu_est=mstep$mu_est,
              sigma_est=mstep$sigma_est)
}










