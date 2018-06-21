#' @title Initialize grid points for estimation
#'
#' @export
initialize_grids <- function(Y, grids=100,
                             method.grid=c("pca", "uniform"), ...) {

  len <- (2*pi)/(2*grids)
  theta_grids <- seq(len, (2*pi)-(len), length.out=grids)

  library(circular)
  if (method.grid=="pca") {
    pc_res <- prcomp(t(Y), scale = TRUE)
    grid_approx <- coord2rad(cbind(pc_res$x[,1], pc_res$x[,2]))
    grid_approx <- as.numeric(grid_approx)
    names(grid_approx) <- colnames(Y)

    theta_initial <- rep(0, length(grid_approx))
    theta_initial_ind <- rep(0, length(grid_approx))
    names(theta_initial) <- names(grid_approx)
    names(theta_initial_ind) <- names(grid_approx)

    for (i in 1:length(grid_approx)) {
      theta_initial_ind[i] <- which.min(pmin(abs(theta_grids-grid_approx[i]),
                                             abs(theta_grids-(2*pi-grid_approx[i]))))
      theta_initial <- theta_grids[theta_initial_ind]
    }
  }

  if (method.grid=="uniform") {
    theta_initial <- theta_grids
  }

  return(theta_initial)
}


#' #' @title Initialize cell time estimates for the nonparametric approach
#' #'
#' #' @description Estimate cell times using PC1 and PC2. If \code{nbins} is specified,
#' #'   then group the estimated cell times to \code{nbins} groups.
#' #'
#' #' @param Y gene by sample matrix
#' #' @param nbins Y is grouped to nbins between 0 to 2pi.
#' #'
#' #' @export
#' initialize_cell_times <- function(Y,
#'                                   method.initialize.theta=c("pca", "uniform"), ...) {
#'   library(circular)
#'   if (method.initialize.theta=="pca") {
#'     pc_res <- prcomp(t(Y), scale = TRUE)
#'     theta <- coord2rad(cbind(pc_res$x[,1], pc_res$x[,2]))
#'     theta <- as.numeric(theta)
#'     names(theta) <- colnames(Y)
#'   }
#'   if (method.initialize.theta=="uniform") {
#'     len <- 2*pi/ncol(Y)/2
#'     theta <- seq(len, 2*pi-len, length.out=nbins)
#'     names(theta) <- colnames(Y)
#'     return(theta)
#'   }
#'
#'   return(theta)
#' }


#' @title Estimate parameters for the cyclial ordering using nonparametric smoothing
#'
#' @description Conditioned on observed cell times, estimate the cyclical trend
#'   for each gene, and compute the likelihood of the observed cell times.
#'
#' @param Y gene by sample expression matrix (log2CPM).
#' @param theta observed cellt times
cycle.npreg.mstep <- function(Y, theta, method.trend=c("trendfilter",
                                                       "npcirc.nw",
                                                       "npcirc.ll"),
                              polyorder=3,
                              ncores=12, ...) {
      library(NPCirc)
      library(genlasso)
      library(assertthat)
      G <- nrow(Y)
      N <- ncol(Y)

      if (!assert_that(all.equal(names(theta), colnames(Y)))) {
        Y_ordered <- Y[,match(names(theta), colnames(Y))]
        ord <- order(theta)
        theta_ordered <- theta[ord]
        Y_ordered <- Y_ordered[,ord]
      } else {
        ord <- order(theta)
        theta_ordered <- theta[ord]
        Y_ordered <- Y[,ord]
      }

      # for each gene, estimate the cyclical pattern of gene expression
      # conditioned on the given cell times
      fit <- mclapply(1:G, function(g) {
        #    print(g)
        y_g <- Y_ordered[g,]

        if (method.trend=="npcirc.nw") {
          fit_g <- kern.reg.circ.lin(theta_ordered, y_g, method = "NW",
                                     t=theta_ordered)
          fun_g <- approxfun(x=as.numeric(fit_g$x), y=fit_g$y, rule=2)
          mu_g <- fun_g(theta_ordered)
        }
        if (method.trend=="npcirc.ll") {
          fit_g <- kern.reg.circ.lin(theta_ordered, y_g, method = "LL",
                                     t=theta_ordered)
          fun_g <- approxfun(x=as.numeric(fit_g$x), y=fit_g$y, rule=2)
          mu_g <- fun_g(theta_ordered)
        }
        if (method.trend=="trendfilter") {
          fit_g <- fit.trendfilter.generic(yy=y_g, polyorder = polyorder)
          fun_g <- approxfun(x=as.numeric(theta_ordered),
                             y=as.numeric(fit_g$trend.yy), rule=2)
          mu_g <- fit_g$trend.yy
        }
        sigma_g <- sqrt(sum((y_g-mu_g)^2)/N)

        list(y_g =y_g,
             mu_g=mu_g,
             sigma_g=sigma_g,
             fun_g=fun_g)
      }, mc.cores = ncores)

      sigma_est <- sapply(fit, "[[", "sigma_g")
      names(sigma_est) <- rownames(Y_ordered)

      mu_est <- do.call(rbind, lapply(fit, "[[", "mu_g"))
      colnames(mu_est) <- colnames(Y_ordered)
      rownames(mu_est) <- rownames(Y_ordered)

      funs <- sapply(fit, "[[", "fun_g")
      names(funs) <- rownames(Y_ordered)

      return(list(Y = Y_ordered,
                  theta = theta_ordered,
                  mu_est = mu_est,
                  sigma_est = sigma_est,
                  funs = funs))
    }



#' @title log-likelihood of nonparametric smoothing
#'
#' @param Y gene by sample expression matrix
#' @param mu_est gene by sample matrix of expected mean
#' @param sigma_est vector of standard errors for each gene
#'
#' @export
#' @title log-likelihood of nonparametric smoothing
#'
#' @param Y gene by sample expression matrix
#' @param mu_est gene by sample matrix of expected mean
#' @param sigma_est vector of standard errors for each gene
#'
#' @export
cycle.npreg.loglik <- function(Y, mu_est, sigma_est,
                               funs_est,
                               grids=100,
                               method.grid=c("pca", "uniform"), ...) {


  N <- ncol(Y)
  G <- nrow(Y)

  theta_choose <- initialize_grids(Y, grids=grids, method.grid=method.grid)
  loglik_per_cell_by_celltimes <- matrix(0, N, grids)


  for (n in 1:N) {
    # for each cell, sum up the loglikelihood for each gene
    # at the observed cell times
    loglik_per_cell <- do.call(rbind, lapply(1:G, function(g) {
      dnorm(Y[g,n], funs_est[[g]](theta_choose), sigma_est[g], log = TRUE)
    }))
    loglik_per_cell <- colSums(loglik_per_cell)

    loglik_per_cell_by_celltimes[n,] <- loglik_per_cell
  }

  # use max likelihood to assign samples
  prob_per_cell_by_celltimes <- matrix(0, N, grids)
  for (n in 1:N) {
    # print(n)
    maxll <- max(exp(loglik_per_cell_by_celltimes)[n,], na.rm=T)
    if (maxll == 0) {
      prob_per_cell_by_celltimes[n,] <- rep(0, grids)
    } else {
      prob_per_cell_by_celltimes[n,] <- exp(loglik_per_cell_by_celltimes)[n,]/maxll
    }
  }
  cell_times_samp_ind <- sapply(1:N, function(n) {
    if (max(prob_per_cell_by_celltimes[n,], na.rm=T)==0) {
      sample(1:grids, 1, replace=F)
    } else {
      which.max(prob_per_cell_by_celltimes[n,])
    }
  })
  cell_times_est <- sapply(1:N, function(n) {
    theta_choose[cell_times_samp_ind[n]]
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







#' @title Estimate cell cycle ordering in the current sample
#'
#' @param update T/F to update cell times
#'
#' @export
cycle.npreg.insample <- function(Y, theta,
                                 ncores=12,
                                 polyorder=3,
                                 method.trend=c("npcirc.nw", "npcirc.ll", "trendfilter"),
                                 ...) {

  # order data by initial cell times
  G <- nrow(Y)
  N <- ncol(Y)
  theta_ordered_initial <- theta[order(theta)]
  Y_ordered <- Y[,order(theta)]

  # initialize mu and sigma
  initial_mstep <- cycle.npreg.mstep(Y = Y_ordered,
                                     theta = theta_ordered_initial,
                                     polyorder=polyorder,
                                     method.trend=method.trend,
                                     ncores = ncores)

  out <- list(Y=Y,
              theta=theta,
              sigma_est=initial_mstep$sigma_est,
              funs_est=initial_mstep$funs)
  return(out)
}





#' @title Predict test-sample ordering using training lables (no update)
#'
#' @export
cycle.npreg.outsample <- function(Y_test,
                                  sigma_est,
                                  funs_est,
                                  method.trend=c("npcirc.nw", "npcirc.ll", "trendfilter"),
                                  polyorder=3,
                                  method.grid=c("pca", "uniform"),
                                  ncores=12,...) {

  # compute expected cell time for the test samples
  # under mu and sigma estimated from the training samples
  initial_loglik <- cycle.npreg.loglik(Y = Y_test,
                                       sigma_est = sigma_est,
                                       method.grid=method.grid,
                                       funs_est=funs_est)
  updated_estimates <- cycle.npreg.mstep(Y = Y_test,
                                     theta = initial_loglik$cell_times_est,
                                     method.trend = method.trend,
                                     polyorder=polyorder,
                                     ncores = ncores)

  out <- list(Y=Y_test,
              cell_times_est=initial_loglik$cell_times_est,
              loglik_est=initial_loglik$loglik_est,
              Y_reordered=updated_estimates$Y,
              cell_times_reordered=updated_estimates$theta,
              mu_reordered=updated_estimates$mu_est,
              sigma_reordered=updated_estimates$sigma_est,
              funs_reordered=updated_estimates$funs)
  return(out)
}




#' @title Predict test-sample ordering using training lables (no update)
#'
#' @export
cycle.npreg.unsupervised <- function(Y_test,
                                  sigma_est,
                                  funs_est,
                                  method.trend=c("npcirc.nw", "npcirc.ll", "trendfilter"),
                                  method.initialize.theta=c("pca", "uniform"),
                                  ncores=12,
                                  maxiter=10,
                                  tol=1, verbose=TRUE,...) {

  # compute expected cell time for the test samples
  # under mu and sigma estimated from the training samples
  initial_loglik <- cycle.npreg.loglik(Y = Y_test,
                                       sigma_est = sigma_est,
                                       method.initialize.theta=method.initialize.theta,
                                       funs_est=funs_est, outsample=T)
  updated_estimates <- cycle.npreg.mstep(Y = Y_test,
                                         theta = initial_loglik$cell_times_est,
                                         method.trend = method.trend,
                                         polyorder=polyorder,
                                         ncores = ncores)
  theta_update <- updated_estimates$theta
  names(theta_update) <- colnames(updated_estimates$Y)

  out <- list(Y_ordered=updated_estimates$Y,
              cell_times_update=theta_update,
              loglik_update=initial_loglik$loglik_est,
              mu_update=updated_estimates$mu_est,
              sigma_update=updated_estimates$sigma_est,
              funs_update=updated_estimates$funs,
              sigma_update=sigma_est,
              funs_update=funs_est)
  return(out)
}
