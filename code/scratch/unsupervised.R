#' @title unsupervised method for learning cell times
#'
#' @export
cycle.npreg.unsupervised <- function(Y, theta,
                                 ncores=12,
                                 method.trend=c("npcirc.nw", "npcirc.ll", "trendfilter"),
                                 grids=100,
                                 maxiter=30, verbose=TRUE, tol=1,
                                 ...) {

  # order data by initial cell times
  G <- nrow(Y)
  N <- ncol(Y)
  theta_ordered_initial <- theta[order(theta)]
  Y_ordered <- Y[,order(theta)]

  # initialize mu and sigma
  initial_mstep <- cycle.npreg.mstep(Y = Y_ordered,
                                     theta = theta_ordered_initial,
                                     method.trend=method.trend,
                                     ncores = ncores)

  # compute log-likelihood under initial mu and sigma
  initial_estep <- cycle.npreg.loglik(Y = initial_mstep$Y,
                                      sigma_est = initial_mstep$sigma_est,
                                      funs_est=initial_mstep$funs,
                                      method.grid="pca",
                                      method.type="unsupervised")

  sigma_est_previous <- initial_mstep$sigma_est
  funs_est_previous <- initial_mstep$funs
  Y_previous <- initial_mstep$Y
  cell_times_previous <- initial_estep$cell_times_est
  loglik_previous <- initial_estep$loglik_est

  iter <- 0
  while(TRUE) {
    current_mstep <- cycle.npreg.mstep(Y = Y_previous,
                                       theta = cell_times_previous,
                                       method.trend = method.trend,
                                       ncores = ncores)
    current_loglik <- cycle.npreg.loglik(Y = current_mstep$Y,
                                sigma_est = current_mstep$sigma_est,
                                funs_est = current_mstep$funs,
                                method.grid="pca",
                                method.type="unsupervised")

    sigma_est_current <- current_mstep$sigma_est
    funs_est_current <- current_mstep$funs
    Y_current <- current_mstep$Y
    cell_times_current <- current_loglik$cell_times_est
    loglik_current <- current_loglik$loglik_est

    if (verbose) message("log-likelihood:", loglik_current)
    eps <- loglik_current - loglik_previous

    # loop out if converged
    if (!(eps > tol & iter < maxiter)) break

    iter <- iter + 1
    loglik_previous <- loglik_current
    Y_previous <- Y_current
    sigma_est_previous <- sigma_est_current
    cell_times_previous <- cell_times_current
    funs_est_previous <- funs_est_current
    }

  out <- list(Y_ordered=Y_current,
              loglik_est=loglik_current,
              cell_times_est=cell_times_current,
             # mu_est=mu_est_current,
              sigma_est=sigma_est_current,
              funs_est=funs_est_current)
  return(out)
}

