training_predict_test <- function(log2cpm_normed, theta, pdata, which_training_sample,
                                  which_test_sample, which_genes, ncores) {

  # get cyclical expression patterns from the training samples
  log2cpm_normed_train_topX <- log2cpm[which_genes, which_training_samples]
  theta_train <- theta[which_training_sample]

  fit.train <- cycle.npreg.insample(Y = log2cpm_normed_train_topX,
                                    theta = theta_train,
                                    polyorder=2,
                                    ncores=ncores,
                                    method.trend="trendfilter")

  # predict cell times
  log2cpm_normed_test_topX <- log2cpm[which_genes, which_test_samples]
  theta_test <- theta[which_test_samples]
  pdata_test <- pdata[which_test_samples,]

  fit.test <- cycle.npreg.outsample(Y_test=log2cpm_normed_test_topX,
                                          sigma_est=fit.train$sigma_est,
                                          funs_est=fit.train$funs_est,
                                          theta_prior=fit.train$theta,
                                          method.grid = "uniform",
                                          method.trend="trendfilter",
                                          polyorder=2,
                                          ncores=ncores)

  out <- fit.test

  # tidy up the results
  out$ref_time <- theta_test
  out$pred_time <- with(out, cell_times_est[match(names(ref_time),
                                                  names(cell_times_est))])
  out$pred_time_shift <- with(out, rotation(ref_time, pred_time))
  out$diff_time <- with(out, circ_dist(pred_time_shift, ref_time))
  out$dapi <- pdata_test$dapi.median.log10sum.adjust[
    match(names(theta_test), rownames(pdata_test))]
  out$gfp <- pdata_test$gfp.median.log10sum.adjust[
    match(names(theta_test), rownames(pdata_test))]
  out$rfp <- pdata_test$rfp.median.log10sum.adjust[
    match(names(theta_test), rownames(pdata_test))]

  out <- list(fit.train=fit.train,
              fit.test=fit.test)
}







