validation <- function(log2cpm_normed, log2cpm,
                       fdata, which_genes,
                       which_training_sample,
                       which_validation_sample,
                       seurat.genes) {

  # supervised method -----------------------------------------------

  # get cyclical expression patterns from the training samples
  log2cpm_normed_train_topX <- log2cpm_normed[which_genes, which_training_samples]
  theta_train <- theta[which_training_sample]

  fit.train <- cycle.npreg.insample(Y = log2cpm_normed_train_topX,
                                    theta = theta_train,
                                    polyorder=2,
                                    ncores=ncores,
                                    method.trend="trendfilter")

  # predict cell times
  log2cpm_normed_test_topX <- log2cpm_normed[which_genes, which_validation_sample]
  theta_test <- theta[which_validation_sample]
  pdata_test <- pdata[which_validation_sample,]

  fit.supervised <- cycle.npreg.outsample(Y_test=log2cpm_normed_test_topX,
                                    sigma_est=fit.train$sigma_est,
                                    funs_est=fit.train$funs_est,
                                    theta_prior=fit.train$theta,
                                    method.grid = "uniform",
                                    method.trend="trendfilter",
                                    polyorder=2,
                                    ncores=ncores)

  # Seurat -----------------------------------------------
  log2cpm_test <- log2cpm_[, which_test_samples]
  symbols <- fdata$name[match(rownames(Y_test), rownames(fdata))]
  rownames(log2cpm_test) <- symbols

  fit.seurat <- run_seurat(Y=log2cpm_test,
                           s.genes=seurat.genes$s.genes,
                           g2m.genes=seurat.genes$g2m.genes,
                           n.bin=25,
                           seed.use=1, random.seed=1)

  seurat.pca <- prcomp(cbind(fit.seurat$G2M, fit.seurat$S), scale=TRUE)
  seurat.cell_times_est <- as.numeric(coord2rad(cbind(seurat.pca$x[,1],seurat.pca$x[,2])))
  names(seurat.cell_times_est) <- colnames(Y_test_seurat)
  fit.seurat$cell_times_est <- seurat.cell_times_est

  fit.seurat <- as.list(fit.seurat)

  out <- list(fit.supervised=fit.supervised,
              fit.seurat=fit.seurat)

  for (i in 1:length(out)) {
    out[[i]]$ref_time <- theta_test
    out[[i]]$pred_time <- with(out[[i]], cell_times_est[match(names(ref_time),
                                                    names(cell_times_est))])
    out[[i]]$pred_time_shift <- with(out[[i]], rotation(ref_time, pred_time))
    out[[i]]$diff_time <- with(out[[i]], circ_dist(pred_time_shift, ref_time))
    out[[i]]$dapi <- pdata_test$dapi.median.log10sum.adjust[
      match(names(theta_test), rownames(pdata_test))]
    out[[i]]$gfp <- pdata_test$gfp.median.log10sum.adjust[
      match(names(theta_test), rownames(pdata_test))]
    out[[i]]$rfp <- pdata_test$rfp.median.log10sum.adjust[
      match(names(theta_test), rownames(pdata_test))]
  }

  return(out)
}



