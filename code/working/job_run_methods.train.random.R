#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


ncores <- as.numeric(args[1])
ngenes <- as.numeric(args[2])


dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"code/working/run_methods.R"))


data_training <- readRDS(file=file.path(dir, "data/results/data_training.rds"))
data_withheld <-readRDS(file=file.path(dir, "data/results/data_withheld.rds"))

sig.genes <- readRDS(file=file.path(dir,
    "output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.476.rds"))

seurat.genes <- readLines(con = file.path(dir,
                                          "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])


fold_indices <- readRDS(file=file.path(dir, "data/results/fold_indices.rds"))

# get genes not in the top X list
set.seed(99)
ii_genes_include <- sample((ngenes+1):nrow(data_training$log2cpm.nonvalid), ngenes)


Y_normed_randomX <- data_training$log2cpm.quant.nonvalid[
  rownames(data_training$log2cpm.quant.nonvalid) %in% rownames(sig.genes)[ii_genes_include], ]
fdata_randomX <- data_training$fdata[rownames(data_training$fdata) %in% rownames(sig.genes)[ii_genes_include],]

# Y <- data_training$log2cpm.nonvalid[
#   rownames(data_training$log2cpm.nonvalid) %in% rownames(sig.genes)[1:ngenes], ]



fits <- vector("list", length(fold_indices))
for (run in 1:length(fold_indices)) {
  print(run)

  Y_train_random_fold <- Y_normed_randomX[,fold_indices[[run]]$train]
  theta_train_fold <- data_training$theta.nonvalid[fold_indices[[run]]$train]
  fit.train <- cycle.npreg.insample(Y = Y_train_random_fold,
                                    theta = theta_train_fold,
                                    polyorder=2,
                                    ncores=ncores,
                                    method.trend="trendfilter")

  Y_test_fold <- data_training$log2cpm.nonvalid[,fold_indices[[run]]$test]
  Y_test_normed_fold <- Y_normed_randomX[,fold_indices[[run]]$test]
  theta_test_fold <- data_training$theta.nonvalid[fold_indices[[run]]$test]
  pdata_test_fold <- data_training$pdata.nonvalid[fold_indices[[run]]$test,]

  fits[[run]] <- run_methods(Y_test=Y_test_fold,
                             Y_test_normed=Y_test_normed_fold,
                             theta_test=theta_test_fold,
                             training_model=fit.train,
                             seurat.genes=seurat.genes,
                             pdata_test=pdata_test_fold,
                             fdata=fdata_randomX,
                             ncores=ncores, maxiter=30)
}
names(fits) <- paste0("fold.", 1:length(fold_indices))

saveRDS(fits,
        file=file.path(dir, paste0("data/results/results_train_random",ngenes,".rds")))





