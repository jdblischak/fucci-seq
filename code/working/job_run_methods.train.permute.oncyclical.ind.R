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
fold <- as.numeric(args[3])
ind <- args[4]


dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"code/working/run_methods.R"))
source(file.path(dir,"peco/R/fit.trendfilter.generic.R"))

data_permute_training <- readRDS(file=file.path(dir,
                            paste0("data/results/ind_",ind, "_data_permute_training.rds")))
data_training <- readRDS(file=file.path(dir,
                          paste0("data/results/ind_",ind, "_data_training.rds")))


seurat.genes <- readLines(con = file.path(dir,
                  "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])

cyclical_genes <- readRDS(file=file.path(dir, paste0(
  paste0("data/results/triple_ind_",ind,"_data_training_cyclical_genes.fold.",fold,".rds"))))


which_genes <- rownames(cyclical_genes)[order(cyclical_genes$pve,
                                              decreasing = T)[1:ngenes]]

fold_indices <- readRDS(file=file.path(dir,
                             paste0("data/results/ind_",ind,"_fold_indices.rds")))

print(ind)
print(fold)

Y_train_normed_fold_topX <- data_training$log2cpm.quant.nonvalid[
      which(rownames(data_training$log2cpm.quant.nonvalid) %in% which_genes),
      fold_indices[[fold]]$train]
theta_train_fold <- data_training$theta.nonvalid[fold_indices[[fold]]$train]

fit.train <- cycle.npreg.insample(Y = Y_train_normed_fold_topX,
                                  theta = theta_train_fold,
                                  polyorder=2,
                                  ncores=ncores,
                                  method.trend="trendfilter")

Y_test_normed_fold_topX <- data_permute_training$log2cpm.quant.nonvalid[
  which(rownames(data_training$log2cpm.quant.nonvalid) %in% which_genes),
  fold_indices[[fold]]$test]
Y_test_fold <- data_permute_training$log2cpm.nonvalid[,fold_indices[[fold]]$test]
theta_test_fold <- data_permute_training$theta.nonvalid[fold_indices[[fold]]$test]
pdata_test_fold <- data_permute_training$pdata.nonvalid[fold_indices[[fold]]$test,]

fit.test <- run_methods(Y_test=Y_test_fold,
                     Y_test_normed=Y_test_normed_fold_topX,
                     theta_test=theta_test_fold,
                     training_model=fit.train,
                     seurat.genes=seurat.genes,
                     pdata_test=pdata_test_fold,
                     fdata=data_training$fdata,
                     ncores=ncores, maxiter=30)

out <- list(fit.train=fit.train,
            fit.test=fit.test)
#}
#names(fits) <- paste0("fold.", 1:length(fold_indices))
saveRDS(out,
        file=file.path(dir,
                       paste0("data/results/ind_",ind,"_results_train_permute_oncyclical.fold.",fold,".top",ngenes,".rds")))



