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
ind <- as.character(args[3])


dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"code/working/run_methods.R"))

data_training <- readRDS(file.path(dir,paste0("data/results/ind_",
                                ind, "_data_training.rds")))
data_withheld <- readRDS(file.path(dir,paste0("data/results/ind_",
                                ind, "_data_withheld.rds")))

fits_all <- readRDS(file.path(dir,"output/npreg-trendfilter-quantile.Rmd/fit.quant.rds"))
genes_all <- names(fits_all)[order(sapply(fits_all,"[[",3), decreasing=T)]

# cyclical_genes <- readRDS(file=file.path(dir, paste0(
#   paste0("data/results/finalizing/ind_",ind,"_cyclical_genes.rds"))))
# which_genes <- rownames(cyclical_genes)[order(cyclical_genes$pve,
#                                               decreasing = T)[1:ngenes]]

which_genes <- genes_all[1:ngenes]

seurat.genes <- readLines(con = file.path(dir,
                                          "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])



print(ind)
print(ngenes)

Y_train_normed_fold_topX <- data_training$log2cpm.quant.nonvalid[
  which(rownames(data_training$log2cpm.quant.nonvalid) %in% which_genes), ]
theta_train_fold <- data_training$theta.nonvalid

fit.train <- cycle.npreg.insample(Y = Y_train_normed_fold_topX,
                                  theta = theta_train_fold,
                                  polyorder=2,
                                  ncores=ncores,
                                  method.trend="trendfilter")

Y_test_normed_topX <- data_withheld$log2cpm.quant.valid[
  which(rownames(data_withheld$log2cpm.quant.valid) %in% which_genes),]
Y_test <- data_withheld$log2cpm.valid
theta_test <- data_withheld$theta.valid
pdata_test <- data_withheld$pdata.valid

fit.test <- run_methods(Y_test=Y_test,
                        Y_test_normed=Y_test_normed_topX,
                        theta_test=theta_test,
                        training_model=fit.train,
                        seurat.genes=seurat.genes,
                        pdata_test=pdata_test,
                        fdata=data_training$fdata,
                        ncores=ncores, maxiter=30)

out <- list(fit.train=fit.train,
            fit.test=fit.test)
#}
#names(fits) <- paste0("fold.", 1:length(fold_indices))
saveRDS(out,
        file=file.path(dir,
                       paste0("data/results/finalizing/ind_",ind,"_results_overallcyclical.top",ngenes,".rds")))





