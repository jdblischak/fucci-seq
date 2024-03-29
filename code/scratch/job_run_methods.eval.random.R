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

# get genes not in the top X list
set.seed(99)
ii_genes_include <- sample((ngenes+1):nrow(data_training$log2cpm.nonvalid), ngenes)

# make prediction parameters
Y_train_randomX <- with(data_training,
    log2cpm.quant.nonvalid[
      rownames(log2cpm.quant.nonvalid) %in% rownames(log2cpm.quant.nonvalid)[ii_genes_include], ])

training_randomX <- cycle.npreg.insample(Y = Y_train_randomX,
                                        theta = data_training$theta.nonvalid,
                                        polyorder=2,
                                        ncores=ncores,
                                        method.trend="trendfilter")


seurat.genes <- readLines(con = file.path(dir,
                                          "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])

results_eval_randomX <- run_methods(Y_test=data_withheld$log2cpm.valid,
                 Y_test_normed=data_withheld$log2cpm.quant.valid,
                 theta_test=data_withheld$theta.valid,
                 training_model=training_randomX,
                 seurat.genes=seurat.genes,
                 pdata_test=data_withheld$pdata.valid,
                 fdata=data_withheld$fdata,
                 ncores=ncores, maxiter=30)

saveRDS(results_eval_randomX,
        file=file.path(dir, paste0("data/results/results_eval_random",ngenes,".rds")))

