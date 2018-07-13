#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


ncores <- as.numeric(args[1])
nlist <- as.numeric(args[2])


dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"code/working/run_methods.R"))

data_training <- readRDS(file=file.path(dir, "data/results/data_training.rds"))
data_withheld <-readRDS(file=file.path(dir, "data/results/data_withheld.rds"))

data_permute_training <- readRDS(file=file.path(dir, "data/results/data_permute_training.rds"))
data_permute_withheld <-readRDS(file=file.path(dir, "data/results/data_permute_withheld.rds"))

sig.genes <- readRDS(file=file.path(dir,
    "data/results/results_topgenes.rds"))

sig.genes_top <- sig.genes[[nlist]]

num_genes <- length(sig.genes_top)

# make prediction parameters
Y_train_topX <- data_training$log2cpm.quant.nonvalid[
  rownames(data_training$log2cpm.quant.nonvalid) %in% sig.genes_top, ]
training_topX <- cycle.npreg.insample(Y = Y_train_topX,
                                        theta = data_training$theta.nonvalid,
                                        polyorder=2,
                                        ncores=15,
                                        method.trend="trendfilter")

# Y_train_seurat <- data_training$log2cpm.quant.nonvalid[
#   rownames(data_training$log2cpm.quant.nonvalid) %in% seurat.genes, ]
# training_seurat <- cycle.npreg.insample(Y = Y_train_topX,
#                                       theta = data_training$theta.nonvalid,
#                                       polyorder=2,
#                                       ncores=15,
#                                       method.trend="trendfilter")
seurat.genes <- readLines(con = file.path(dir,
                                          "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])

results_eval_topX <- run_methods(Y_test=data_permute_withheld$log2cpm.valid,
                 Y_test_normed=data_permute_withheld$log2cpm.quant.valid,
                 theta_test=data_permute_withheld$theta.valid,
                 training_model=training_topX,
                 seurat.genes=seurat.genes,
                 pdata_test=data_permute_withheld$pdata.valid,
                 fdata=data_withheld$fdata,
                 ncores=ncores, maxiter=30)

saveRDS(results_eval_topX,
        file=file.path(dir,
                       paste0("data/results/results_eval_permute_oncyclical.top",num_genes,".rds")))

