#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


ncores <- as.numeric(args[1])
fold <- as.numeric(args[2])
ind <- as.character(args[3])

dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"peco/R/fit.trendfilter.generic.R"))
source(file.path(dir,"peco/R/get.cyclical.R"))


data_training <- readRDS(file=file.path(dir,
                         paste0("data/results/triple_ind_",ind,"_data_training_alge.rds")))

fold_indices <- readRDS(file=file.path(dir,
                         paste0("data/results/ind_", ind, "_fold_indices.rds")))

print(ind)

print(fold)

Y_train_normed_fold <- data_training$log2cpm.quant.nonvalid[ ,fold_indices[[fold]]$train]
theta_train_fold <- data_training$theta.nonvalid[fold_indices[[fold]]$train]

cyclical_genes <- get.cyclical(Y=Y_train_normed_fold,
                               theta=theta_train_fold,
                               polyorder=2, ncores=ncores)


saveRDS(cyclical_genes,
        file=file.path(dir,
             paste0("data/results/triple_alge_ind_",ind,"_data_training_cyclical_genes.fold.",fold,".rds")))






