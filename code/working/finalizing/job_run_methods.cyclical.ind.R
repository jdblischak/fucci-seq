#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


ncores <- as.numeric(args[1])
ind <- as.character(args[2])

dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"peco/R/fit.trendfilter.generic.R"))
source(file.path(dir,"peco/R/get.cyclical.R"))


print(ind)

data_training <- readRDS(file=file.path(dir,
                         paste0("data/results/ind_",ind,"_data_training.rds")))

# fold_indices <- readRDS(file=file.path(dir,
#                          paste0("data/results/ind_", ind, "_fold_indices.rds")))

#print(fold)

Y_train_normed <- data_training$log2cpm.quant.nonvalid
theta_train <- data_training$theta.nonvalid

cyclical_genes <- get.cyclical(Y=Y_train_normed,
                               theta=theta_train,
                               polyorder=2, ncores=ncores)

saveRDS(cyclical_genes,
        file=file.path(dir,
             paste0("data/results/finalizing/ind_",ind,"_cyclical_genes.rds")))






