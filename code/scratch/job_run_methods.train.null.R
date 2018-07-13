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

dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"code/working/run_methods.R"))
source(file.path(dir,"peco/R/fit.trendfilter.generic.R"))

data_training <- readRDS(file=file.path(dir, "data/results/data_training.rds"))
fold_indices <- readRDS(file=file.path(dir, "data/results/fold_indices.rds"))

nperm <- nrow(data_training$log2cpm.quant.nonvalid)
#nperm <- 2

print(fold)

# choose a gene with low fraction of missing and permute data
set.seed(31)
ii.frac.miss.low <- names(sample(which(
  rowMeans(data_training$log2cpm.nonvalid[,fold_indices[[fold]]$train]==0) < .1),1))

out <- mclapply(1:nperm, function(g) {
  fit.trendfilter.generic(
    sample(data_training$log2cpm.quant.nonvalid[ii.frac.miss.low,
    fold_indices[[fold]]$train]), polyorder = 2)$pve
}, mc.cores=ncores)

out <- do.call(c, out)

saveRDS(out,
        file=file.path(dir, paste0("data/results/results_train.fold.",fold,".null.rds")))





