#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


fold <- args[1]

library(Biobase)
library(parallel)
library(genlasso)



#### Get data
df <- readRDS(file="/project2/gilad/joycehsiao/fucci-seq/data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

log2cpm.all <- log2cpm.all[,order(pdata$theta)]
pdata <- pdata[order(pdata$theta),]

log2cpm.quant <- readRDS("/project2/gilad/joycehsiao/fucci-seq/output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")


# select external validation samples
set.seed(99)
nvalid <- round(ncol(log2cpm.quant)*.15)
ii.valid <- sample(1:ncol(log2cpm.quant), nvalid, replace = F)
ii.nonvalid <- setdiff(1:ncol(log2cpm.quant), ii.valid)

log2cpm.quant.nonvalid <- log2cpm.quant[,ii.nonvalid]
log2cpm.quant.valid <- log2cpm.quant[,ii.valid]
theta.nonvalid <- pdata$theta[ii.nonvalid]
theta.valid <- pdata$theta[ii.valid]

sig.genes <- readRDS("/project2/gilad/joycehsiao/fucci-seq/output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.101.rds")
expr.sig <- log2cpm.quant.nonvalid[rownames(log2cpm.quant.nonvalid) %in% rownames(sig.genes), ]




# set training samples
source("/project2/gilad/joycehsiao/fucci-seq/peco/R/primes.R")
source("/project2/gilad/joycehsiao/fucci-seq/peco/R/partitionSamples.R")
parts <- partitionSamples(1:ncol(log2cpm.quant.nonvalid), runs=5,
                          nsize.each = rep(151,5))
part_indices <- parts$partitions


#### Fit model
source("/project2/gilad/joycehsiao/fucci-seq/peco/R/fit.trendfilter.generic.R")
source("/project2/gilad/joycehsiao/fucci-seq/peco/R/cycle.npreg.R")

#fold.train.pca <- vector("list", 5)
#fold.test.bytrain.pca <- vector("list", 5)
#fold.test.insample.pca <- vector("list", 5)
# fold.train <- vector("list", 5)
# fold.test.bytrain.fucci <- vector("list", 5)
# fold.test.insample.fucci <- vector("list", 5)
fold <- as.numeric(fold)

#for (i in 1:5) {
  Y_train <- expr.sig[,part_indices[[fold]]$train]
  theta_train <- theta.nonvalid[part_indices[[fold]]$train]

  Y_test <- expr.sig[,part_indices[[fold]]$test]
  theta_test <- theta.nonvalid[part_indices[[fold]]$test]

  fold.train <- cycle.npreg.insample(Y = Y_train,
                                     theta = theta_train,
                                     nbins = 100, ncores=12)

  # fold.train <- readRDS(paste0("/project2/gilad/joycehsiao/fucci-seq/output/method-npreg.Rmd/fold.train.",
  #                 "fold.",fold,".rds"))

  fold.test.bytrain.fucci <- cycle.npreg.outsample(Y_test,
                                                  theta_est=fold.train$cell_times_est,
                                                  mu_est=fold.train$mu_est,
                                                  sigma_est=fold.train$sigma_est)

  fold.test.insample.fucci <- cycle.npreg.insample(Y = Y_test,
                                                    theta = theta_test,
                                                    nbins = 100, ncores=12)



saveRDS(fold.train,
        paste0("/project2/gilad/joycehsiao/fucci-seq/output/method-npreg.Rmd/fold.train.",
               "fold.",fold,".rds"))
saveRDS(fold.test.bytrain.fucci,
        paste0("/project2/gilad/joycehsiao/fucci-seq/output/method-npreg.Rmd/fold.test.bytrain.fucci.",
               "fold.",fold,".rds"))
saveRDS(fold.test.insample.fucci,
        paste0("/project2/gilad/joycehsiao/fucci-seq/output/method-npreg.Rmd/fold.test.insample.fucci.",
               "fold.",fold,".rds"))

