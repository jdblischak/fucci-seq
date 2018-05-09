#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
# To run this script on the command line
# sbatch

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


library(Biobase)
library(parallel)
library(genlasso)

dir <- "/project2/gilad/joycehsiao/fucci-seq"

df <- readRDS(file=file.path(dir,"data/eset-final.rds"))
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

log2cpm.all <- log2cpm.all[,order(pdata$theta)]
pdata <- pdata[order(pdata$theta),]

log2cpm.quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))

source(file.path(dir,"code/npreg/npreg.methods.R"))

nperm <- 1000

# choose a gene with high fraction of missing and permute data
set.seed(17)
ii.frac.miss.high <- names(sample(which(rowMeans(log2cpm.all==0) > .8),1))

fit.trend.highmiss <- mclapply(1:nperm, function(g) {
  fit.trendfilter.generic(sample(log2cpm.quant[ii.frac.miss.high,]),
                              polyorder = 2)
}, mc.cores=25)

saveRDS(fit.trend.highmiss,
        file.path(dir,"output/npreg-trendfilter-quantile.Rmd/fit.trend.perm.highmiss.rds"))


# choose a gene with low fraction of missing and permute data
set.seed(31)
ii.frac.miss.low <- names(sample(which(rowMeans(log2cpm.all==0) < .1),1))

fit.trend.lowmiss <- mclapply(1:nperm, function(g) {
  fit.trendfilter.generic(sample(log2cpm.quant[ii.frac.miss.low,]),
                  polyorder = 2)
}, mc.cores=25)

saveRDS(fit.trend.lowmiss,
        file.path(dir,"output/npreg-trendfilter-quantile.Rmd/fit.trend.perm.lowmiss.rds"))

