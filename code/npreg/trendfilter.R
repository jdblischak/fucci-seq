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

source(file.path(dir,"code/npreg/npreg.methods.R"))


ii.detected <- rowMeans(counts > 0) > .95
log2cpm.detected <- log2cpm.all[ii.detected,]
log2cpm.detected.sub <- log2cpm.detected[1:100,]

fit.trend <- mclapply(1:nrow(log2cpm.detected.sub), function(g) {
  fit.trendfilter.includezero(log2cpm.detected.sub[g,],
                              polyorder = 2)
}, mc.cores=10)

saveRDS(fit.trend, file = file.path(dir,"output/npreg-statistic-calibrate.Rmd/fit.trend.rds"))

