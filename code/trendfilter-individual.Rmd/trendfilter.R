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


library(parallel)
library(genlasso)
library(Biobase)

dir <- "/project2/gilad/joycehsiao/fucci-seq"

df_sub <- readRDS(file.path(dir, paste0("output/trendfilter-individual.Rmd/df_",args[1],".rds")))
# log2cpm.quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))
# eset <- readRDS(file.path(dir,"data/eset-final.rds"))
# pdata <- pData(eset)
#
# # order data by cell time
# pdata_ordered <- pdata[order(pdata$theta),]
# df <- log2cpm.quant[,match(rownames(pdata_ordered),colnames(log2cpm.quant))]

# extract data of each individual
# df_NA19098 <- df[,pdata_ordered$chip_id == "NA19098"]
# df_NA18511 <- df[,pdata_ordered$chip_id == "NA18511"]
# df_NA18870 <- df[,pdata_ordered$chip_id == "NA18870"]
# df_NA19101 <- df[,pdata_ordered$chip_id == "NA19101"]
# df_NA18855 <- df[,pdata_ordered$chip_id == "NA18855"]
# df_NA19160 <- df[,pdata_ordered$chip_id == "NA19160"]
# saveRDS(df_NA19098,
#         file.path(dir, "output/trendfilter-individual.Rmd/df_NA19098.rds"))
# saveRDS(df_NA18511,
#         file.path(dir, "output/trendfilter-individual.Rmd/df_NA18511.rds"))
# saveRDS(df_NA18870,
#         file.path(dir, "output/trendfilter-individual.Rmd/df_NA18870.rds"))
# saveRDS(df_NA19101,
#         file.path(dir, "output/trendfilter-individual.Rmd/df_NA19101.rds"))
# saveRDS(df_NA18855,
#         file.path(dir, "output/trendfilter-individual.Rmd/df_NA18855.rds"))
# saveRDS(df_NA19160,
#         file.path(dir, "output/trendfilter-individual.Rmd/df_NA19160.rds"))

source(file.path(dir,"peco/R/fit.trendfilter.generic.R"))

# fit.trend <- lapply(1:6, function(i) {
#   id <- unique(pdata_ordered$chip_id)[i]
#   df_sub <- df[,which(pdata_ordered$chip_id == id)]
# #  fits <- mclapply(1:nrow(log2cpm.quant), function(g) {
#   fits <- mclapply(1:100, function(g) {
#   fit.trendfilter.generic(df_sub[g,], polyorder = 2)
#   }, mc.cores=15)
#   names(fits) <- rownames(df_sub)[1:100]
#   return(fits)
# })
# names(fit.trend) <- unique(pdata_ordered$chip_id)
#

fits <- mclapply(1:nrow(df_sub), function(g) {
  fit.trendfilter.generic(df_sub[g,], polyorder = 2)
}, mc.cores=15)
names(fits) <- rownames(df_sub)


saveRDS(fits,
        file.path(dir,paste0("output/trendfilter-individual.Rmd/fits_",args[1],".rds")))


